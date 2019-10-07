#!/usr/bin/env python3

''' ogr2osm beta

This program takes any vector data understadable by OGR and outputs an OSM file
with that data.

By default tags will be naively copied from the input data. Hooks are provided
so that, with a little python programming, you can translate the tags however
you like. More hooks are provided so you can filter or even modify the features
themselves.

To use the hooks, create a file in the translations/ directory called myfile.py
and run ogr2osm.py -t myfile. This file should define a function with the name
of each hook you want to use. For an example, see the uvmtrans.py file.

The program will use projection metadata from the source, if it has any. If
there is no projection information, or if you want to override it, you can use
-e or -p to specify an EPSG code or Proj.4 string, respectively. If there is no
projection metadata and you do not specify one, EPSG:4326 will be used (WGS84
latitude-longitude)

For additional usage information, run ogr2osm.py --help

Copyright (c) 2012-2013 Paul Norman <penorman@mac.com>, Sebastiaan Couwenberg
<sebastic@xs4all.nl>, The University of Vermont <andrew.guertin@uvm.edu>

Released under the MIT license: http://opensource.org/licenses/mit-license.php

Based very heavily on code released under the following terms:

(c) Iván Sánchez Ortega, 2009
<ivan@sanchezortega.es>
###############################################################################
#  "THE BEER-WARE LICENSE":                                                   #
#  <ivan@sanchezortega.es> wrote this file. As long as you retain this notice #
#  you can do whatever you want with this stuff. If we meet some day, and you #
#  think this stuff is worth it, you can buy me a beer in return.             #
###############################################################################

'''

import argparse
import logging as l
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Tuple, Union

from geom import Feature, Point, Relation, Way
from lxml import etree
from osgeo import ogr, osr

# By default, OGR doesn't raise exceptions, but we want to quit if something isn't working
# Thanks to OGR weirdness, exceptions are not raised on ogr.Open()
ogr.UseExceptions()

l.basicConfig(level=l.DEBUG, format="%(message)s")


class OSMSink:
    """
    Can accept an ogr.DataSource or a path to one in string or Path form
    """

    destspatialref = osr.SpatialReference()
    destspatialref.ImportFromEPSG(4326)

    def __init__(self, sourcepath: Union[str, Path, ogr.DataSource], *args, **kwargs):
        self.never_upload = kwargs.get('never_upload', False)
        self.no_upload_false = kwargs.get('no_upload_false', False)
        self.never_download = kwargs.get('never_download', False)
        self.locked = kwargs.get('locked', False)
        self.add_version = kwargs.get('add_version', False)
        self.add_timestamp = kwargs.get('add_timestamp', False)
        self.significant_digits = kwargs.get('significant_digits', 9)
        self.rounding_digits = kwargs.get('rounding_digits', 7)
        self.sqlquery = kwargs.get('sql_query', None)
        self.element_id_counter = kwargs.get('id', 0)
        self.force_overwrite = kwargs.get('force_overwrite', False)
        if kwargs.get("positive_id", False):
            self.element_id_counter_incr = 1
        else:
            self.element_id_counter_incr = -1

        nomemorycopy = kwargs.get('no_memory_copy', None)
        self.geometries = []
        self.features = []
        self.long_ways_from_polygons = set()
        # Running dict of points created. If a point already has been created,
        # it is referenced from the dict rather than created again.
        self.linestring_points = {}

        # Open OGR source for reading
        # If passed an ogr datasource directly:
        if isinstance(sourcepath, ogr.DataSource):
            self.source = sourcepath
        elif isinstance(sourcepath, str) and re.match('^PG:', sourcepath):
            self.source = ogr.Open(sourcepath, 0)  # 0 means read-only
        else:
            sourcepath = Path(sourcepath)
            self.source = self.get_file_data(sourcepath, nomemorycopy)
        if self.source is None:
            raise RuntimeError

        # Projection
        proj4 = kwargs.get('source_proj4', None)
        sourceepsg = kwargs.get('source_epsg', None)
        if proj4:
            spatialref = osr.SpatialReference().ImportFromProj4(proj4)
        elif sourceepsg and sourceepsg != 4326:
            spatialref = osr.SpatialReference().ImportFromEPSG(sourceepsg)
        else:
            spatialref = None
        if spatialref:
            self.coordtrans = osr.CoordinateTransformation(
                spatialref, self.destspatialref)
        else:
            self.coordtrans = None

    def get_new_id(self) -> int:
        self.element_id_counter += self.element_id_counter_incr
        return self.element_id_counter

    @staticmethod
    def get_file_data(filename: Path, nomemorycopy: bool = None) -> ogr.DataSource:
        ogr_accessmethods = {"vsicurl",
                             "vsicurl_streaming", "vsisubfile", "vsistdin"}
        ogr_filemethods = {"vsisparse", "vsigzip", "vsitar", "vsizip"}
        ogr_unsupported = {"vsimem", "vsistdout"}
        has_unsup = any([m in filename.parts for m in ogr_unsupported])
        if has_unsup:
            l.error("Unsupported OGR access method(s) found: %s." %
                    str(has_unsup))
            raise RuntimeError
        if not any([m in filename.parts for m in ogr_accessmethods]):
            # Not using any ogr_accessmethods
            if not any([m in filename.parts for m in ogr_filemethods]):
                if filename.suffix == '.gz':
                    filename = Path(f'/vsigzip/{str(filename)}')
                elif any([ext in filename.suffixes for ext in {'.tar', '.tgz'}]):
                    filename = Path(f'/vsitar/{str(filename)}')
                elif filename.suffix == '.zip':
                    filename = Path(f'/vsizip/{str(filename)}')
        file_data_source = ogr.Open(str(filename), 0)  # 0 means read-only
        if not nomemorycopy:
            file_data_source = ogr.GetDriverByName(
                'Memory').CopyDataSource(file_data_source, 'memoryCopy')
        return file_data_source

    def parse_data(self):
        l.debug("Parsing data")
        if self.sqlquery:
            layer = self.source.ExecuteSQL(self.sqlquery)
            layer.ResetReading()
            self.parse_layer(self.translations.filter_layer(layer))
        else:
            for i in range(self.source.GetLayerCount()):
                layer = self.source.GetLayer(i)
                layer.ResetReading()
                self.parse_layer(
                    self.translations.filter_layer(layer))

    def parse_layer(self, layer: ogr.Layer):
        if layer is None:
            return
        field_names = self.get_layer_fields(layer)
        if self.coordtrans:
            layer_coordtrans = self.coordtrans
        else:
            spatialref = layer.GetSpatialRef()
            # If source is already in WGS84, don't bother reprojecting
            if spatialref != self.destspatialref:
                layer_coordtrans = osr.CoordinateTransformation(
                    spatialref, self.destspatialref)
            else:
                layer_coordtrans = None

        for j in range(layer.GetFeatureCount()):
            ogrfeature = layer.GetNextFeature()
            self.parse_feature(self.translations.filter_feature(
                ogrfeature, field_names, layer_coordtrans), field_names, layer_coordtrans)

    def parse_feature(self, ogrfeature: ogr.Feature, field_names: list,
                      coordtrans: osr.CoordinateTransformation):
        if ogrfeature is None:
            return
        ogrgeometry = ogrfeature.GetGeometryRef()
        if ogrgeometry is None:
            return
        if coordtrans:
            ogrgeometry.Transform(coordtrans)
        geometries = self.parse_geometry([ogrgeometry])
        for geometry in geometries:
            if geometry is None:
                return
            feature = Feature(self)
            feature.tags = self.get_feature_tags(ogrfeature, field_names)
            feature.geometry = geometry
            geometry.addparent(feature)

            self.translations.filter_feature_post(
                feature, ogrfeature, ogrgeometry)

    def parse_geometry(self, ogrgeometries) -> list:
        returngeometries = []
        for ogrgeometry in ogrgeometries:
            geometryType = ogrgeometry.GetGeometryType()

            if (geometryType == ogr.wkbPoint or
                    geometryType == ogr.wkbPoint25D):
                returngeometries.append(self.parse_point(ogrgeometry))
            elif (geometryType == ogr.wkbLineString or
                  geometryType == ogr.wkbLinearRing or
                  geometryType == ogr.wkbLineString25D):
                # geometryType == ogr.wkbLinearRing25D does not exist
                returngeometries.append(self.parse_line_string(ogrgeometry))
            elif (geometryType == ogr.wkbPolygon or
                  geometryType == ogr.wkbPolygon25D):
                returngeometries.append(self.parse_polygon(ogrgeometry))
            elif (geometryType == ogr.wkbMultiPoint or
                  geometryType == ogr.wkbMultiLineString or
                  geometryType == ogr.wkbMultiPolygon or
                  geometryType == ogr.wkbGeometryCollection or
                  geometryType == ogr.wkbMultiPoint25D or
                  geometryType == ogr.wkbMultiLineString25D or
                  geometryType == ogr.wkbMultiPolygon25D or
                  geometryType == ogr.wkbGeometryCollection25D):
                returngeometries.extend(self.parse_collection(ogrgeometry))
            else:
                l.warning("unhandled geometry, type: " + str(geometryType))
                returngeometries.append(None)

        return returngeometries

    @staticmethod
    def get_layer_fields(layer: ogr.Layer) -> list:
        featureDefinition = layer.GetLayerDefn()
        field_names = []
        fieldCount = featureDefinition.GetFieldCount()
        for j in range(fieldCount):
            field_names.append(featureDefinition.GetFieldDefn(j).GetNameRef())
        return field_names

    def parse_point(self, ogrgeometry) -> Point:
        x = int(round(ogrgeometry.GetX() * 10**self.significant_digits))
        y = int(round(ogrgeometry.GetY() * 10**self.significant_digits))
        geometry = Point(self, x, y)
        return geometry

    def parse_line_string(self, ogrgeometry) -> Way:
        geometry = Way(self)
        # LineString.GetPoint() returns a tuple, so we can't call parse_point on it
        # and instead have to create the point ourself
        for i in range(ogrgeometry.GetPointCount()):
            (x, y, _) = ogrgeometry.GetPoint(i)
            (rx, ry) = (int(round(x*10**self.rounding_digits)),
                        int(round(y*10**self.rounding_digits)))
            (x, y) = (int(round(x*10**self.significant_digits)),
                      int(round(y*10**self.significant_digits)))
            if (rx, ry) in self.linestring_points:
                mypoint = self.linestring_points[(rx, ry)]
            else:
                mypoint = Point(self, x, y)
                self.linestring_points[(rx, ry)] = mypoint
            geometry.points.append(mypoint)
            mypoint.addparent(geometry)
        return geometry

    def parse_polygon(self, ogrgeometry):
        # Special case polygons with only one ring. This does not (or at least
        # should not) change behavior when simplify relations is turned on.
        if ogrgeometry.GetGeometryCount() == 0:
            l.warning("Polygon with no rings?")
        elif ogrgeometry.GetGeometryCount() == 1:
            result = self.parse_line_string(ogrgeometry.GetGeometryRef(0))
            if len(result.points) > self.max_nodes_per_way:
                self.long_ways_from_polygons.add(result)
            return result
        else:
            geometry = Relation(self)
            try:
                exterior = self.parse_line_string(
                    ogrgeometry.GetGeometryRef(0))
                exterior.addparent(geometry)
            except:
                l.warning("Polygon with no exterior ring?")
                return None
            geometry.members.append((exterior, "outer"))
            for i in range(1, ogrgeometry.GetGeometryCount()):
                interior = self.parse_line_string(
                    ogrgeometry.GetGeometryRef(i))
                interior.addparent(geometry)
                geometry.members.append((interior, "inner"))
            return geometry

    def parse_collection(self, ogrgeometry):
        # OGR MultiPolygon maps easily to osm multipolygon, so special case it
        # TODO: Does anything else need special casing?
        geometryType = ogrgeometry.GetGeometryType()
        if (geometryType == ogr.wkbMultiPolygon or
                geometryType == ogr.wkbMultiPolygon25D):
            if ogrgeometry.GetGeometryCount() > 1:
                geometry = Relation(self)
                for polygon in range(ogrgeometry.GetGeometryCount()):
                    exterior = self.parse_line_string(
                        ogrgeometry.GetGeometryRef(polygon).GetGeometryRef(0))
                    exterior.addparent(geometry)
                    geometry.members.append((exterior, "outer"))
                    for i in range(1, ogrgeometry.GetGeometryRef(polygon).GetGeometryCount()):
                        interior = self.parse_line_string(
                            ogrgeometry.GetGeometryRef(polygon).GetGeometryRef(i))
                        interior.addparent(geometry)
                        geometry.members.append((interior, "inner"))
                return [geometry]
            else:
                return [self.parse_polygon(ogrgeometry.GetGeometryRef(0))]
        elif (geometryType == ogr.wkbMultiLineString or
              geometryType == ogr.wkbMultiLineString25D):
            geometries = []
            for linestring in range(ogrgeometry.GetGeometryCount()):
                geometries.append(self.parse_line_string(
                    ogrgeometry.GetGeometryRef(linestring)))
            return geometries
        else:
            geometry = Relation()
            for i in range(ogrgeometry.GetGeometryCount()):
                member = self.parse_geometry(ogrgeometry.GetGeometryRef(i))
                member.addparent(geometry)
                geometry.members.append((member, "member"))
            return [geometry]

    def get_feature_tags(self, ogrfeature: ogr.Feature, field_names: list) -> dict:
        '''
        This function builds up a dictionary with the source data attributes
        and passes them to the filter_tags function, returning the result.
        '''
        tags = {}
        for index, item in enumerate(field_names):
            tags[item] = ogrfeature.GetFieldAsString(index).strip()
        return self.translations.filter_tags(tags)
        # return tags

    def merge_points(self):
        l.debug("Merging points")
        points = [geom for geom in self.geometries if type(geom) == Point]

        # Make list of Points at each location
        l.debug("Making list")
        pointcoords = {}
        for i in points:
            rx = int(
                round(i.x * 10**(self.significant_digits-self.rounding_digits)))
            ry = int(
                round(i.y * 10**(self.significant_digits-self.rounding_digits)))
            if (rx, ry) in pointcoords:
                pointcoords[(rx, ry)].append(i)
            else:
                pointcoords[(rx, ry)] = [i]

        # Use list to get rid of extras
        l.debug("Checking list")
        for (location, pointsatloc) in pointcoords.items():
            if len(pointsatloc) > 1:
                for point in pointsatloc[1:]:
                    for parent in set(point.parents):
                        parent.replacejwithi(pointsatloc[0], point)

    def merge_way_points(self):
        l.debug("Merging duplicate points in ways")
        ways = [geom for geom in self.geometries if type(geom) == Way]

        # Remove duplicate points from ways,
        # a duplicate has the same id as its predecessor
        for way in ways:
            previous = self.element_id_counter
            merged_points = []

            for node in way.points:
                if previous == self.element_id_counter or previous != node.id:
                    merged_points.append(node)
                    previous = node.id

            if len(merged_points) > 0:
                way.points = merged_points

    def split_long_ways(self, max_points_in_way, waysToCreateRelationFor):
        l.debug("Splitting long ways")
        ways = [geom for geom in self.geometries if type(geom) == Way]

        featuresmap = {
            feature.geometry: feature for feature in self.features}

        for way in ways:
            is_way_in_relation = 0 < len(
                [p for p in way.parents if type(p) == Relation])
            if len(way.points) > max_points_in_way:
                way_parts = self.split_way(way, max_points_in_way,
                                           featuresmap, is_way_in_relation)
                if not is_way_in_relation:
                    if way in waysToCreateRelationFor:
                        self.merge_into_new_relation(way_parts)
                else:
                    for rel in way.parents:
                        self.split_way_in_relation(rel, way_parts)

    def split_way(self, way, max_points_in_way, features_map, is_way_in_relation: bool):
        new_points = [way.points[i:i + max_points_in_way]
                      for i in range(0, len(way.points), max_points_in_way - 1)]
        new_ways = [way, ] + [Way(self) for i in range(len(new_points) - 1)]

        if not is_way_in_relation:
            way_tags = features_map[way].tags

            for new_way in new_ways:
                if new_way != way:
                    feat = Feature(self)
                    feat.geometry = new_way
                    feat.tags = way_tags

        for new_way, points in zip(new_ways, new_points):
            new_way.points = points
            if new_way.id != way.id:
                for point in points:
                    point.removeparent(self, way, shoulddestroy=False)
                    point.addparent(new_way)

        return new_ways

    def merge_into_new_relation(self, way_parts):
        new_relation = Relation(self)
        feat = Feature()
        feat.geometry = new_relation
        new_relation.members = [(way, "outer") for way in way_parts]
        for way in way_parts:
            way.addparent(new_relation)

    @staticmethod
    def split_way_in_relation(rel, way_parts):
        way_roles = [m[1] for m in rel.members if m[0] == way_parts[0]]
        way_role = "" if len(way_roles) == 0 else way_roles[0]
        for way in way_parts[1:]:
            rel.members.append((way, way_role))

    def output(self, outputfile: Union[str, Path],  **kwargs):
        """
        Outputs the sink to the specified path

        :param outputfile: string or Path representing the output location
        :param force_overwrite: if True, overwrite a file if it exists
        """
        force_overwrite = kwargs.get('force_overwrite', self.force_overwrite)
        never_upload = kwargs.get('never_upload', self.never_upload)
        no_upload_false = kwargs.get('no_upload_false', self.no_upload_false)
        never_download = kwargs.get('never_download', self.never_download)
        locked = kwargs.get('locked', self.locked)
        add_version = kwargs.get('add_version', self.add_version)
        add_timestamp = kwargs.get('add_timestamp', self.add_timestamp)
        significant_digits = kwargs.get(
            'significant_digits', self.significant_digits)

        # Promote string to Path if neccessary, has no effect if already a Path
        outputfile = Path(outputfile)

        if force_overwrite:
            write_mode = 'w'
        else:
            write_mode = 'x'

        l.debug("Outputting XML")
        # First, set up a few data structures for optimization purposes
        nodes = [geom for geom in self.geometries if type(geom) == Point]
        ways = [geom for geom in self.geometries if type(geom) == Way]
        relations = [
            geom for geom in self.geometries if type(geom) == Relation]
        featuresmap = {
            feature.geometry: feature for feature in self.features}

        # Open up the output file with the system default buffering
        with outputfile.open(write_mode, buffering=-1) as f:

            dec_string = '<?xml version="1.0"?>\n<osm version="0.6" generator="uvmogr2osm"'
            if never_upload:
                dec_string += ' upload="never"'
            elif not no_upload_false:
                dec_string += ' upload="false"'
            if never_download:
                dec_string += ' download="never"'
            if locked:
                dec_string += ' locked="true"'
            dec_string += '>\n'
            f.write(dec_string)

            # Build up a dict for optional settings
            attributes = {}
            if add_version:
                attributes.update({'version': '1'})

            if add_timestamp:
                attributes.update(
                    {'timestamp': datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')})

            for node in nodes:
                xmlattrs = {'visible': 'true', 'id': str(node.id), 'lat': str(
                    node.y*10**-significant_digits), 'lon': str(node.x*10**-significant_digits)}
                xmlattrs.update(attributes)

                xmlobject = etree.Element('node', xmlattrs)

                if node in featuresmap:
                    for (key, value) in featuresmap[node].tags.items():
                        tag = etree.Element('tag', {'k': key, 'v': value})
                        xmlobject.append(tag)
                f.write(etree.tostring(xmlobject, encoding='unicode'))
                f.write('\n')

            for way in ways:
                xmlattrs = {'visible': 'true', 'id': str(way.id)}
                xmlattrs.update(attributes)

                xmlobject = etree.Element('way', xmlattrs)

                for node in way.points:
                    nd = etree.Element('nd', {'ref': str(node.id)})
                    xmlobject.append(nd)
                if way in featuresmap:
                    for (key, value) in featuresmap[way].tags.items():
                        tag = etree.Element('tag', {'k': key, 'v': value})
                        xmlobject.append(tag)
                f.write(etree.tostring(xmlobject, encoding='unicode'))
                f.write('\n')

            for relation in relations:
                xmlattrs = {'visible': 'true', 'id': str(relation.id)}
                xmlattrs.update(attributes)

                xmlobject = etree.Element('relation', xmlattrs)

                for (member, role) in relation.members:
                    member = etree.Element(
                        'member', {'type': 'way', 'ref': str(member.id), 'role': role})
                    xmlobject.append(member)

                tag = etree.Element('tag', {'k': 'type', 'v': 'multipolygon'})
                xmlobject.append(tag)
                if relation in featuresmap:
                    for (key, value) in featuresmap[relation].tags.items():
                        tag = etree.Element('tag', {'k': key, 'v': value})
                        xmlobject.append(tag)

                f.write(etree.tostring(xmlobject, encoding='unicode'))
                f.write('\n')

            f.write('</osm>')


def setup(args: dict) -> Tuple[dict, argparse.ArgumentParser]:
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    parser.add_argument('source',
                        type=str, metavar="INPUT",
                        help="OGR datasource for the data")
    parser.add_argument("-t", "--translation", dest="translation_method",
                        metavar="TRANSLATION",
                        help="Select the attribute-tags translation method. See " +
                        "the translations/ directory for valid values.")
    parser.add_argument("-o", "--output", dest="output_file", metavar="OUTPUT",
                        help="Set destination .osm file name and location.")
    parser.add_argument("-e", "--epsg", dest="source_epsg", metavar="EPSG_CODE", type=int,
                        help="EPSG code of source file. Do not include the " +
                        "'EPSG:' prefix. If specified, overrides projection " +
                        "from source metadata if it exists.")
    parser.add_argument("-p", "--proj4", dest="source_proj4", metavar="PROJ4_STRING",
                        help="PROJ.4 string. If specified, overrides projection " +
                        "from source metadata if it exists.")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true")
    parser.add_argument("-d", "--debug-tags", dest="debug_tags", action="store_true",
                        help="Output the tags for every feature parsed.")
    parser.add_argument("-f", "--force", dest="force_overwrite", action="store_true",
                        help="Force overwrite of output file.")
    parser.add_argument("--significant-digits",  dest="significant_digits", type=int,
                        help="Number of decimal places for coordinates", default=9)
    parser.add_argument("--rounding-digits",  dest="rounding_digits", type=int,
                        help="Number of decimal places for rounding", default=7)
    parser.add_argument("--no-memory-copy", dest="no_memory_copy", action="store_true",
                        help="Do not make an in-memory working copy")
    parser.add_argument("--no-upload-false", dest="no_upload_false", action="store_true",
                        help="Omit upload=false from the completed file to " +
                        "supress JOSM warnings when uploading.")
    parser.add_argument("--never-download", dest="never_download", action="store_true",
                        help="Prevent JOSM from downloading more data to this file.")
    parser.add_argument("--never-upload", dest="never_upload", action="store_true",
                        help="Completely disables all upload commands for this file in JOSM, " +
                        "rather than merely showing a warning before uploading.")
    parser.add_argument("--locked", dest="locked", action="store_true",
                        help="Prevent any changes to this file in JOSM, " +
                        "such as editing or downloading, and also prevents uploads. " +
                        "Implies upload=\"never\" and download=\"never\".")
    parser.add_argument("--id", dest="id", type=int, default=0,
                        help="ID to start counting from for the output file. Defaults to 0.")
    parser.add_argument("--idfile", dest="idfile", type=str, default=None,
                        help="Read ID to start counting from from a file.")
    parser.add_argument("--split-ways", dest="max_nodes_per_way", type=int, default=1800,
                        help="Split ways with more than the specified number of nodes. " +
                        "Defaults to 1800. Any value below 2 - do not split.")
    parser.add_argument("--saveid", dest="saveid", type=str, default=None,
                        help="Save last ID after execution to a file.")
    # Positive IDs can cause big problems if used inappropriately so hide the help for this
    parser.add_argument("--positive-id", dest="positive_id", action="store_true",
                        help=argparse.SUPPRESS)
    # Add version attributes. Again, this can cause big problems so surpress the help
    parser.add_argument("--add-version", dest="add_version", action="store_true",
                        help=argparse.SUPPRESS)
    # Add timestamp attributes. Again, this can cause big problems so surpress the help
    parser.add_argument("--add-timestamp", dest="add_timestamp", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--sql", dest="sql_query", type=str, default=None,
                        help="SQL query to execute on a PostgreSQL source")

    # parser.set_defaults(source_epsg=None, source_proj4=None, verbose=False,
    #                     debug_tags=False,
    #                     translation_method=None, output_file=None,
    #                     force_overwrite=False, no_upload_false=False,
    #                     never_download=False, never_upload=False,
    #                     locked=False)

    # Parse and process arguments
    options = vars(parser.parse_args(args))

    if options["source_epsg"]:
        try:
            options["source_epsg"] = int(options["source_epsg"])
        except ValueError:
            parser.error(
                "EPSG code must be numeric (e.g. '4326', not 'epsg:4326')")

    if options["locked"] and (options["never_upload"] or options["never_download"]):
        l.info("When specifying a locked file, it is not neccessary to also "
               "specify the never download/never upload options.")

    source_is_database = bool(re.match('^PG:', options["source"]))

    if options["output_file"]:
        options["output_file"] = Path(options["output_file"]).resolve()
    elif source_is_database:
        parser.error(
            "ERROR: An output file must be explicitly specified when using a database source")
    else:
        # If no output specified, use input name, in working directory, with .osm extension
        newname = Path(options["source"]).with_suffix(".osm").name
        options["output_file"] = Path.cwd() / newname

    if options["sql_query"] and not source_is_database:
        parser.error(
            "ERROR: You must use a database source when specifying a query with --sql")

    if not options["force_overwrite"] and options["output_file"].exists():
        parser.error("ERROR: output file '%s' exists" %
                     (options["output_file"]))

    l.info("Preparing to convert '%s' to '%s'." %
           (options["source"], options["output_file"]))

    # Projection
    if not options["source_proj4"] and not options["source_epsg"]:
        l.info(
            "Will try to detect projection from source metadata, or fall back to EPSG:4326")
    elif options["source_proj4"]:
        l.info("Will use the PROJ.4 string: " + options["source_proj4"])
    elif options["source_epsg"]:
        l.info("Will use EPSG:" + str(options["source_epsg"]))

    if options["idfile"]:
        try:
            with open(options["idfile"], 'r') as ff:
                options["id"] = int(ff.readline(20))
        except OSError:
            l.exception("Couldn't read id file %s." % options["idfile"])
        else:
            l.info("Starting counter value '%d' read from file '%s'."
                   % (options["id"], options["idfile"]))

    return options, parser


def main(args: dict):
    # Parse args and return as dict, along with parser for errors
    options, parser = setup(args)

    # Create memory object for data destination
    try:
        sink = OSMSink(options["source"], **options)
    except RuntimeError:
        # TODO: Useful CLI error message here
        parser.error("Could not parse OGR data source.")

    # TODO: Move this into OSMSink?
    # Stuff needed for locating translation methods
    if options["translation_method"]:
        # add dirs to path if necessary
        path = Path(options["translation_method"])
        # (root, ext) = os.path.splitext(options["translation_method"])
        if path.exists() and path.suffix == ".py":
            # user supplied translation file directly
            # sys.path.insert(0, os.path.dirname(root))
            sys.path.insert(0, path.parent)
        else:
            # first check translations in the subdir translations of cwd
            sys.path.insert(0, Path.getcwd() / "translations")
            # then check subdir of script dir
            sys.path.insert(1, Path(__file__).parent /
                            "translations").resolve()
            # (the cwd will also be checked implicityly)

        # strip .py if present, as import wants just the module name
        if path.suffix == '.py':
            options["translation_method"] = path.name

        try:
            sink.translations = __import__(
                options["translation_method"], fromlist=[''])
        except ImportError:
            parser.error("Could not load translation method '%s'. Translation "
                         "script must be in your current directory, or in the "
                         "translations/ subdirectory of your current or ogr2osm.py "
                         "directory. The following directories have been considered: %s"
                         % (options["translation_method"], str(sys.path)))
        except SyntaxError as e:
            parser.error("Syntax error in '%s'. Translation script is malformed:\n%s"
                         % (options["translation_method"], e))

        l.info("Successfully loaded '%s' translation method ('%s')."
               % (options["translation_method"], Path(sink.translations.__file__).resolve()))
    else:
        import types
        sink.translations = types.ModuleType("translationmodule")
        l.info("Using default translations")

    default_translations = {
        'filter_layer': lambda layer: layer,
        'filter_feature': lambda feature, field_names, reproject: feature,
        'filter_tags': lambda tags: tags,
        'filter_feature_post': lambda feature, field_names, reproject: feature,
        'pre_output_transform': lambda geometries, features: None,
    }

    for k, v in default_translations.items():
        if hasattr(sink.translations, k) and getattr(sink.translations, k):
            l.debug("Using user " + k)
        else:
            l.debug("Using default " + k)
            setattr(sink.translations, k, v)

    # Main flow
    sink.parse_data()
    sink.merge_points()
    sink.merge_way_points()
    if options["max_nodes_per_way"] >= 2:
        sink.split_long_ways(options["max_nodes_per_way"],
                             sink.long_ways_from_polygons)
    sink.translations.pre_output_transform(sink.geometries, sink.features)
    # Save data to file
    sink.output(
        options["output_file"]
    )

    # Save last used id to file
    if options["saveid"]:
        with open(options["saveid"], 'wb') as ff:
            ff.write(str(sink.element_id_counter))
        l.info("Wrote element_id_counter '%d' to file '%s'"
               % (sink.element_id_counter, options["saveid"]))


if __name__ == '__main__':
    main(sys.argv[1:])

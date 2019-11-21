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
from typing import Union

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

    def __init__(self, sourcepath: Union[str, Path, ogr.DataSource], **kwargs):
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
        self.translation_method = kwargs.get('translation_method', None)
        self.max_nodes_per_way = kwargs.get('max_nodes_per_way', 1800)

        self.saveid = kwargs.get('saveid', None)
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
        if proj4 or (sourceepsg and sourceepsg != 4326):
            spatialref = osr.SpatialReference()
            if proj4:
                spatialref.ImportFromProj4(proj4)
            else:
                spatialref.ImportFromEPSG(sourceepsg)
            self.coordtrans = osr.CoordinateTransformation(
                spatialref, self.destspatialref)
        else:
            self.coordtrans = None

        self.translation_method_setup()

        # Main flow
        self.parse_data()
        self.merge_points()
        self.merge_way_points()
        if self.max_nodes_per_way >= 2:
            self.split_long_ways(self.max_nodes_per_way,
                                 self.long_ways_from_polygons)
        self.translations.pre_output_transform(self.geometries, self.features)

    def get_new_id(self) -> int:
        self.element_id_counter += self.element_id_counter_incr
        return self.element_id_counter

    def translation_method_setup(self):
        # Stuff needed for locating translation methods
        if self.translation_method:
            # add dirs to path if necessary
            path = Path(self.translation_method)
            # (root, ext) = os.path.splitext(self.translation_method)
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
                self.translation_method = path.name

            try:
                self.translations = __import__(
                    self.translation_method, fromlist=[''])
            except ImportError:
                # parser.error("Could not load translation method '%s'. Translation "
                #             "script must be in your current directory, or in the "
                #             "translations/ subdirectory of your current or ogr2osm.py "
                #             "directory. The following directories have been considered: %s"
                #             % (self.translation_method, str(sys.path)))
                pass
            except SyntaxError:
                # parser.error("Syntax error in '%s'. Translation script is malformed:\n%s"
                #             % (self.translation_method, e))
                pass

            l.info("Successfully loaded '%s' translation method ('%s')."
                   % (self.translation_method, Path(self.translations.__file__).resolve()))
        else:
            import types
            self.translations = types.ModuleType("translationmodule")
            l.info("Using default translations")

        default_translations = {
            'filter_layer': lambda layer: layer,
            'filter_feature': lambda feature, field_names, reproject: feature,
            'filter_tags': lambda tags: tags,
            'filter_feature_post': lambda feature, field_names, reproject: feature,
            'pre_output_transform': lambda geometries, features: None,
        }

        for k, v in default_translations.items():
            if hasattr(self.translations, k) and getattr(self.translations, k):
                l.debug("Using user " + k)
            else:
                l.debug("Using default " + k)
                setattr(self.translations, k, v)

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
                    # Can't use slash operator here because filename has a root element
                    # (the first slash) which prevents prepending to the Path
                    # Converting to str and then assembling is easier
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
            if ogrgeometry.HasCurveGeometry():  # OSM can't have curved geometry
                l.info('Converting %s type to linear type',
                       ogrgeometry.GetGeometryName())
                ogrgeometry = ogrgeometry.GetLinearGeometry()

            geometry_type = ogrgeometry.GetGeometryType()

            if geometry_type in {ogr.wkbPoint, ogr.wkbPoint25D}:
                returngeometries.append(self.parse_point(ogrgeometry))
            elif geometry_type in {ogr.wkbLineString, ogr.wkbLinearRing,
                                   # ogr.wkbLinearRing25D does not exist
                                   ogr.wkbLineString25D}:
                returngeometries.append(self.parse_line_string(ogrgeometry))
            elif geometry_type in {ogr.wkbPolygon, ogr.wkbPolygon25D}:
                returngeometries.append(self.parse_polygon(ogrgeometry))
            elif geometry_type in {ogr.wkbMultiPoint, ogr.wkbMultiLineString,
                                   ogr.wkbMultiPolygon, ogr.wkbGeometryCollection,
                                   ogr.wkbMultiPoint25D, ogr.wkbMultiLineString25D,
                                   ogr.wkbMultiPolygon25D, ogr.wkbGeometryCollection25D}:
                returngeometries.extend(self.parse_collection(ogrgeometry))
            else:
                l.warning("unhandled geometry, type: %s",
                          ogrgeometry.GetGeometryName())
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

    def parse_point(self, ogrgeometry: ogr.Geometry) -> Point:
        x = int(round(ogrgeometry.GetX() * 10**self.significant_digits))
        y = int(round(ogrgeometry.GetY() * 10**self.significant_digits))
        geometry = Point(self, x, y)
        return geometry

    def parse_line_string(self, ogrgeometry: ogr.Geometry) -> Way:
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

    def parse_polygon(self, ogrgeometry: ogr.Geometry) -> Union[Way, Relation]:
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

    def parse_collection(self, ogrgeometry: ogr.Geometry):
        # OGR MultiPolygon maps easily to osm multipolygon, so special case it
        # TODO: Does anything else need special casing?
        geometry_type = ogrgeometry.GetGeometryType()
        if geometry_type in {ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D}:
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
        elif geometry_type in {ogr.wkbMultiLineString, ogr.wkbMultiLineString25D}:
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
        points = [geom for geom in self.geometries if isinstance(geom, Point)]

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
        ways = [geom for geom in self.geometries if isinstance(geom,  Way)]

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

    def split_long_ways(self, max_points_in_way: int, waysToCreateRelationFor: set):
        l.debug("Splitting long ways")
        ways = [geom for geom in self.geometries if isinstance(geom, Way)]

        featuresmap = {
            feature.geometry: feature for feature in self.features}

        for way in ways:
            is_way_in_relation = 0 < len(
                [p for p in way.parents if isinstance(p, Relation)])
            if len(way.points) > max_points_in_way:
                way_parts = self.split_way(way, max_points_in_way,
                                           featuresmap, is_way_in_relation)
                if not is_way_in_relation:
                    if way in waysToCreateRelationFor:
                        self.merge_into_new_relation(way_parts)
                else:
                    for rel in way.parents:
                        self.split_way_in_relation(rel, way_parts)

    def split_way(self, way, max_points_in_way: int, features_map: dict, is_way_in_relation: bool):
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
                    point.removeparent(self, way, should_destroy=False)
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
        # The options set at instance creation can be overriden with kwargs
        force_overwrite = kwargs.get('force_overwrite', self.force_overwrite)
        never_upload = kwargs.get('never_upload', self.never_upload)
        no_upload_false = kwargs.get('no_upload_false', self.no_upload_false)
        never_download = kwargs.get('never_download', self.never_download)
        locked = kwargs.get('locked', self.locked)
        add_version = kwargs.get('add_version', self.add_version)
        add_timestamp = kwargs.get('add_timestamp', self.add_timestamp)
        significant_digits = kwargs.get(
            'significant_digits', self.significant_digits)
        saveid = kwargs.get('saveid', self.saveid)
        # This does nothing at the moment, but could be used in the future
        pretty_print = kwargs.get('pretty_print', False)

        # Promote string to Path if neccessary, has no effect if already a Path
        outputfile = Path(outputfile)

        if force_overwrite:
            write_mode = 'wb'
        else:
            write_mode = 'xb'

        l.debug("Outputting XML")
        # First, set up a few data structures for optimization purposes
        nodes = [geom for geom in self.geometries if isinstance(geom, Point)]
        ways = [geom for geom in self.geometries if isinstance(geom, Way)]
        relations = [
            geom for geom in self.geometries if isinstance(geom, Relation)]
        featuresmap = {
            feature.geometry: feature for feature in self.features}

        # Open up the output file with the system default buffering
        with outputfile.open(write_mode, buffering=-1) as f:
            with etree.xmlfile(f, encoding='utf-8', buffered=False) as xf:
                xf.write_declaration()
                root_attrib = {}
                if never_upload:
                    root_attrib["upload"] = "never"
                elif not no_upload_false:
                    root_attrib["upload"] = "false"
                if never_download:
                    root_attrib["download"] = "never"
                if locked:
                    root_attrib["locked"] = "true"
                with xf.element('osm', version="0.6", generator="uvmogr2osm", **root_attrib):
                    # Build up a dict for optional settings
                    attributes = {}
                    if add_version:
                        attributes['version'] = '1'

                    if add_timestamp:
                        attributes['timestamp'] = datetime.utcnow().strftime(
                            '%Y-%m-%dT%H:%M:%SZ')

                    for node in nodes:
                        xmlattrs = {'visible': 'true', 'id': str(node.id), 'lat': str(
                            node.y*10**-significant_digits), 'lon': str(node.x*10**-significant_digits)}
                        xmlattrs.update(attributes)

                        xmlobject = etree.Element('node', xmlattrs)

                        if node in featuresmap:
                            for (key, value) in sorted(featuresmap[node].tags.items()):
                                tag = etree.Element(
                                    'tag', {'k': key, 'v': value})
                                xmlobject.append(tag)
                        xf.write(xmlobject)

                    for way in ways:
                        xmlattrs = {'visible': 'true', 'id': str(way.id)}
                        xmlattrs.update(attributes)

                        xmlobject = etree.Element('way', xmlattrs)

                        for node in way.points:
                            nd = etree.Element('nd', {'ref': str(node.id)})
                            xmlobject.append(nd)
                        if way in featuresmap:
                            for (key, value) in sorted(featuresmap[way].tags.items()):
                                tag = etree.Element(
                                    'tag', {'k': key, 'v': value})
                                xmlobject.append(tag)
                        xf.write(xmlobject)

                    for relation in relations:
                        xmlattrs = {'visible': 'true', 'id': str(relation.id)}
                        xmlattrs.update(attributes)

                        xmlobject = etree.Element('relation', xmlattrs)

                        for (member, role) in relation.members:
                            member = etree.Element(
                                'member', {'type': 'way', 'ref': str(member.id), 'role': role})
                            xmlobject.append(member)

                        tag = etree.Element(
                            'tag', {'k': 'type', 'v': 'multipolygon'})
                        xmlobject.append(tag)
                        if relation in featuresmap:
                            for (key, value) in sorted(featuresmap[relation].tags.items()):
                                tag = etree.Element(
                                    'tag', {'k': key, 'v': value})
                                xmlobject.append(tag)
                        xf.write(xmlobject)

        # Save last used id to file
        if saveid:
            with open(saveid, 'wb') as ff:
                ff.write(str(self.element_id_counter))
            l.info("Wrote element_id_counter '%d' to file '%s'"
                   % (self.element_id_counter, saveid))


def setup(args: list) -> dict:
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
                        help="EPSG code of source file." +
                        "If specified, overrides projection " +
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
            # Check if there was an epsg prefix and parse out the number if so
            if re.match(r'epsg', options['source_epsg'], 'i'):
                options['source_epsg'] = int(
                    re.search(r'\d{4}', options['source_epsg']).group())
            else:
                parser.error(
                    "Couldn't parse an EPSG code from the given input.")

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

    return options


def main():
    # Parse args and return as dict, along with parser for errors
    options = setup(sys.argv[1:])

    # Create memory object for data destination
    try:
        sink = OSMSink(options["source"], **options)
    except RuntimeError:
        # TODO: Useful CLI error message here
        # parser.error("Could not parse OGR data source.")
        l.exception("Could not parse OGR data source")
        sys.exit(2)

    # Save data to file
    sink.output(options["output_file"])


if __name__ == '__main__':
    main()

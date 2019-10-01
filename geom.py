# -*- coding: utf-8 -*-
#
# Copyright (c) 2012-2013 Paul Norman
# <penorman@mac.com>
# Released under the MIT license: http://opensource.org/licenses/mit-license.php

# Classes


class Geometry(object):

    def __init__(self, osm_sink):
        self.id = getNewID(osm_sink)
        self.parents = set()
        osm_sink.geometries.append(self)

    def replacejwithi(self, i, j):
        pass

    def addparent(self, parent):
        self.parents.add(parent)

    def removeparent(self, osm_sink, parent, shoulddestroy=True):
        self.parents.discard(parent)
        if shoulddestroy and len(self.parents) == 0:
            osm_sink.geometries.remove(self)

# Helper function to get a new ID


def getNewID(parent):
    parent.element_id_counter += parent.element_id_counter_incr
    return parent.element_id_counter


class Point(Geometry):
    def __init__(self, parent, x, y):
        Geometry.__init__(self, parent)
        self.x = x
        self.y = y

    def replacejwithi(self, i, j):
        pass


class Way(Geometry):
    def __init__(self, parent):
        Geometry.__init__(self, parent)
        self.points = []

    def replacejwithi(self, i, j):
        self.points = [i if x == j else x for x in self.points]
        j.removeparent(self)
        i.addparent(self)


class Relation(Geometry):
    def __init__(self, parent):
        Geometry.__init__(self, parent)
        self.members = []

    def replacejwithi(self, i, j):
        self.members = [(i, x[1]) if x[0] == j else x for x in self.members]
        j.removeparent(self)
        i.addparent(self)


class Feature(object):
    # features = []

    def __init__(self, parent):
        self.geometry = None
        self.tags = {}
        parent.features.append(self)

    def replacejwithi(self, i, j):
        if self.geometry == j:
            self.geometry = i
        j.removeparent(self)
        i.addparent(self)

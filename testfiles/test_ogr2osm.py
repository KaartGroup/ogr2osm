import unittest
from pathlib import Path

import ogr2osm
from xmldiff import main

TESTDIR = Path('testfiles').resolve()
DIFF_OPTIONS = {
    'uniqueattrs': [('node', 'lat'), ('tag', 'k')]
}
# TESTDIR = Path.cwd() / 'testfiles'
# TESTDIR = Path.cwd()

# class TestSetup(unittest.TestCase):
#     def setUp(self):
#         self.maxDiff = None


class TestParseArgs(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_parse_args(self):
        goldoptions = {
            "add_timestamp": False,
            "add_version": False,
            "debug_tags": False,
            "encoding": "utf-8",
            "force_overwrite": True,
            "id": 0,
            "idfile": None,
            "locked": False,
            "max_nodes_per_way": 1800,
            "never_download": True,
            "never_upload": False,
            "no_memory_copy": False,
            "no_upload_false": False,
            "output_file": Path('test1.osm').resolve(),
            "positive_id": False,
            "rounding_digits": 7,
            "saveid": None,
            "significant_digits": 9,
            "source": "shapefiles/test1.shp",
            "source_epsg": 3857,
            "source_proj4": None,
            "sql_query": None,
            "translation_method": None,
            "verbose": False
        }
        testargs = ['shapefiles/test1.shp',
                    '--never-download', '-f', '-e', '3857']
        testoptions = ogr2osm.setup(testargs)
        self.assertEqual(testoptions, goldoptions)

# usage

# test1


class TestOutput(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test1(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'
        goldfile = TESTDIR / 'test1.xml'
        outputfile = TESTDIR / 'test1.osm'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True)

        sink.output(outputfile)

        diff = main.diff_files(str(goldfile), str(outputfile), DIFF_OPTIONS)

        # diff should be empty if there are not substantial differences between gold and output
        self.assertFalse(diff)


# duplicatefile


class TestDuplicateFile(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_duplicate_file(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'

        sink = ogr2osm.OSMSink(testfile)

        with self.assertRaises(FileExistsError):
            sink.output(TESTDIR / 'test1.xml')

# force


class TestForce(unittest.TestCase):
    # This test should get changed because it affects behavior of other tests
    def setUp(self):
        return super().setUp()

    def test_force(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'

        sink = ogr2osm.OSMSink(testfile, force_overwrite=True)

        # This should work without raising an exception
        sink.output(TESTDIR / 'test1.xml')

# nomemorycopy


class TestNoMemoryCopy(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_no_memory_copy(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True, no_memory_copy=True)

        # This should work without raising an exception
        sink.output(TESTDIR / 'test1.xml')

# positiveid


class TestPositiveID(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_positive_id(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'
        goldfile = TESTDIR / 'positiveid.xml'
        outputfile = TESTDIR / 'positiveid.osm'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True, positive_id=True)

        sink.output(outputfile)

        diff = main.diff_files(str(goldfile), str(outputfile), DIFF_OPTIONS)

        # diff should be empty if there are not substantial differences between gold and output
        self.assertFalse(diff)

# version


class TestVersion(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_version(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'
        goldfile = TESTDIR / 'version.xml'
        outputfile = TESTDIR / 'version.osm'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True, add_version=True)

        sink.output(outputfile)

        diff = main.diff_files(str(goldfile), str(outputfile), DIFF_OPTIONS)

        # diff should be empty if there are not substantial differences between gold and output
        self.assertFalse(diff)

# timestamp


class TestTimestamp(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_timestamp(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True, add_timestamp=True)

        # This should work without raising an exception
        sink.output(TESTDIR / 'test1.xml')

# utf8


class TestUTF8(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_utf8(self):
        testfile = TESTDIR / 'shapefiles/sp_usinas.shp'
        goldfile = TESTDIR / 'utf8.xml'
        outputfile = TESTDIR / 'utf8.osm'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True)

        sink.output(outputfile)

        diff = main.diff_files(str(goldfile), str(outputfile), DIFF_OPTIONS)

        # diff should be empty if there are not substantial differences between gold and output
        self.assertFalse(diff)

# japanese


class TestJapanese(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_version(self):
        testfile = TESTDIR / 'shapefiles/japanese.shp'
        goldfile = TESTDIR / 'japanese.xml'
        outputfile = TESTDIR / 'japanese.osm'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True, encoding='shift_jis')

        sink.output(outputfile)

        diff = main.diff_files(str(goldfile), str(outputfile), DIFF_OPTIONS)

        # diff should be empty if there are not substantial differences between gold and output
        self.assertFalse(diff)

# duplicatewaynodes


class TestDuplicateWayNodes(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_duplicate_way_nodes(self):
        testfile = TESTDIR / 'duplicate-way-nodes.gml'
        goldfile = TESTDIR / 'duplicate-way-nodes.xml'
        outputfile = TESTDIR / 'duplicate-way-nodes.osm'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True)

        sink.output(outputfile)

        diff = main.diff_files(str(goldfile), str(outputfile), DIFF_OPTIONS)

        # diff should be empty if there are not substantial differences between gold and output
        self.assertFalse(diff)

# require_output_file_when_using_db_source


# class TestRequireOutputFileWhenUsingDBSource(unittest.TestCase):
#     def setUp(self):
#         return super().setUp()

#     def test_require_output_file_when_using_db_source(self):
#         teststring = "PG:dbname=test"


# require_db_source_for_sql_query


# class TestRequireDbSourceForSqlQuery(unittest.TestCase):
#     def setUp(self):
#         return super().setUp()

#     def test_require_db_source_for_sql_query(self):
#         testfile = TESTDIR / 'shapefiles/test1.shp'

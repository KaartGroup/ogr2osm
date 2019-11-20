# from tempfile import NamedTemporaryFile
import tempfile
import unittest
from pathlib import Path

import ogr2osm
from lxml import etree
from xmldiff import main, formatting


TESTDIR = Path('testfiles').resolve()
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

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True)

        _, outputfile = tempfile.mkstemp()
        outputfile = Path(outputfile).with_suffix('.osm')
        sink.output(outputfile)

        # gold_tree = etree.parse(str(goldfile))
        # output_tree = etree.parse(str(outputfile))

        # diff = main.diff_trees(gold_tree, output_tree)

        # print(diff)
        # self.assertEqual(gold_tree, output_tree)
        # diff = main.diff_files(str(outputfile), str(outputfile),
        #                        formatter=formatting.XMLFormatter())
        with outputfile.open('rb') as outputfile_stream, goldfile.open('rb') as goldfile_stream:
            self.assertEqual(outputfile_stream, goldfile_stream)


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

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True, positive_id=True)
        outputfile = tempfile.mkstemp()
        sink.output(tempfile[1])
        with outputfile.open('rb') as outputfile_stream, goldfile.open('rb') as goldfile_stream:
            self.assertEqual(outputfile_stream, goldfile_stream)

# version


class TestVersion(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_version(self):
        testfile = TESTDIR / 'shapefiles/test1.shp'
        goldfile = TESTDIR / 'version.xml'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True, add_version=True)

        with NamedTemporaryFile('rb') as outputfile, goldfile.open('rb') as goldfile_stream:
            sink.output(outputfile.name)
            self.assertEqual(outputfile, goldfile_stream)

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

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True)

        with NamedTemporaryFile('rb') as outputfile, goldfile.open('rb') as goldfile_stream:
            sink.output(outputfile.name)
            self.assertEqual(outputfile, goldfile_stream)

# japanese


# class TestJapanese(unittest.TestCase):
#     def setUp(self):
#         return super().setUp()

#     def test_version(self):
#         testfile = TESTDIR / 'shapefiles/japanese.shp'
#         goldfile = TESTDIR / 'japanese.xml'

#         sink = ogr2osm.OSMSink(
#             testfile, force_overwrite=True, encoding=shift_jis)

#         with NamedTemporaryFile('wb') as outputfile, goldfile.open('wb') as goldfile_stream:
#             sink.output(outputfile.name)
#             self.assertEqual(outputfile, goldfile_stream)

# duplicatewaynodes


class TestDuplicateWayNodes(unittest.TestCase):
    def setUp(self):
        return super().setUp()

    def test_duplicate_way_nodes(self):
        testfile = TESTDIR / 'duplicate-way-nodes.gml'
        goldfile = TESTDIR / 'duplicate-way-nodes.xml'

        sink = ogr2osm.OSMSink(
            testfile, force_overwrite=True)

        _, outputfile = tempfile.mkstemp()
        outputfile = Path(outputfile).with_suffix('.osm')
        sink.output(outputfile)

        with outputfile.open('rb') as outputfile_stream, goldfile.open('rb') as goldfile_stream:
            self.assertEqual(outputfile_stream, goldfile_stream)

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

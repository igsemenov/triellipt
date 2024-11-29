# -*- coding: utf-8 -*-
"""Tests MSH parser.
"""
import unittest
from triellipt.mshread import mshparser

MSH_SAMPLE = """
$Nodes
$EndNodes
$Elements
$EndElements
"""


class TestMSHSectionFactory(unittest.TestCase):

    SAMPLE = '$Nodes\n-----\n$EndNodes'

    @classmethod
    def setUpClass(cls):
        cls.FACTORY = cls.factory(cls.SAMPLE)

    @classmethod
    def factory(cls, content):
        return mshparser.MSHSectionFactory.with_content(content)

    def test_prefix(self):
        assert self.FACTORY.prefix == '$Nodes'

    def test_suffix(self):
        assert self.FACTORY.suffix == '$EndNodes'

    def test_get_name(self):
        assert self.FACTORY.get_name() == 'Nodes'

    def test_get_body(self):
        assert self.FACTORY.get_body() == '-----'


class TestMSHParser(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.PARSER = cls.parser(MSH_SAMPLE)

    @classmethod
    def parser(cls, content):
        return mshparser.MSHParser.with_content(content)

    @property
    def sections(self):
        return self.PARSER.getsections()

    def test_find_edges(self):
        for obj in self.PARSER.find_edges():
            assert hasattr(obj, 'group'), obj

    def test_sections_keys(self):
        for key in self.sections.keys():
            assert key in ['Nodes', 'Elements']

    def test_sections_values(self):
        for val in self.sections.values():
            assert val.body == ''


if __name__ == '__main__':
    unittest.main()

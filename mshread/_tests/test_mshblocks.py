# -*- coding: utf-8 -*-
"""Test MSH data blocks.
"""
import unittest
from triellipt.mshread import mshblocks

NODES_BODY_SAMPLE = """
0 0 0 0
0 0 0 2
1
2
1 1 1
2 2 2
"""

ELEMENTS_BODY_SAMPLE = """
0 0 0 2
0 0 1 1
0 3 4
0 0 2 2
1 4 5 6
2 7 8 9
"""


class TestNodesDecoder(unittest.TestCase):
    """Basic decoder test on nodes as an example.
    """

    @classmethod
    def setUpClass(cls):
        cls.DECODER = cls.decoder(NODES_BODY_SAMPLE)

    @classmethod
    def decoder(cls, content):
        return mshblocks.NodesDecoder.with_content(content)

    @property
    def iterator(self):
        return self.DECODER.get_iterator()

    def test_read_numdata(self):
        assert self.DECODER.read_numdata().numdata.size == 16

    def test_get_iterator(self):
        assert self.DECODER.get_iterator().numdata.size == 12

    def test_read_enty_spec(self):
        assert self.iterator.read_enty_spec().dtype.name == 'int32'

    def test_run_iterator(self):
        assert self.DECODER.run_iterator()[0].tolist() == [
            1, 1, 1, 2, 2, 2
        ]

    def test_getdata(self):
        assert self.DECODER.getdata().shape[1] == 2
        assert self.DECODER.getdata().dtype.name == 'float64'


class TestElementsDecoder(unittest.TestCase):
    """Specific decoder test for elements.
    """

    @classmethod
    def setUpClass(cls):
        cls.DECODER = cls.decoder(ELEMENTS_BODY_SAMPLE)

    @classmethod
    def decoder(cls, content):
        return mshblocks.ElementsDecoder.with_content(content)

    def test_run_iterator(self):
        assert self.DECODER.run_iterator()[0] is None
        assert self.DECODER.run_iterator()[1].tolist() == [
            1, 4, 5, 6, 2, 7, 8, 9
        ]

    def test_getdata(self):
        assert self.DECODER.getdata().shape[1] == 3
        assert self.DECODER.getdata().dtype.name == 'int32'


if __name__ == '__main__':
    unittest.main()

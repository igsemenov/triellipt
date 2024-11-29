# -*- coding: utf-8 -*-
"""Defines MSH data blocks.
"""
from abc import ABC, abstractmethod
import numpy as np

MSHBlocksError = type(
    'MSHBlocksError', (Exception,), {}
)


class MSHBlock(ABC):
    """Base data block.
    """

    NAME = None

    @classmethod
    def from_sections(cls, sections: dict):

        section = sections.get(cls.NAME)

        if section is not None:
            return cls.from_section_body(section.body)

        raise MSHBlocksError(
            f"cannot find '{cls.NAME}' in sections"
        )

    @classmethod
    @abstractmethod
    def from_section_body(cls, content):
        """Fetches block data from the section content.
        """


class MSHNodes(MSHBlock):
    """MSH nodes.
    """

    NAME = 'Nodes'

    @classmethod
    def from_section_body(cls, content):
        return NodesDecoder.with_content(content).getdata()


class MSHElements(MSHBlock):
    """MSH elements.
    """

    NAME = 'Elements'

    @classmethod
    def from_section_body(cls, content):
        return ElementsDecoder.with_content(content).getdata()


class MSHDecoder(ABC):
    """Base decoder of MSH data.
    """

    DATA_RANK = None
    DATA_TYPE = None

    def __init__(self, content=None):
        self.content = content
        self.numdata = None

    @classmethod
    def with_content(cls, content):
        return cls(content)

    @property
    def content_split(self):
        return str.split(self.content)

    def getdata(self):
        """Returns block data extracted from content.
        """

        out = self.run_iterator()
        dat = self.run_postproc(out)

        return dat

    @abstractmethod
    def run_postproc(self, itrout: list):
        """Runs postprocessing of iterator output.
        """

    def run_iterator(self):
        return list(
            self.get_iterator()
        )

    def read_numdata(self):

        if not self.content:
            raise MSHBlocksError(
                "block decoder is undefined, content is empty"
            )

        self.numdata = np.loadtxt(
            self.content_split, dtype=np.dtype(self.DATA_TYPE)
        )

        return self

    def get_iterator(self):

        self.read_numdata()

        if self.numdata.size < 4:
            raise MSHBlocksError(
                "corrupted block content, cannot read section spec"
            )

        _, self.numdata = np.split(self.numdata, [4])

        return iter(self)

    def read_enty_spec(self):

        if self.numdata.size < 4:
            raise MSHBlocksError(
                "corrupted block content, cannot read entity spec"
            )

        enty_spec, self.numdata = np.split(self.numdata, [4])
        return enty_spec.astype(int)

    @abstractmethod
    def make_enty_vars(self, enty_spec):
        """Computes entity parameters.
        """

    def read_enty_data(self, enty_vars):

        size_tags = enty_vars['size-tags']
        size_data = enty_vars['size-data']

        body_size = size_tags + size_data

        if self.numdata.size < body_size:
            raise MSHBlocksError(
                "corrupted block content, cannot read entity data"
            )

        enty_body, self.numdata = np.split(self.numdata, [body_size])
        _, enty_data = np.split(enty_body, [size_tags])

        return enty_data.copy('C')

    def __iter__(self):
        return self

    def __next__(self):

        if self.numdata.size == 0:
            raise StopIteration("data exhausted")

        enty_spec = self.read_enty_spec()
        enty_vars = self.make_enty_vars(enty_spec)
        enty_data = self.read_enty_data(enty_vars)

        if enty_vars['data-rank'] == self.DATA_RANK:
            return enty_data
        return None

    def filt_nones(self, data):
        return [
            val for val in data if val is not None
        ]


class NodesDecoder(MSHDecoder):
    """Gets MSH nodes.
    """

    DATA_RANK = 3
    DATA_TYPE = 'float'

    def make_enty_vars(self, enty_spec):

        *_, num_nodes = enty_spec

        data_rank = 3
        size_tags = num_nodes
        size_data = num_nodes * data_rank

        return {
            'size-tags': size_tags,
            'size-data': size_data,
            'data-rank': data_rank
        }

    def run_postproc(self, itrout: list):

        data = self.filt_nones(itrout)
        data = self.make_table(data)

        return data

    def make_table(self, itrout_no_nones: list):

        rank = self.DATA_RANK

        data = [
            v.reshape(v.size // rank, rank) for v in itrout_no_nones
        ]

        data = np.vstack(data)
        data = data[:, [0, 1]]  # index 2 omitted (z coordinate)

        return data.copy('C')


class ElementsDecoder(MSHDecoder):
    """Gets MSH elements.
    """

    DATA_RANK = 4
    DATA_TYPE = 'int'

    NODES_IN_ELM_TYPE = {
        1: 2,
        2: 3,
        3: 4,
        15: 1
    }

    @property
    def elm_types(self):
        return self.NODES_IN_ELM_TYPE.keys()

    def make_enty_vars(self, enty_spec):

        *_, elm_type, num_elms = enty_spec

        if elm_type not in self.elm_types:
            raise MSHBlocksError(
                f"got unknown element type {elm_type}"
            )

        elm_tag_size = 1
        nodes_in_elm = self.NODES_IN_ELM_TYPE[elm_type]

        data_rank = elm_tag_size + nodes_in_elm

        size_tags = 0
        size_data = num_elms * data_rank

        return {
            'size-tags': size_tags,
            'size-data': size_data,
            'data-rank': data_rank
        }

    def run_postproc(self, itrout: list):

        data = self.filt_nones(itrout)
        data = self.make_table(data)

        return data

    def make_table(self, itrout_no_nones: list):

        rank = self.DATA_RANK

        data = [
            v.reshape(v.size // rank, rank) for v in itrout_no_nones
        ]

        data = np.vstack(data)

        data = data[:, [1, 2, 3]]  # index 0 omitted (MSH element tag)
        data = data - 1  # because MSH numbering starts from 1

        return data.astype(int)

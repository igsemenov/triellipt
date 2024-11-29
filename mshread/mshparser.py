# -*- coding: utf-8 -*-
"""Fetches sections of MSH files.
"""
import re
import itertools as itr


def parsemsh(content):
    """Parses an MSH file.

    Parameters
    ----------
    content : str
        Content of MSH file.

    Returns
    -------
    dict
        Name-to-content map of MSH sections.

    """
    parser = MSHParser.with_content(content)
    sections = parser.getsections()
    return sections


MSHSectionError = type(
    'MSHSectionError', (Exception,), {}
)

MSHParserError = type(
    'MSHParserError', (Exception,), {}
)


class MSHSection:
    """Section of an MSH file. 
    """

    def __init__(self, name=None, body=None):
        self.name = name
        self.body = body

    @classmethod
    def from_data(cls, name, body):
        return cls(name, body)

    @classmethod
    def from_content(cls, content):
        return MSHSectionFactory().with_content(content).get_section()


class TextProcessor:

    def __init__(self, content=None):
        self.content = content

    @classmethod
    def with_content(cls, content):
        return cls(content)

    @classmethod
    def get_text_body(cls, text, pref, suff):

        rest = text.removeprefix(pref)
        body = rest.removesuffix(suff)

        return body.strip()


class MSHSectionFactory(TextProcessor):

    RE_AFFIXES = {
        'prefix': r'^\$[a-zA-Z]{1,}',
        'suffix': r'\$End[a-zA-z]{1,}$'
    }

    def __init__(self, content=None):
        self.content = content

    @classmethod
    def with_content(cls, content):
        return cls(content)

    def get_section(self):

        name, body = self.find_attrs()
        sectionobj = self.from_attrs(name, body)

        return sectionobj

    def from_attrs(self, name, body):
        return MSHSection.from_data(name, body)

    def find_attrs(self):
        yield self.get_name()
        yield self.get_body()

    def get_name(self):
        return str.lstrip(self.prefix, '$')

    def get_body(self):
        return TextProcessor.get_text_body(
            self.content, self.prefix, self.suffix
        )

    @property
    def prefix(self):
        return self.get_affix('prefix')

    @property
    def suffix(self):
        return self.get_affix('suffix')

    def get_affix(self, key):

        match = re.search(
            self.RE_AFFIXES[key], self.content
        )

        if match is not None:
            return match.group()

        raise MSHSectionError(
            f"got content with an undefined {key}"
        )


class MSHParser(TextProcessor):
    """Parser of MSH files.
    """

    RE_SECTION_EDGE = r'\$[a-zA-Z]{3,}'

    def getsections(self) -> dict:

        if not self.content:
            return {}

        snippets = self.find_snippets()
        sections = self.from_snippets(snippets)

        return sections

    def from_snippets(self, snippets):

        sections = list(
            map(MSHSection.from_content, snippets)
        )

        return {
            o.name: o for o in sections
        }

    def find_snippets(self):

        edges = self.find_edges()
        snippets = self.from_edges(edges)

        return snippets

    def find_edges(self):
        return list(
            re.finditer(self.RE_SECTION_EDGE, self.content)
        )

    def from_edges(self, edges):

        if len(edges) % 2 != 0:
            raise MSHParserError(
                "corrupted MSH file, unmatched section ends found"
            )

        pairs = self.make_starts_ends(edges)
        snippets = self.from_starts_ends(pairs)

        return snippets

    def make_starts_ends(self, edges):

        starts = itr.islice(edges, 0, None, 2)
        ends = itr.islice(edges, 1, None, 2)

        return zip(starts, ends)

    def from_starts_ends(self, starts_ends):
        return list(
            itr.starmap(self.fetch_snippet, starts_ends)
        )

    def fetch_snippet(self, start, end):
        return self.content[
            start.start():end.end()
        ]

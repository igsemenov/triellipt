# -*- coding: utf-8 -*-
"""Documentation maker — modules.
"""
import docspyer
import triellipt

DOCPATH = 'docs/sources'

MODULES = [
    triellipt.geom,
    triellipt.mesher,
    triellipt.trimesh,
    triellipt.mshread
]

config = {
    'docsname': 'triellipt',
    'hostname': 'triellipt',
    'modrefs': True,
    'clsverbs': 2
}

docspyer.docmods(
    MODULES, DOCPATH, **config
)

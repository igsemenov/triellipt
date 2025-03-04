# -*- coding: utf-8 -*-
"""Documentation maker.
"""
import runpy
import docspyer

CFGPATH = 'docs/configs'
SRCPATH = 'docs/sources'
DOCPATH = 'docs/build'

with open(f'{CFGPATH}/logo.html', 'r', encoding='utf-8') as logo_file:
    logo = logo_file.read()
logo += '<p id="logo-title">triellipt</p>'

config = {
    'doclogo': logo,
    'swaplinks': True,
    'codeblocks': True,
    'extracss': f'{CFGPATH}/theme.css'
}

runpy.run_path('_docmods.py')

docspyer.builddocs(SRCPATH, DOCPATH, **config)

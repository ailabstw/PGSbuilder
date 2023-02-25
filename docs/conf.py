# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'PGSBuilder'
copyright = '2023, Taiwan AI Labs'
author = 'Ko-Han Lee'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'myst_parser',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

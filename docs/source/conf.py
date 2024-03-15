# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MuDirac documentation'
copyright = '2024, Simone Sturniolo'
author = 'Simone Sturniolo'
release = 'November 16, 2020'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    #'sphinx.ext.autodoc',
    #'breathe',
    'sphinx_rtd_theme',
    # Other extensions...
]


templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_context = {
    "display_github": True, # Integrate GitHub
    "github_user": "muon-spectroscopy-computational-project", # Username
    "github_repo": "mudirac", # Repo name
    "github_version": "main", # Version
    "conf_py_path": "/docs/source", # Path in the checkout to the docs root
}
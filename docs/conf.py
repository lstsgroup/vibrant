# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'vibrant'
copyright = '2024, Ekin Esme Bas'
author = 'Ekin Esme Bas'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx_wagtail_theme", 'myst_parser']

exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_wagtail_theme'
html_theme_options = dict(
    project_name = "vibrant",    
    github_url = "https://github.com/lstsgroup/vibrant/tree/main/src",
    footer_links = ",".join([
        "Golze Group|https://golzegroup.org/",
        "Contact|https://golzegroup.org/contact/",
    ]),
    logo = "all_img/logo_vibrant_dark.svg",
    logo_alt = "vibrant",
    logo_height = 100,
    logo_url = "/",
    #logo_width = 45,
)
html_static_path = ['_static']

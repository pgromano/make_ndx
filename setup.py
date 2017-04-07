from setuptools import setup, Extension, find_packages
import numpy

setup(
    name = 'make_ndx',
    author = 'Pablo Romano',
    author_email = 'promano@uoregon.edu',
    description = 'Python API indexing topology files',
    version = '0.1',
    url = 'https://github.com/pgromano/make_ndx',

    packages = ['make_ndx'],
    install_requires=[
        'numpy',
        'parmed'
    ],
    ext_modules=extensions,
    zip_safe = False,
)

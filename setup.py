"""Setuptools for building the extension package"""
from setuptools import setup, Extension
import numpy.distutils.misc_util

PACKAGE_NAME = 'ar'

setup(
    name=PACKAGE_NAME,
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    version='0.0.1',
    author='Rhys Ulerich',
    author_email='rhys.ulerich@gmail.com',
    license='MPL2',
    classifiers=['Development Status :: 3 - Alpha'],
    ext_modules=[
        Extension(
            PACKAGE_NAME,
            depends=['ar.hpp'],
            sources=['ar-python.cpp'],
            include_dirs=[
                numpy.distutils.misc_util.get_numpy_include_dirs()
            ],
        )
    ],
)

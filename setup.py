"""Setuptools for building the extension package"""
from setuptools import setup, Extension
import numpy

PACKAGE_NAME = 'ar'

setup(
    ext_modules=[
        Extension(
            PACKAGE_NAME,
            depends=['ar.hpp'],
            sources=['ar-python.cpp'],
            include_dirs=[numpy.get_include()],
        )
    ],
)

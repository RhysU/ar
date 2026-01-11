"""Setuptools for building the extension package"""
from setuptools import setup, Extension
import numpy
import os

PACKAGE_NAME = 'ar'

# Get version from environment variable set by Makefile
version = os.environ.get('ARSEL_VERSION', 'unknown')

setup(
    ext_modules=[
        Extension(
            PACKAGE_NAME,
            depends=['ar.hpp'],
            sources=['ar-python.cpp'],
            include_dirs=[numpy.get_include()],
            extra_compile_args=[f'-DARSEL_VERSION="{version}"'],
        )
    ],
)

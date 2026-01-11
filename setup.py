"""Setuptools for building the extension package"""
from setuptools import setup, Extension
import numpy
import subprocess

PACKAGE_NAME = 'ar'

# Get version from git describe, similar to Makefile
try:
    version = subprocess.check_output(
        ['git', 'describe', '--always', '--dirty'],
        stderr=subprocess.DEVNULL,
        text=True
    ).strip()
except (subprocess.CalledProcessError, FileNotFoundError):
    version = 'unknown'

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

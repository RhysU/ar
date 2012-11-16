from distutils.core import setup, Extension

import numpy.distutils.misc_util

ar = Extension('ar', sources = ['ar-python.cpp'])

setup(name = 'ar',
      url = 'http://rhysu.github.com/ar/',
      ext_modules = [ar],
      include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs())

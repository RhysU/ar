import distutils.core
import numpy.distutils.misc_util

distutils.core.setup(name = 'ar',
    url = 'http://rhysu.github.com/ar/',
    ext_modules = [
        distutils.core.Extension(
            'ar',
            depends = ['ar.hpp'],
            sources = ['ar-python.cpp']
        )
    ],
    include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
)

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [Extension("cstvenant", ["cstvenant.pyx"],
                        extra_compile_args=["-std=c99"])]

setup(ext_modules=cythonize(extensions, include_path=[numpy.get_include()]))

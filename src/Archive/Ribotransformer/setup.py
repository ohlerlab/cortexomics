from setuptools import setup
from Cython.Build import cythonize
import numpy


setup(
    name='arraychoose',
    ext_modules=cythonize("arraychoose.pyx"),
    include_dirs=[numpy.get_include()],
    zip_safe=False
)
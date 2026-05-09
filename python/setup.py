from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension ("pyatsat",
              ['pyatsat.pyx', 'PythonStream.cpp'],
              language = 'c++',
              extra_compile_args = ["-std=c++11", "-fopenmp"],
              libraries = ['atsat'],
              include_dirs = ['../planet', '../lib', '../tinyxml2'],
              library_dirs = ['..'],
              extra_link_args = ['-fopenmp'],
              )]

setup(
    name = 'pyatsat',
    ext_modules = cythonize(extensions, language_level = "2"),
    )

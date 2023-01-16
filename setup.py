from setuptools import setup, Extension

setup(name='kmeansspk',
      version='1.0',
      description='kmeans final',
      ext_modules=[Extension('kmeansspk', sources=['spkmeans.c', 'spkmeansmodule.c'])])
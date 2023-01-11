from setuptools import setup

setup(name='bacmapping',
      version='0.1',
      description='A set of tools for mapping and understanding bacterial artificial chromosomes',
      url='http://github.com/Ewinden/BACMapping',
      author='Eamon Winden',
      author_email='ewinden@wisc.edu',
      license='',
      packages=['bacmapping'],
      install_requires=[
     	'pandas',
	'biopython',
	'multiprocess',
	'matplotlib'
	],
      zip_safe=False)

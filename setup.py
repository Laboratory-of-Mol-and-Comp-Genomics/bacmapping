from setuptools import setup

setup(name='bacmapping',
      version='1.0',
      description='A set of tools for mapping and understanding bacterial artificial chromosomes',
      url='http://github.com/Ewinden/bacmapping',
      author='Eamon Winden',
      author_email='ewinden@wisc.edu',
      license='MIT',
      packages=['bacmapping'],
      install_requires=[
     	'pandas==1.5.2',
	'biopython==1.80',
	'multiprocess==0.70.14',
	'matplotlib==3.6.3'
	],
      zip_safe=False)

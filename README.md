# bacmapping
This is a set of tools used in production and exploration of bacterial artificial chromosome restriction maps.

## Dataset location
This dataset is kept at gigaDB at , it can also be built locally using the scripts in installing

## Installation

We recommend creating a new enviroment

```bash
conda create -n bacmapping
conda activate bacmapping
```

This package requires:
- numpy
- pandas  1.5.2
- biopython  1.80
- matplotlib  3.6.3
- multiprocess  0.70.14

version control is probably unnecessary but included all the same

```bash
conda install -c conda-forge pandas=1.5.2 biopython=1.80 multiprocess=0.70.14 matplotlib=3.6.3
```

Then clone and install this github repository on Linux

```bash
git clone https://github.com/Ewinden/bacmapping/
cd bacmapping
pip install .
```
## Example
If you wish to go through the functions with an example in a jupyter notebook, look in the Examples folder for the files "bacmappingexample.ipynb" or "Ch19_bacmapping_example.ipynb". Open any this notebooks and use the bacmapping environment.


## Functions
### Main pipeline functions
- bacm.getNewClones(download = True, onlyType=True, vtype='BAC', chunk_size = 5000, email='example@website.com')
  - Download files from the FTP server
    - download is a boolean whether to download sequence files which are required, but may exist locally already
    - onlyType allows the user to only download details and sequence for one type of clone, ie BACs
    - vtype determines the type of clone
    - chunk_size determines the amount of lines to read into pandas at once, larger is faster but requires more memory
    - email should be set to the user email to allow the NIH servers to know who is accessing the database
    
- mapSequencedClones(cpustouse=1, maxcuts=50, chunk_size=500)
  - Maps all the clones which are insert sequenced
    - cpustouse determines the number of cores to use when running multiprocessing
    - maxcuts determines the maximum number of cuts with a specific enzyme before the map is truncated, to save space
    - chunk_size determines the amount of lines to read into pandas at once, larger is faster but requires more memory

- mapPlacedClones(include_libraries=True, cpustouse=1, maxcuts=50, chunk_size=500)
  - Maps all the clones which are end-sequenced and mapped to the reference genome
    - include_libraries determines whether or not to use only libraries referenced in vtype in getNewClones
    - cpustouse determines the number of cores to use when running multiprocessing
    - maxcuts determines the maximum number of cuts with a specific enzyme before the map is truncated, to save space
    - chunk_size determines the amount of lines to read into pandas at once, larger is faster but requires more memory
        
### Functions for statistics regarding the database
- countPlacedBACs()
  - returns the number of placed BACs and library identifiers
    
- getCoverage(chunk_size = 5000)
  - determines coverage of the primary assembly for the genome
    - chunk_size determines the amount of lines to read into pandas at once, larger is faster but requires more memory
    
- getAverageLength()
  - returns average length of mapped clones, removing extreme size cases
    
- getSequencedClonesStats(include_libraries=True)
  - returns statistics for the sequenced clone database
    - include_libraries determines whether or not to use the list generated from narrowDownLibraries

### Exploring the library
- getRightIsoschizomer(enzyme)
  - Given an enzyme name, returns the enzyme name and Bio.restriction class which corresponds to the isoschizomer which is in the database
    - enzyme is a string of the common name of an enzyme included in the list, such as "HindIII"
    
- drawMap(name, enzyme, circular=False)
  - Given the name of a BAC and an enzyme, draws a map
    - name is a string containing the common name of a BAC, including library, such as "RP11-168H2"
    - enzyme is a string of the common name of an enzyme included in the list, such as "HindIII"
    - circular is whether to return a circular or linear plot

- getSequenceFromName(name)
  - Given the name of a BAC, tries to return the sequence of that insert
    - name is a string containing the common name of a BAC, including library, such as "RP11-168H2"

- getSequenceFromLoc(chrom,start,end)
  - Given a chromosome, start and end location, returns sequence of that location
    - chrom is the chromosome to pull from
    - start is the location of the first base in the chromosome sequence
    - end is the location of the base after the last base in the chromosome sequence
    
- getMapFromName(name)
  - Given the name of a BAC, tries to return the all the restriction maps for that name
    - name is a string containing the common name of a BAC, including library, such as "RP11-168H2"

- getMapsFromLoc(chrom,start,end, inclusive=True)
  - Given a chromosome, start and end location, returns all the maps in that region
    - chrom is the chromosome to pull from
    - start is the location of the first base in the chromosome sequence
    - end is the location of the base after the last base in the chromosome sequence
    - inclusive determines whether BACs that overlap the start and end locations but are not wholely inside the range should be included (True means yes)
        
- getMaps(name)
  - Given the name of a BAC, returns all restriction maps
    - name is a string containing the common name of a BAC, including library, such as "RP11-168H2"
   
- getRestrictionMap(name,enzyme)
  - Given the name of a BAC and an enzyme, returns the cut locations
    - name is a string containing the common name of a BAC, including library, such as "RP11-168H2"
    - enzyme is a string of the common name of an enzyme included in the list, such as "HindIII"

- findPairsFromName(name, longestoverlap, shortestoverlap,)
    - Given a name and overlap parameters, finds all pairs of BACs where enzymes which linearize both BACs and results in ends with overlap parameters
    - name is a string containing the common name of a BAC, including library, such as "RP11-168H2"
    - longestoverlap determines the longest distance between the two enzymes, or the overlap if doing gibson synthesis
    - shortestoverlap determines the shortest distance between the two enzymes, or the overlap if doing gibson synthesis; setting this below zero allows for the same enzyme to cut both
   
- makePairs(cpustouse=1,longestoverlap=200,shortestoverlap=20)
  - Generates a database of pairs of BACs which have overlaps generated by restriction enzymes that linearize the BACs
    - cpustouse determines the number of cores to use when running multiprocessing
    - longestoverlap determines the longest distance between the two enzymes, or the overlap if doing gibson synthesis
    - shortestoverlap determines the shortest distance between the two enzymes, or the overlap if doing gibson synthesis; setting this below zero allows for the same enzyme to cut both

### Internal functions
- getRow(name)
  - Given the name of a BAC, tries to return the set of restriction maps for that BAC
    - name is a string containing the common name of a BAC, including library, such as "RP11-168H2"

- splitAttributesWithMids(Ser,middles)
  - Required to split the weird way details are given
    
- getCuts(row)
  - Given a row from the mapping functions, returns all the restriction maps

- onlySingleCutters(row)
  - Given a row of restriction maps, returns the enzymes that only cut once

- findPairs(row)
  - Given a row for a specific BAC as well as overlap and other details, finds possible BACs with acceptable overlap and restriction sites

- splitAttributes(ser)
  - Degenerated, but perhaps useful, function to split features for a BAC
    
- openSeqgetCuts(row)    
  - Given a row, opens the sequence and returns all the restriction maps. This has been deprecated in the primary path, but could be useful

- narrowDownLibraries()
  - Incorporated into getNewClones, generates a file with only BAC libraries. Can be modified to include other libraries

- makeIndexFiles(loc)
  - makes index files for all BACs in a given library for easy finding
    - loc is the location of the library file used

## Use
Ensure you are on a computer/ server that can handle a large throughput and can be left for some time to download/ process everything as well as some space to save the database.
In python, the main pipeline is run as

```python
import bacmapping as bmap
from time import perf_counter

times = []
times.append(perf_counter())
print('starting downloads')
bmap.getNewClones(download = True, onlyType=True, vtype='BAC', chunk_size = 5000, email='user@github.com/ewinden/bmapapping')
print('This took ' + str(perf_counter()-times[-1]) + ' seconds')
times.append(perf_counter())
print('starting sequenced maps')
bmap.mapSequencedClones(cpustouse=4,chunk_size=5000)
print('This took ' + str(perf_counter()-times[-1]) + ' seconds')
times.append(perf_counter())
print('starting placed maps')
bmap.mapPlacedClones(cpustouse=4,chunk_size=5000)
print('This took ' + str(perf_counter()-times[-1]) + ' seconds')

#the fundamental database has been established, the next lines describe the database

times.append(perf_counter())
print('starting stats')
bmap.countPlacedBACs()
bmap.getCoverage()
bmap.getAverageLength()
bmap.getSequencedClonesStats()
print('This took ' + str(perf_counter()-times[-1]) + ' seconds')
times.append(perf_counter())
bmap.makePairs(cpustouse=8)
print('This took ' + str(perf_counter()-times[-1]) + ' seconds')
times.append(perf_counter())
print(times)
```


Further examples on the use of Bacmapping can be found in the Jupyter notebook located in the examples folder.

## Anatomy of the database
Once the database is generated, or on the gigaDB database, restriction maps are saved first by type (insert-sequenced in sequenced or end-sequenced in placed) then library, then chromosome. Within the database are an initial set of folders and csv files, the csv files include statistics regarding the database. The details folder contains data directly downloaded from CloneDB and modified to move headers out of the files. The sequences folder contains all of the sequences included in the downloaded CloneDB files. Each chromosome csv is actually a tsv where each line is an individual BAC along with details for the BAC and then the set of cut sites for each enzyme included. Regarding cut sites, entries are either empty ([]) meaning there are no cutsites, overflow meaning the enzyme cuts more than 50 times in the BAC, or a python formatted list ([145,2352,6546]) where each cut site is included in the list. 


## Authors
This project was produced in the Laboratory of Molecular and Computational Genomics at University of Wisconsin- Madison
    
Eamon Winden - ewinden@wisc.edu         
Alejandro Vasquez-Echeverri provided guidance in git, multiprocessing, python, ssh...

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.

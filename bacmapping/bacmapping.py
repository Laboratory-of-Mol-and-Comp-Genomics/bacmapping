import os
import pandas as pd
import csv
import numpy as np
from ftplib import FTP
from Bio import SeqIO
from Bio import Entrez
from Bio import Restriction as rst
import multiprocessing
import ast
import matplotlib.pyplot as plt
import math
from itertools import repeat

#Download files from the FTP server
def getNewClones():
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Setup names
    cloneacstate = 'clone_acstate_'+taxid+'.out'
    librarys = 'library_'+taxid+'.out'
    ucname = version + '.unique_concordant.gff'    

    #Setup folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    if os.path.isdir(clonesDetails) == False:
        os.mkdir(clonesDetails)
    if os.path.isdir(clonesSequences) == False:
        os.mkdir(clonesSequences)

    #Login to NCBI FTP to download details files
    with FTP("ftp.ncbi.nih.gov") as ftp:
        ftp.login()
        ftp.cwd('repository/clone/reports/Homo_sapiens/')
        filenames = ftp.nlst()
        ucst = [x for x in filenames if ucname in x]
        downloads = ucst + [librarys,cloneacstate]
        for f in downloads:
            p = os.path.join(clonesDetails,f)
            if os.path.exists(p) == False:
                ftp.retrbinary("RETR " + f ,open(clonesDetails + f, 'wb').write)
                
    #Get accession list
    cloneacstatepath = os.path.join(clonesDetails,cloneacstate)
    libraryspath = os.path.join(clonesDetails,librarys)
    ucstpaths= [os.path.join(clonesDetails,x) for x in ucst]
    
    #Get accessions of sequenced clones from clone_acstate
    seqdclones = pd.read_csv(cloneacstatepath, sep='\t')
    finseqclones = seqdclones[(seqdclones['Stdn']=='Y') & (seqdclones['CloneState']=='fin')]
    finseqnames = finseqclones['CloneName'].unique()
    finseqaccs = finseqclones['Accession'].unique()
    finseqclones.to_csv(os.path.join(clonesDetails,'clone_acstate_'+taxid+'_onlyfinished.out'),sep='\t',index=False)
    
    #Split unique_concordant files into a header and details
    for uc in ucstpaths:
        ucinfo = uc + '.info.txt'
        if os.path.exists(ucinfo) == True:
            continue
        lines = open(uc).readlines()
        with open(ucinfo, 'w') as ui:
            ui.writelines(lines[:7])
        with open(uc, 'w') as ci:
            ci.writelines(['seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes']+lines[7:])
            
    #Get accessions for placed clones
    allplacedaccs = []
    for uc in ucstpaths:
        uccur = pd.read_csv(uc, sep='\t')
        plasecaccs = uccur['seqid'].unique()
        [allplacedaccs.append(x) for x in list(plasecaccs) if x not in allplacedaccs]
        
    #Make superfile and save accessions lists
    allaccs = allplacedaccs + list(finseqaccs)
    with open(os.path.join(clonesSequences,'Accessions.csv'), 'w') as accessions:
        accessions.writelines('\n'.join(allaccs) + '\n')
    with open(os.path.join(clonesSequences,'PlacedAccessions.csv'), 'w') as accessions:
        accessions.writelines('\n'.join(allplacedaccs) + '\n')
    with open(os.path.join(clonesSequences,'SequencedAccessions.csv'), 'w') as accessions:
        accessions.writelines('\n'.join(finseqaccs) + '\n')
        
    #Download all sequences
    Entrez.email = "student@wisc.edu"  # Always tell NCBI who you are
    for accession in allaccs:
        save = os.path.join(clonesSequences,accession+'.fasta')
        net_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        out_handle = open(save, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()

#Get only BAC libraries        
def narrowDownLibraries():
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Files
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    
    #set up important files to map and get details
    libraryDetails = os.path.join(clonesDetails,'library_'+taxid+'.out')
    librariesToInclude = os.path.join(clonesDetails, 'includedLibraries.csv')
    
    #Read, shorten, write
    fullp = pd.read_csv(libraryDetails, sep='\t')
    cut = fullp[fullp['vtype'] == 'BAC']
    uselibs = cut['libabbr']
    uselibs.to_csv(librariesToInclude, index = False, header = False)

#using the list of spots to split, cut up attributes row        
def splitAttributesWithMids(Ser,middles):
    new = str(Ser).split(';')
    clipped = [x[middles[i]+1:] for i,x in enumerate(new)]
    df = pd.Series(clipped)
    return(df)    

#Given a row, gets the sequence and cuts it
def getCuts(row):
    clonesSequences = row[7]
    maxcuts = row[6]
    accPath = os.path.join(clonesSequences,row[4]+'.fasta')
    if os.path.isfile(accPath) == False:
        return('NoA')
    record_iter = SeqIO.parse(open(accPath), "fasta")
    record = list(record_iter)[0]
    seq = record.seq[int(row[2]):int(row[3])]
    cuts = shortComm.search(seq)
    for key in cuts:
        val = cuts[key]
        if type(val) == int:
            lengthv = 1
        else:
            lengthv = len(val)
        if lengthv > maxcuts:
            cuts[key] = 'NA'
    chromosome = row[1]
    if row[1] == '0':
        de = str(record.description)
        chromosome = de[de.find('chromosome ')+11:de.find(',')]
    cuts['Name'] = row[0]
    cuts['Chrom'] = chromosome
    cuts['Start'] = row[2]
    cuts['End'] = row[3]
    cuts['Accession'] = row[4]
    cuts['Library'] = row[5]
    return(cuts)

#Map the clones detailed in clone_acstate_taxid 
def mapSequencedClones(include_libraries=True, cpustouse=1, maxcuts=50, chunk_size=500):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
     #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'sequenced')
    enzfile = os.path.join(cwd,'shortComm.csv')
    if os.path.isdir(clonesMapsBase) == False:
        os.mkdir(clonesMapsBase)
    if os.path.isdir(clonesMaps) == False:
        os.mkdir(clonesMaps)
    
    #set up important files to map and get details
    sequencedClonesDetails = os.path.join(clonesDetails,'clone_acstate_'+taxid+'_onlyfinished.out')
    included_libraries = os.path.join(clonesDetails,'includedLibraries.csv')
    if include_libraries == True:
        with open(included_libraries) as incl:
            reader = csv.reader(incl)
            included_list = [row[0] for row in reader]
            
    #Open a shortened enzyme file, without isoschizomers
    with open(enzfile,"r") as rit:
        rite = csv.reader(rit, lineterminator = '\n', delimiter=",")
        for row in rite:
            shortC = row
    global shortComm
    shortComm = rst.RestrictionBatch(shortC)
    knames = ['Name','Library','Chrom','Start','End','Accession'] + list(shortComm)
    
    sequencedClones = pd.read_csv(sequencedClonesDetails, sep='\t')
    sequencedClones['start'] = '0'
    if include_libraries == True:
        sequencedClones = sequencedClones[sequencedClones['LibAbbr'].isin(included_list)]
    row = zip(sequencedClones['CloneName'], sequencedClones['Chrom'], sequencedClones['start'],
              sequencedClones['SeqLen'], sequencedClones['Accession'], sequencedClones['LibAbbr'],
             repeat(maxcuts), repeat(clonesSequences))
    with multiprocessing.Pool(cpustouse) as p:
        for result in p.imap_unordered(getCuts, row):
            if result == 'NoA':
                continue
            folder = os.path.join(clonesMaps, result['Library'])
            if os.path.isdir(folder) == False:
                os.mkdir(folder)
            file = os.path.join(folder, str(result['Chrom'])+'.csv')
            if os.path.isfile(file) == False:
                with open(file, "a") as fwri:
                    writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter=",", fieldnames = knames)
                    writer.writeheader()
            with open(file, "a") as fwri:
                writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter=",", fieldnames = knames)
                writer.writerow(result)

#Map the end-sequenced clones based on placement details
def mapPlacedClones(include_libraries=True, cpustouse=1, maxcuts=50, chunk_size=500):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    enzfile = os.path.join(cwd,'shortComm.csv')
    if os.path.isdir(clonesMapsBase) == False:
        os.mkdir(clonesMapsBase)
    if os.path.isdir(clonesMaps) == False:
        os.mkdir(clonesMaps)
    
    #set up important files to map and get details
    placedClonesDetailFs = [x for x in os.listdir(clonesDetails) if x[-3:] == 'gff']
    included_libraries = os.path.join(clonesDetails,'includedLibraries.csv')
    if include_libraries == True:
        with open(included_libraries) as incl:
            reader = csv.reader(incl)
            included_list = [row[0] for row in reader]
            
    #Open a shortened enzyme file, without isoschizomers
    with open(enzfile,"r") as rit:
        rite = csv.reader(rit, lineterminator = '\n', delimiter=",")
        for row in rite:
            shortC = row
    global shortComm
    shortComm = rst.RestrictionBatch(shortC)
    knames = ['Name','Library','Chrom','Start','End','Accession'] + list(shortComm)
    
    #get the location to return to if necesary
    startline = 0
    startlibf = os.path.join((clonesMaps),'current_library.txt')
    startchunkf = os.path.join(clonesMaps,'current_chunk.txt')
    if os.path.isfile(startlibf) == True and os.path.isfile(startchunkf) == True:
        with open(startlibf, 'r') as sli:
            startlib = sli.read()
        startlibn = placedClonesDetailFs.index(startlib)
        with open(startchunkf, 'r') as chi:
            startchunk = chi.read()
            startchunk = int(startchunk)

        first = True
    else:
        startlibn = 0
        first = False

    #make a list based on the clone libraries to include
    cPlacedClonesDetailFs= placedClonesDetailFs.copy()
    if include_libraries == True:
        cPlacedClonesDetailFs=[]
        for pCD in placedClonesDetailFs[startlibn:]:
            library = pCD[:pCD.find('.')]
            if library in included_list:
                cPlacedClonesDetailFs.append(pCD)

    #loop through libaries, loop through chunks clones, process them
    print(str(len(cPlacedClonesDetailFs)) + ' libraries')
    for pCD in cPlacedClonesDetailFs[startlibn:]:
        library = pCD[:pCD.find('.')]
        folder = os.path.join(clonesMaps, library)
        with open(os.path.join(startlibf), 'w') as libri:
            libri.write(pCD)
        if os.path.isdir(folder) == False:
            os.mkdir(folder)
        fpCD = os.path.join(clonesDetails,pCD)
        stahp = False
        for placedClones in pd.read_csv(fpCD, sep='\t', chunksize=1): #load first line to get the weird attributes list sorted
            if placedClones.empty:
                stahp = True
                break
            firstline = str(placedClones.iloc[0]['attributes']).split(';')
            cnames = [x[:x.find('=')] for x in firstline]
            middles = [x.find('=') for x in firstline]
            break
        if stahp == True:
            continue
        chunknum = 0
        if first == True:
            first = False
            startline = chunk_size * startchunk
        for placedClones in pd.read_csv(fpCD, sep='\t', chunksize=chunk_size, skiprows = startline): #load chunk_size many lines
            startline = 0
            print(str(library)+', '+str(chunknum) + ', ' + str(len(placedClones)))
            with open(os.path.join(startchunkf), 'w') as chri:
                chri.write(str(chunknum))
            chunknum+=1
            newcols = placedClones['attributes'].apply(splitAttributesWithMids, middles=middles)
            newcols.columns = cnames
            placedClones = pd.concat([placedClones,newcols], axis=1)
            placedClones['Library'] = library
            accs = (placedClones['seqid']).unique()
            splitPlacedClones = [placedClones[placedClones['seqid']==x] for x in accs] #split clones into separate dfs by accession
            for oneSeqPlacedClones in splitPlacedClones:
                accession = oneSeqPlacedClones.iloc[0]['seqid']
                accPath = os.path.join(clonesSequences,accession+'.fasta')
                if os.path.isfile(accPath) == False:
                    continue
                record_iter = list(SeqIO.parse(open(accPath), "fasta"))
                currentseq = record_iter[0]
                de = currentseq.description
                curChrom = de[de.find('chromosome ')+11:de.find(',')]
                row = zip(oneSeqPlacedClones['Name'],repeat(curChrom),oneSeqPlacedClones['start'],
                          oneSeqPlacedClones['end'],oneSeqPlacedClones['seqid'],repeat(library),
                         repeat(maxcuts), repeat(clonesSequences))
                with multiprocessing.Pool(cpustouse) as p:
                    for result in p.imap_unordered(getCuts, (row)):
                        if result == 'NoA':
                            continue
                        file = os.path.join(folder, curChrom+'.csv')
                        if os.path.isfile(file) == False:
                            with open(file, "w") as fwri:
                                writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter="\t", fieldnames = knames)
                                writer.writeheader()
                        with open(file, "a") as fwri:
                            writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter="\t", fieldnames = knames)
                            writer.writerow(result)
    os.remove(startlibf)
    os.remove(startchunkf)

#Draw a circular or linear map of an enzyme
def drawMap(row,enzyme,circular = True):
    nums = ast.literal_eval(list(r[enzyme])[0])
    totallength = int(row['End']) - int(row['Start'])
    partnums = np.array(nums)/totallength
    if circular == False:
        figure, axes = plt.subplots()
        axes.set_aspect( 1 )
        plt.axis('off')
        plt.plot([0, 1],[0.5,0.5],color='k',linestyle='-',linewidth=2)
        plt.text(0,0.53,'0')
        plt.text(1,.53,str(totallength))
        axes.set_xlim(-0.1,1.1)
        axes.set_ylim(0.25,.75)
        plt.axis('off')
        plt.plot([0, 0],[0.475,0.525],color='k',linestyle='-',linewidth=2)
        plt.plot([1, 1],[0.475,0.525],color='k',linestyle='-',linewidth=2)
        for v in range(0,totallength,2000):
            piv = v/totallength
            if v%10000 == 0:
                plt.plot([piv, piv], [0.48, .52], color='grey', linestyle='-', linewidth=1)
            else:
                plt.plot([piv, piv], [0.49, .51], color='grey', linestyle='-', linewidth=1)
        flip = -1
        for v in nums:
            piv = v/totallength
            flip = flip * -1
            plt.plot([piv, piv], [.475, .525], color='blue', linestyle='-', linewidth=2)
            rotate = 45 * flip
            plt.text(piv + (flip*0.01), 0.49 + (flip * 0.05),str(v), rotation = rotate, rotation_mode = 'anchor', transform_rotates_text=True)
    else:
        figure, axes = plt.subplots() 
        cc = plt.Circle(( 0 , 0 ), 1 , fill= False) 
        axes.set_aspect( 1 )
        axes.add_patch(cc)
        axes.set_xlim(-1.5,1.5)
        axes.set_ylim(-1.5,1.5)
        plt.axis('off')
        plt.plot([0, 0], [0.8, 1], color='grey', linestyle='-', linewidth=2)
        for v in range(0,totallength,2000):
            piv = (v/totallength) * 2 * math.pi
            x = math.sin(piv)
            y = math.cos(piv)
            if v%10000 == 0:
                plt.plot([0.9*x, x], [0.9* y, y], color='grey', linestyle='-', linewidth=1)
            else:
                plt.plot([0.95*x, x], [0.95* y, y], color='grey', linestyle='-', linewidth=1)
        for v in nums:
            piv = (v/totallength) * 2 * math.pi
            x = math.sin(piv)
            y = math.cos(piv)
            plt.plot([.95*x, 1.1*x], [.95*y, 1.1*y], color='blue', linestyle='-', linewidth=2)
            rotate = (math.degrees(piv)-90)*-1
            plt.text(1.15*x,1.15*y,str(v), rotation = rotate, rotation_mode = 'anchor', transform_rotates_text=True)
    return(plt)

#Count BACs
def countPlacedBACs():
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    enzfile = os.path.join(cwd,'shortComm.csv')
    
    #Count
    totalnumber = 0
    libnumbers = []
    libraries = os.listdir(clonesMaps)
    for library in libraries:
        libnumber = 0
        libpath = os.path.join(clonesMaps, library)
        if os.path.isdir(libpath) == False:
            continue
        chroms = os.listdir(libpath)
        for chrom in chroms:
            chrompath = os.path.join(libpath,chrom)
            locmaps = pd.read_csv(chrompath, sep = '\t', low_memory=False)
            chromnumber = len(locmaps)
            libnumber+=chromnumber
        libnumbers.append([library,libnumber])
        totalnumber += libnumber
    libnumbers.append(['total',totalnumber])
    outputf = os.path.join(cwd, 'counts.csv')

    with open (outputf,'w') as file:
        w = csv.writer(file)
        w.writerows(libnumbers)

def splitAttributes(ser):
        middles = [x.find('=') for x in ser]
        new = str(ser).split(';')
        clipped = [x[middles[i]+1:] for i,x in enumerate(new)]
        df = pd.Series(clipped)
        return(df)        

#Get coverage of the genome
def getCoverage(include_libraries=True):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    enzfile = os.path.join(cwd,'shortComm.csv')
    
    placedAccessions = os.path.join(clonesSequences,'PlacedAccessions.csv')
    libraryDetails = os.path.join(clonesDetails,'library_'+taxid+'.out')
    sequencedClonesDetails = os.path.join(clonesDetails,'clone_acstate_'+taxid+'_onlyfinished.out')
    placedClonesDetailFs = [x for x in os.listdir(clonesDetails) if x[-3:] == 'gff']
    included_libraries = os.path.join(clonesDetails,'includedLibraries.csv')
    if include_libraries == True:
        with open(included_libraries) as incl:
            reader = csv.reader(incl)
            included_list = [row[0] for row in reader]

    if include_libraries == True:
        cPlacedClonesDetailFs=[]
        for pCD in placedClonesDetailFs:
            library = pCD[:pCD.find('.')]
            if library in included_list:
                cPlacedClonesDetailFs.append(pCD)

    with open (placedAccessions, 'r') as fi:
        allaccs = [x.rstrip() for x in fi]
    allchrs = ['0'] * len(allaccs)
    accnames = ['0'] * len(allaccs)

    for pCD in cPlacedClonesDetailFs:
        library = pCD[:pCD.find('.')]
        fpCD = os.path.join(clonesDetails,pCD)
        placedClones = pd.read_csv(fpCD, sep='\t')
        placedClones['start'] = placedClones['start'].astype(int)
        placedClones['end'] = placedClones['end'].astype(int)
        if placedClones.empty:
            continue
        newcols = placedClones['attributes'].apply(splitAttributes)
        newcols.columns = [x[:x.find('=')] for x in str(placedClones.iloc[0]['attributes']).split(';')]
        placedClones = pd.concat([placedClones,newcols], axis=1)
        placedClones['Library'] = library
        accs = (placedClones['seqid']).unique()
        splitPlacedClones = [placedClones[placedClones['seqid']==x] for x in accs]
        for oneSeqPlacedClones in splitPlacedClones:
            accession = oneSeqPlacedClones.iloc[0]['seqid']
            accPath = os.path.join(clonesSequences,accession+'.fasta')
            if os.path.isfile(accPath) == False:
                continue
            accloc = allaccs.index(accession)
            if len(allchrs[accloc]) == 1:
                record_iter = SeqIO.parse(open(accPath), "fasta")
                currentseq = list(record_iter)[0]
                de = currentseq.description
                curChrom = de[de.find('chromosome ')+11:de.find(',')]
                accnames[accloc] = curChrom
                seq = currentseq.seq
                seqlen = len(seq)
                allchrs[accloc] = np.asarray([0]*seqlen)
            for clone in oneSeqPlacedClones.iterrows():
                s= int(clone[1]['start'])
                e= int(clone[1]['end'])
                rep = [1]*(e-s)
                loc = list(range(s,e))
                allchrs[accloc][loc]=rep

    outputf = os.path.join(cwd,'coverage.csv')
    with open (outputf,'w') as file:
        w = csv.writer(file)
        w.writerow('Accession','Chromosome','BasesCovered','TotalBases')
        for i in range(len(allaccs)):
            if len(allchrs[i]) == 1:
                continue
            w.writerow([allaccs[i],accnames[i],sum(allchrs[i]),len(allchrs[i])])
            
#Get average length
def getAverageLength(include_libraries=True):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    enzfile = os.path.join(cwd,'shortComm.csv')
    
    placedClonesDetailFs = [x for x in os.listdir(clonesDetails) if x[-3:] == 'gff']
    included_libraries = os.path.join(clonesDetails,'includedLibraries.csv')
    if include_libraries == True:
        with open(included_libraries) as incl:
            reader = csv.reader(incl)
            included_list = [row[0] for row in reader]

    cPlacedClonesDetailFs= placedClonesDetailFs.copy()
    if include_libraries == True:
        cPlacedClonesDetailFs=[]
        for pCD in placedClonesDetailFs:
            library = pCD[:pCD.find('.')]
            if library in included_list:
                cPlacedClonesDetailFs.append(pCD)

    averages = []
    for pCD in cPlacedClonesDetailFs:
        library = pCD[:pCD.find('.')]
        fpCD = os.path.join(clonesDetails,pCD)
        placedClones = pd.read_csv(fpCD, sep='\t')
        placedClones['start'] = placedClones['start'].astype(int)
        placedClones['end'] = placedClones['end'].astype(int)
        placedClones['length'] = placedClones['end'] - placedClones['start']
        average = placedClones['length'].mean()
        averages.append([library, average])

    outputf = os.path.join(cwd,'averagelength.csv')
    with open (outputf,'w') as file:
        w = csv.writer(file)
        w.writerows(averages)
        
#Get statistics for sequenced clones
def getSequencedClonesStats(include_libraries=True):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    enzfile = os.path.join(cwd,'shortComm.csv')
    
    libraryDetails = os.path.join(clonesDetails,'library_'+taxid+'.out')
    sequencedClonesDetails = os.path.join(clonesDetails,'clone_acstate_'+taxid+'_onlyfinished.out')
    sequencedClones = pd.read_csv(sequencedClonesDetails, sep='\t')
    sequencedClones['SeqLen'] = sequencedClones['SeqLen'].astype(int)

    included_libraries = os.path.join(clonesDetails,'includedLibraries.csv')
    if include_libraries == True:
        with open(included_libraries) as incl:
            reader = csv.reader(incl)
            included_list = [row[0] for row in reader]

    liblist = sequencedClones['LibAbbr'].unique()
    cliblist = list(set(liblist) & set(included_list))
    splitty = [sequencedClones[sequencedClones['LibAbbr']==x] for x in cliblist]
    comlist = []
    for seC in splitty:
        comlist.append([seC.iloc[0]['LibAbbr'],seC['SeqLen'].mean(),len(seC)])

    outputf = os.path.join(cwd,'sequencedStats.csv')
    with open (outputf,'w') as file:
        w = csv.writer(file)
        w.writerows(comlist)
    
def onlySingleCutters(row):
    desc = row.iloc[:6]
    vals = row.iloc[6:]
    vals = vals.dropna()
    vals = vals.map(lambda x: ast.literal_eval(x))
    lvals = [len(x)==1 for x in vals.to_numpy()]
    shvals = vals[lvals].map(lambda x: int(x[0])+int(desc['Start']))
    return(shvals)

def findPairs(row):
    shortestoverlap = row[-4]
    longestoverlap = row[-3]
    locmaps = row[-2]
    cloneslistlength = row[-1]
    row = row[0]
    index = row[0]
    row = row[1]
    singles = onlySingleCutters(row)
    nx = index
    pairs = pd.DataFrame(columns = ['Name1','Start1','End1','Enzyme1','Site1',
                                    'Name2','Star2t','End2','Enzyme2','Site2'])
    while True:
        nx+=1
        if nx > cloneslistlength - 2:
            break
        testrow = locmaps.iloc[nx]
        if int(row['End']) < int(testrow['Start']):
            break
        testsingles = onlySingleCutters(testrow)
        for site in singles.items():
            nearsites = testsingles[testsingles > (site[1] - longestoverlap)]
            nearsites = nearsites[nearsites < site[1] - shortestoverlap]
            if len(nearsites) == 0:
                break
            for nearsite in nearsites.items():
                pairs.loc[len(pairs.index)] = [row['Name'],row['Start'],row['End'],site[0],site[1],
                                            testrow['Name'],testrow['Start'],testrow['End'],nearsite[0],nearsite[1]]
    return(pairs)    

def makePairs(cpustouse=1,longestoverlap=200,shortestoverlap=20):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    clonesPairs = os.path.join(cwd,'pairs')
    enzfile = os.path.join(cwd,'shortComm.csv')
    
    if os.path.isdir(clonesPairs) == False:
        os.mkdir(clonesPairs)

    libraries = os.listdir(clonesMaps)
    for library in libraries:
        libpath = os.path.join(clonesMaps, library)
        if os.path.isdir(libpath) == False:
            continue
        chroms = os.listdir(libpath)
        for chrom in chroms:
            chrompath = os.path.join(libpath,chrom)
            locmaps = pd.read_csv(chrompath, sep = '\t', low_memory=False)
            locmaps['Start'] = locmaps['Start'].astype(int)
            locmaps = locmaps.sort_values(by=['Start'])
            locmaps = locmaps.reset_index(drop=True)
            cloneslistlength = len(locmaps)
            chromsave = os.path.join(clonesPairs,library)
            if os.path.isdir(chromsave) == False:
                os.mkdir(chromsave)
            savepath = os.path.join(chromsave,chrom[:-4] + '_pairs.csv')
            pairs = pd.DataFrame(columns = ['Name1','Start1','End1','Enzyme1','Site1',
                                        'Name2','Star2t','End2','Enzyme2','Site2'])
            row = zip(locmaps.iterrows(),repeat(shortestoverlap), repeat(longestoverlap), repeat(locmaps),repeat(cloneslistlength))
            with multiprocessing.Pool(cpustouse) as p:
                for result in p.imap_unordered(findPairs, row):
                    pd.concat([pairs,result],ignore_index=True, axis=0)
            pairs.to_csv(savepath, sep='\t')
            
def getSequenceFromName(name):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
     #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    sequencedMaps = os.path.join(clonesMapsBase,'sequenced')
    placedMaps = os.path.join(clonesMapsBase,'placed')
    enzfile = os.path.join(cwd,'shortComm.csv')
    
    #set up important sequenced lists
    sequencedClonesDetails = os.path.join(clonesDetails,'clone_acstate_'+taxid+'.out')
    sequencedClones = pd.read_csv(sequencedClonesDetails, sep='\t')
    
    #set up important placed lists
    placedClonesLibsFiles = [x for x in os.listdir(clonesDetails)]
    ucname = version + '.unique_concordant.gff'
    ucst = [x for x in placedClonesLibsFiles if ((ucname in x) & ('.info' not in x))]
    placedLibs = dict([((x[:x.find('.')]),os.path.join(clonesDetails,x)) for x in ucst])
    
    #Open a shortened enzyme file, without isoschizomers
    with open(enzfile,"r") as rit:
        rite = csv.reader(rit, lineterminator = '\n', delimiter=",")
        for row in rite:
            shortC = row
    global shortComm
    shortComm = rst.RestrictionBatch(shortC)
    
    #Figure out the name
    findLib = name[:name.find('-')]
    seqdClone = sequencedClones[sequencedClones['CloneName'] == name]
    if len(seqdClone) > 0:
        seqdCloneLine = seqdClone.iloc[0]
        seqAcc = seqdCloneLine['Accession']
        accPath = os.path.join(clonesSequences,seqAcc+'.fasta')
        if os.path.isfile(accPath) == False:
            return('clone accession not found')
        record_iter = SeqIO.parse(open(accPath), "fasta")
        record = list(record_iter)[0]
        seq = record.seq
    else:
        placedFile = placedLibs[findLib]
        uccur = pd.read_csv(placedFile, sep='\t')
        firstline = str(uccur.iloc[0]['attributes']).split(';')
        cnames = [x[:x.find('=')] for x in firstline]
        middles = [x.find('=') for x in firstline]
        newcols = uccur['attributes'].apply(splitAttributesWithMids, middles=middles)
        newcols.columns = cnames
        placedClones = pd.concat([uccur,newcols], axis=1)
        placedCloneLine = placedClones[placedClones['Name']==name]
        if len(placedCloneLine) > 0:
            placedCloneLine = placedCloneLine.iloc[0]
        else:
            return('clone not found')
        accPath = os.path.join(clonesSequences,placedCloneLine['seqid']+'.fasta')
        if os.path.isfile(accPath) == False:
            return('clone accession not found')
        record_iter = SeqIO.parse(open(accPath), "fasta")
        record = list(record_iter)[0]
        seq = record.seq[placedCloneLine['start']:placedCloneLine['end']]
    return(seq)

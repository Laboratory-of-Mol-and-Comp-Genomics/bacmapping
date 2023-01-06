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
    #finseqclones = finseqclones.reset_index()
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
        
        
def narrowDownLibraries():
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
def getCuts(arg):
    row = arg[0]
    maxcuts = arg[1]
    accPath = os.path.join(clonesFastas,row[4]+'.fasta')
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
    
def mapSequencedClones(includeLibraries=True, cpustouse=1, maxcuts=50, chunksize=500):
     #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMaps = os.path.join(cwd,'maps')
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
    knames = ['Name','Start','End','Accession'] + list(shortComm)
    
    sequencedClones = pd.read_csv(sequencedClonesDetails, sep='\t')
    sequencedClones['start'] = '0'
    if include_libraries == True:
        sequencedClones = sequencedClones[sequencedClones['LibAbbr'].isin(included_list)]
    row = zip(sequencedClones['CloneName'], sequencedClones['Chrom'], sequencedClones['start'],
              sequencedClones['SeqLen'], sequencedClones['Accession'], sequencedClones['LibAbbr'])
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
    
def mapPlacedClones(includeLibraries=True, cpustouse=1, maxcuts=50, chunksize=500):
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMaps = os.path.join(cwd,'maps')
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
    knames = ['Name','Start','End','Accession'] + list(shortComm)
    
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
            newcols = placedClones['attributes'].apply(splitAttributesWithMids, args=(middles))
            newcols.columns = cnames
            placedClones = pd.concat([placedClones,newcols], axis=1)
            placedClones['Library'] = library
            accs = (placedClones['seqid']).unique()
            splitPlacedClones = [placedClones[placedClones['seqid']==x] for x in accs] #split clones into separate dfs by accession
            for oneSeqPlacedClones in splitPlacedClones:
                accession = oneSeqPlacedClones.iloc[0]['seqid']
                accPath = os.path.join(clonesFastas,accession+'.fasta')
                if os.path.isfile(accPath) == False:
                    continue
                record_iter = SeqIO.parse(open(accPath), "fasta")
                currentseq = list(record_iter)[0]
                de = currentseq.description
                curChrom = de[de.find('chromosome ')+11:de.find(',')]
                row = zip(oneSeqPlacedClones['Name'],oneSeqPlacedClones['start'],
                          oneSeqPlacedClones['end'],oneSeqPlacedClones['seqid'])
                with multiprocessing.Pool(cpustouse) as p:
                    for result in p.imap_unordered(getCuts, (row,50)):
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
import os
import pandas as pd
import csv
import numpy as np
from ftplib import FTP
from Bio import SeqIO
from Bio import Entrez
from Bio import Restriction as rst
from multiprocessing import Pool
from multiprocessing import Manager
import ast
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
from itertools import repeat
from shutil import rmtree

global cseq
global shortComm
shortC = ['AanI', 'AarI', 'AbaSI', 'AbsI', 'Acc36I', 'AccB7I', 
             'AccII', 'AccIII', 'AclI', 'AclWI', 'AcuI', 'AdeI', 
             'AfiI', 'AflIII', 'AgsI', 'AjiI', 'AjuI', 'AloI', 
             'AluBI', 'Alw21I', 'Alw26I', 'Aor51HI', 'ApaLI', 
             'ArsI', 'AscI', 'Asp718I', 'AsuHPI', 'AvaI', 'AxyI', 
             'BaeGI', 'BaeI', 'BamHI', 'BanII', 'BarI', 'BbvCI', 
             'BccI', 'BceAI', 'BcgI', 'BcnI', 'BglI', 'BglII', 
             'BlnI', 'BlpI', 'Bme18I', 'BmeRI', 'BmiI', 'BmsI', 
             'BmuI', 'BplI', 'BpmI', 'Bpu10I', 'BpuEI', 'BsaJI', 
             'BsaWI', 'BsaXI', 'Bse118I', 'Bse1I', 'Bse8I', 'BseCI', 
             'BseMI', 'BseMII', 'BseRI', 'BseXI', 'BsgI', 'Bsh1285I', 
             'BslFI', 'BsmBI', 'Bsp120I', 'BspACI', 'BspT107I', 
             'BspTNI', 'BsrBI', 'BssHII', 'BssMI', 'Bst1107I', 
             'Bst2BI', 'Bst4CI', 'BstAFI', 'BstAPI', 'BstAUI', 
             'BstC8I', 'BstDEI', 'BstDSI', 'BstPAI', 'BstSNI', 
             'BstV2I', 'BstXI', 'BstYI', 'BstZI', 'BsuI', 'BsuRI', 
             'BtgZI', 'BtsI', 'BtsIMutI', 'BtuMI', 'CaiI', 'CciI', 
             'Cfr13I', 'CpoI', 'CseI', 'CspAI', 'CspCI', 'CviJI', 
             'DraI', 'DrdI', 'EaeI', 'Eam1104I', 'EciI', 'Ecl136II', 
             'Eco147I', 'Eco32I', 'EcoO109I', 'EcoO65I', 'EcoRI', 
             'EcoT22I', 'EheI', 'FaiI', 'FalI', 'FauI', 'FauNDI', 
             'FbaI', 'FblI', 'FokI', 'FspAI', 'FspBI', 'FspEI', 
             'FspI', 'GluI', 'HaeII', 'HapII', 'HhaI', 'HindII', 
             'HindIII', 'HinfI', 'HpaI', 'Hpy166II', 'Hpy188I', 
             'Hpy188III', 'Hpy99I', 'HpyAV', 'HpyCH4V', 'HpyF10VI', 
             'Hsp92I', 'Hsp92II', 'KflI', 'LguI', 'LmnI', 'LpnPI', 
             'MabI', 'MaeIII', 'MauBI', 'MboII', 'MfeI', 'MhlI', 
             'MluCI', 'MluI', 'MmeI', 'MnlI', 'MroNI', 'Msp20I', 
             'MspA1I', 'MspJI', 'MssI', 'MteI', 'Mva1269I', 'NcoI', 
             'NheI', 'NmeAIII', 'NotI', 'NspI', 'NspV', 'OliI', 
             'PacI', 'PaeI', 'PasI', 'PciI', 'PcsI', 'PflFI', 'PfoI', 
             'Ple19I', 'PmaCI', 'Ppu21I', 'PshBI', 'PspFI', 'PspGI', 
             'PspLI', 'PspPPI', 'PspXI', 'PsrI', 'PstI', 'PvuII', 
             'RgaI', 'RigI', 'RsaNI', 'SacII', 'SalI', 'SchI', 
             'ScrFI', 'SetI', 'SfcI', 'SfiI', 'Sfr274I', 'SgeI', 
             'SgrAI', 'SgrDI', 'SmiI', 'SmiMI', 'SmlI', 'SpeI', 
             'SrfI', 'Sse8387I', 'SspI', 'StyI', 'TaiI', 'TaqI', 
             'TaqII', 'TatI', 'TauI', 'TfiI', 'Tru9I', 'TseI', 
             'Tsp45I', 'TspDTI', 'TspGWI', 'TspRI', 'XagI', 'XapI', 
             'XbaI', 'XcmI', 'XmaI', 'XmnI', 'ZraI', 'ZrmI']
shortComm = rst.RestrictionBatch(shortC)


#Download files from the FTP server
def getNewClones(download = True, onlyType = True, vtype = 'BAC', chunk_size=5000, email = 'user@github.com/ewinden/bacmapping'):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    species = 'Homo_sapiens'
    
    #Setup names
    cloneacstate = 'clone_acstate_'+taxid+'.out'
    librarys = 'library_'+taxid+'.out'
    ucname = version + '.unique_'  

    #Setup folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesDetailsRepaired = os.path.join(clonesDetails,'repaired')
    clonesDetailsReordered = os.path.join(clonesDetails,'reordered')
    clonesDetailsInfo = os.path.join(clonesDetails,'info')

    for f in [clonesDetails, clonesSequences, clonesDetailsRepaired,clonesDetailsReordered, clonesDetailsInfo]:
        os.makedirs(f, exist_ok=True)

    #Login to NCBI FTP to download details files
    libraryDetails = os.path.join(clonesDetails,librarys)
    cloneacstatepath = os.path.join(clonesDetails,cloneacstate)
    with FTP("ftp.ncbi.nih.gov") as ftp:
        ftp.login()
        ftp.cwd('repository/clone/reports/'+species+'/')
        for f, p in [[librarys, libraryDetails], [cloneacstate, cloneacstatepath]]:
            if os.path.exists(p) == False:
                ftp.retrbinary("RETR " + f ,open(p, 'wb').write)
        filenames = ftp.nlst()
        ucst = [x for x in filenames if ucname in x]

        #Narrow down
        if onlyType==True:
            librariesToInclude = os.path.join(clonesDetails, 'includedLibraries.csv')
            fullp = pd.read_csv(libraryDetails, sep='\t')
            cut = fullp[fullp['vtype'] == vtype]
            uselibs = cut['libabbr']
            uselibs.to_csv(librariesToInclude, index = False, header = False)
            ucst = [x for x in ucst if x[:x.find('.')] in set(uselibs)]

        for f in ucst:
            p = os.path.join(clonesDetails,f)
            if os.path.exists(p) == False:
                ftp.retrbinary("RETR " + f ,open(p, 'wb').write)

    ucst = os.listdir(clonesDetails)
                
    #Get accession list
    ucstpaths= [os.path.join(clonesDetails,x) for x in ucst]
    ucstlibs = list(set([uc[:uc.find('.')] for uc in ucst if ('unique' in uc)]))

    #Split unique files into a header and details, save fixed files and by accession
    nucpaths = []
    fcnames= ['seqid','source','type','start','end','score','strand','phase','attributes']
    for lib in ucstlibs:
        ucs = [x for x in ucstpaths if lib+'.' in x]
        ucu = ucs[0][:ucs[0].find('unique')]+'unique.gff'
        with open(ucu, 'w') as ci:
            ci.writelines(['seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n'])
        ucurepaired = os.path.join(clonesDetailsRepaired,lib+'_repaired.gff')
        nucpaths.append(ucu)
        for uc in ucs:
            ucinfo = os.path.join(clonesDetailsInfo, uc[uc.find(lib):] + '.info.txt')
            stahp = False
            for header in pd.read_csv(uc, sep='\t', header=None, chunksize=7):
                header.to_csv(ucinfo)
                break
            for tlines in pd.read_csv(uc, sep='\t', skiprows=7, chunksize=1, names = fcnames):
                if tlines.empty:
                    stahp = True
                    break
                atts = tlines['attributes'].item().split(';')
                cnames = [x[:x.find('=')] for x in atts]
                middles = [x.find('=') for x in atts]
                newcols = tlines['attributes'].apply(splitAttributesWithMids, middles=middles)
                newcols.columns = cnames
                flinerep = pd.concat([tlines,newcols], axis = 1)
                flinerep['Library'] = lib
                tlines.to_csv(ucu, mode='a', index = False, header = False, sep = '\t')
                if os.path.exists(ucurepaired) == False:
                    flinerep.to_csv(ucurepaired, mode='w', index = False, header = True, sep = '\t')
                else:
                    flinerep.to_csv(ucurepaired, mode='a', index = False, header = False, sep = '\t')
                ucureorder = os.path.join(clonesDetailsReordered,flinerep['seqid'].item())
                if os.path.exists(ucureorder) == False:
                    flinerep.to_csv(ucureorder, mode='w', index = False, header = True, sep = '\t')
                else:
                    flinerep.to_csv(ucureorder, mode='a', index = False, header = False, sep = '\t')
                break
            if stahp == True:
                continue
            for tlines in pd.read_csv(uc, sep='\t', skiprows=8, chunksize=chunk_size, names = fcnames):
                tlines.to_csv(ucu, mode='a', index=False, header=False, sep='\t')
                newcols = tlines['attributes'].apply(splitAttributesWithMids, middles=middles)
                newcols.columns = cnames
                tlinesrep = pd.concat([tlines,newcols], axis=1)
                tlinesrep['Library'] = lib
                tlinesrep.to_csv(ucurepaired, mode='a', index=False, header=False, sep='\t')
                accessions = tlinesrep['seqid'].unique()
                for acc in accessions:
                    tlinesacc = tlinesrep[tlinesrep['seqid'] == acc]
                    ucureorder = os.path.join(clonesDetailsReordered,acc)
                    if os.path.exists(ucureorder) == False:
                        tlinesacc.to_csv(ucureorder, mode='w', index = False, header = True, sep = '\t')
                    else:
                        tlinesacc.to_csv(ucureorder, mode='a', index = False, header = False, sep = '\t')

    #Get accessions of sequenced clones from clone_acstate
    seqdclones = pd.read_csv(cloneacstatepath, sep='\t')
    finseqclones = seqdclones[(seqdclones['Stdn']=='Y') & (seqdclones['CloneState']=='fin')]
    finseqnames = finseqclones['CloneName'].unique()
    finseqaccs = finseqclones['Accession'].unique()
    finseqclones.to_csv(os.path.join(clonesDetails,'clone_acstate_'+taxid+'_onlyfinished.out'),sep='\t',index=False)

    #Get accessions for placed clones
    allplacedaccs = []
    for uc in nucpaths:
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
    if download == True:
        Entrez.email = email  # Always tell NCBI who you are
        for accession in allaccs: #allplacedaccs for only placed
            save = os.path.join(clonesSequences,accession+'.fasta')
            if os.path.isfile(save) == True:
                continue
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

#splits attributes by hard work
def splitAttributes(ser):
        middles = [x.find('=') for x in ser]
        new = str(ser).split(';')
        clipped = [x[middles[i]+1:] for i,x in enumerate(new)]
        df = pd.Series(clipped)
        return(df)        

#Given a row, gets the sequence and cuts it
def openSeqgetCuts(row):
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

#Given a row, cuts at the globally loaded sequence
def getCuts(row):
    clonesSequences = row[7]
    maxcuts = row[6]
    record = cseq
    lengthtest = int(row[2]) - int(row[3])
    if 25000 <= lengthtest <= 350000:
        seq = record.seq[int(row[3]):int(row[2])]
    elif -350000 <= lengthtest <= -25000:
        seq = record.seq[int(row[2]):int(row[3])]
    else:
        return('NoA')
    cuts = shortComm.search(seq)
    for key in cuts:
        val = cuts[key]
        if type(val) == int:
            lengthv = 1
        else:
            lengthv = len(val)
        if lengthv > maxcuts:
            cuts[key] = 'overflow'
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

#make in the index files for easy searching
def makeIndexFiles(loc):
    indset = pd.DataFrame()
    files = os.listdir(loc)
    for csv in files:
        if csv[-3:] != 'csv' or csv == 'index.csv':
            continue
        fullp = pd.read_csv(os.path.join(loc,csv),sep='\t', low_memory=False)
        shorp = fullp[['Name', 'Chrom', 'Accession']]
        indset = pd.concat([indset, shorp], axis=0)
    indset.to_csv(os.path.join(loc,'index.csv'))

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
    knames = ['Name','Library','Chrom','Start','End','Accession'] + list(shortComm)
    
    sequencedClones = pd.read_csv(sequencedClonesDetails, sep='\t')
    sequencedClones['start'] = '0'
    if include_libraries == True:
        sequencedClones = sequencedClones[sequencedClones['LibAbbr'].isin(included_list)]
    row = zip(sequencedClones['CloneName'], sequencedClones['Chrom'], sequencedClones['start'],
              sequencedClones['SeqLen'], sequencedClones['Accession'], sequencedClones['LibAbbr'],
             repeat(maxcuts), repeat(clonesSequences))
    p = Pool(cpustouse)
    for result in p.imap_unordered(openSeqgetCuts, row):
        if result == 'NoA':
            continue
        folder = os.path.join(clonesMaps, result['Library'])
        if os.path.isdir(folder) == False:
            os.mkdir(folder)
        file = os.path.join(folder, str(result['Chrom'])+'.csv')
        if os.path.isfile(file) == False:
            with open(file, "a") as fwri:
                writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter='\t', fieldnames = knames)
                writer.writeheader()
        with open(file, "a") as fwri:
            writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter="\t", fieldnames = knames)
            writer.writerow(result)
    p.close()
    p.join()
    for lib in os.listdir(clonesMaps):
        makeIndexFiles(os.path.join(clonesMaps,lib))

#Map the end-sequenced clones based on placement details
def mapPlacedClones(cpustouse=1, maxcuts=50, chunk_size=500):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    global cseq

    
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesAccessions = os.path.join(clonesDetails,'reordered')
    clonesRepaired = os.path.join(clonesDetails,'repaired')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')

    if os.path.isdir(clonesMapsBase) == False:
        os.mkdir(clonesMapsBase)
    if os.path.isdir(clonesMaps) == False:
        os.mkdir(clonesMaps)
            
    #Open a shortened enzyme file, without isoschizomers
    knames = ['Name','Library','Chrom','Start','End','Accession'] + list(shortComm)
    accs = os.listdir(clonesAccessions)

    #Get a list of libraries
    libraries = [x[:x.find('_')] for x in os.listdir(clonesRepaired)]
    for lib in libraries:
        os.makedirs(os.path.join(clonesMaps, lib), exist_ok=True)

    #loop through libaries, loop through chunks clones, process them
    for accession in accs:
        accPath = os.path.join(clonesAccessions, accession)
        seqPath = os.path.join(clonesSequences, accession + '.fasta')
        for placedClones in pd.read_csv(accPath, sep='\t', chunksize=chunk_size):
            record_iter = list(SeqIO.parse(open(seqPath), "fasta"))
            currentseq = record_iter[0]
            de = currentseq.description
            curChrom = de[de.find('chromosome ')+11:de.find(',')]
            cseq = currentseq
            row = zip(placedClones['Name'],repeat(curChrom),placedClones['start'],
                        placedClones['end'],placedClones['seqid'],placedClones['Library'],
                        repeat(maxcuts), repeat(clonesSequences))
            p = Pool(cpustouse)
            for result in p.imap_unordered(getCuts, (row)):
                if result == 'NoA':
                    continue
                file = os.path.join(os.path.join(clonesMaps, result['Library']), curChrom+'.csv')
                if os.path.isfile(file) == False:
                    with open(file, "w") as fwri:
                        writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter="\t", fieldnames = knames)
                        writer.writeheader()
                with open(file, "a") as fwri:
                    writer = csv.DictWriter(fwri, lineterminator = '\n', delimiter="\t", fieldnames = knames)
                    writer.writerow(result)
            p.close()
            p.join()

    for lib in os.listdir(clonesMaps):
        makeIndexFiles(os.path.join(clonesMaps,lib))

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
            if os.path.isfile(chrompath) == False:
                continue
            locmaps = pd.read_csv(chrompath, sep = '\t', low_memory=False)
            chromnumber = len(locmaps)
            libnumber+=chromnumber
        libnumbers.append([library,libnumber])
        totalnumber += libnumber
    libnumbers.append(['total',totalnumber])
    outputf = os.path.join(cwd, 'counts.csv')

    with open (outputf,'w') as file:
        w = csv.writer(file)
        w.writerow(['Library', 'Count'])
        w.writerows(libnumbers)

#Get coverage of the genome
def getCoverage():
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    placedAccessions = os.path.join(clonesSequences,'PlacedAccessions.csv')
    placedClonesReordered = os.path.join(clonesDetails, 'reordered')

    libs = [x for x in os.listdir(clonesMaps)]
    dtitles = ['accession', 'chromosome'] + libs + ['total', 'length']
    libscoverage =  pd.DataFrame(columns = dtitles)
    outputf = os.path.join(cwd,'coverage.csv')
    libscoverage.to_csv(outputf, mode='w', index=False)

    for ni, f in enumerate(os.listdir(placedClonesReordered)):
        seqPath = os.path.join(clonesSequences,f + '.fasta')
        cloPath = os.path.join(placedClonesReordered, f)
        record_iter = SeqIO.parse(open(seqPath), "fasta")
        currentseq = list(record_iter)[0]
        de = currentseq.description
        curChrom = de[de.find('chromosome ')+11:de.find(',')]
        seq = currentseq.seq
        seqlen = len(seq)
        totalcoverage = np.asarray([0]*seqlen)
        sect = pd.read_csv(cloPath, sep='\t')
        sumline = [0]*len(libs)
        for ind, lib in enumerate(libs):
            coverage = np.asarray([0]*seqlen)
            insect = sect[sect['Library'] == lib]
            for n, line in insect.iterrows():
                coverage[int(line['start']):int(line['end'])] = 1
            sumline[ind] = sum(coverage)
            totalcoverage += coverage
        totalcoverage[totalcoverage > 0] = 1
        cover = [f,curChrom] + sumline + [sum(totalcoverage), seqlen]
        covert = {dtitles[i]:cover[i] for i in range(len(dtitles))}
        outline = pd.Series(covert).to_frame(ni).T
        outline.to_csv(outputf, mode='a', index=False, header=False)
            
#Get average length
def getAverageLength():
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    detailsRepaired = os.path.join(clonesDetails,'repaired')
    libs = os.listdir(detailsRepaired)

    averages = []

    for lib in libs:
        libpath = os.path.join(detailsRepaired, lib)
        placedClones = pd.read_csv(libpath, sep='\t')
        placedClones['start'] = placedClones['start'].astype(int)
        placedClones['end'] = placedClones['end'].astype(int)
        placedClones['length'] = placedClones['end'] - placedClones['start']
        placedClones = placedClones[placedClones['length'] < 350000]
        placedClones = placedClones[placedClones['length'] > 25000]
        average = placedClones['length'].mean()
        averages.append([lib, average])

    outputf = os.path.join(cwd,'averagelength.csv')
    with open (outputf,'w') as file:
        w = csv.writer(file)
        w.writerow(['Library', 'Average Length (bp)'])
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
        w.writerow(['Library', 'Average Length (bp)', 'Number of BACs used'])
        w.writerows(comlist)

#return a row given a name
def getRow(name):
    #Set taxid and version (most recent human)
    version = '118'
    taxid = '9606'
    
     #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')

    #set up important sequenced lists
    sequencedClonesDetails = os.path.join(clonesDetails,'clone_acstate_'+taxid+'.out')
    sequencedClones = pd.read_csv(sequencedClonesDetails, sep='\t')
    
    #set up important placed lists
    #placedClonesLibsFiles = [x for x in os.listdir(clonesDetails)]
    #ucname = version + '.unique.gff'
    #ucst = [x for x in placedClonesLibsFiles if ((ucname in x) & ('.info' not in x))]
    #placedLibs = dict([((x[:x.find('.')]),os.path.join(clonesDetails,x)) for x in ucst])

    findLib = name[:name.find('-')]
    seqdClone = sequencedClones[sequencedClones['CloneName'] == name]
    if len(seqdClone) > 0:
        seqdCloneLine = seqdClone.iloc[0]
        row = [seqdCloneLine,'sequenced']
    else:
        #placedFile = placedLibs[findLib]
        placedFile = os.path.join(os.path.join(clonesDetails,'repaired'), findLib + '_repaired.gff')
        placedClones = pd.read_csv(placedFile, sep='\t')
        #firstline = str(uccur.iloc[0]['attributes']).split(';')
        #cnames = [x[:x.find('=')] for x in firstline]
        #middles = [x.find('=') for x in firstline]
        #newcols = uccur['attributes'].apply(splitAttributesWithMids, middles=middles)
        #newcols.columns = cnames
        #placedClones = pd.concat([uccur,newcols], axis=1)
        placedCloneLine = placedClones[placedClones['Name']==name]
        if len(placedCloneLine) > 0:
            if len(placedCloneLine) > 1:
                placedCloneLine = placedCloneLine.iloc[0]
        else:
            raise NameError('clone not found')
        row = [placedCloneLine,'placed']
    return(row)

#get insert sequence from row
def getSequence(row, local):
#Prep folders
    cwd = os.getcwd()
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')

    if local == 'sequenced':
        seqAcc = row['Accession']
        accPath = os.path.join(clonesSequences,seqAcc+'.fasta')
        if os.path.isfile(accPath) == False:
            raise NameError('clone accession not found')
        record_iter = SeqIO.parse(open(accPath), "fasta")
        record = list(record_iter)[0]
        seq = record.seq
    else:
        accPath = os.path.join(clonesSequences,row['seqid']+'.fasta')
        if os.path.isfile(accPath) == False:
            raise NameError('clone accession not found')
        record_iter = SeqIO.parse(open(accPath), "fasta")
        record = list(record_iter)[0]
        seq = record.seq[row['start']:row['end']]
    return(seq)

#get single enzyme, checking for isoschizomers
def getRightIsoschizomer(enzyme):
    renzyme = eval('rst.'+enzyme)
    isoschizomers = renzyme.isoschizomers()
    if enzyme not in shortC:
        for en in isoschizomers:
            if str(en) in shortC:
                enzyme = str(en)
                renzyme = en
    return(enzyme,renzyme)

#returns only the single cutting enzymes for a specific BAC row 
def onlySingleCutters(row):
    desc = row.iloc[:6]
    vals = row.iloc[6:]
    vals = vals.where(vals.str[0] == '[').dropna()
    vals = vals.map(lambda x: ast.literal_eval(x))
    lvals = [len(x)==1 for x in vals.to_numpy()]
    shvals = vals[lvals].map(lambda x: int(x[0])+int(desc['Start']))
    return(shvals)

#given a line + data, finds all the pairs for that BAC
def findPairs(cline):
    pairs = pd.DataFrame(columns = ['Name1','Start1','End1','Enzyme1','Site1',
                                    'Name2','Star2t','End2','Enzyme2','Site2'])
    if len(cline[0]) < 2:
        return(pairs)
    
    shortestoverlap = cline[1]
    longestoverlap = cline[2]
    row = cline[0].iloc[0]
    tests = cline[0].iloc[1:]
    singles = onlySingleCutters(row)
    if len(singles) < 1:
        return(pairs)

    for i in range(len(tests)):
        testrow = tests.iloc[i]
        testsingles = onlySingleCutters(testrow)
        if len(testsingles) < 1:
            continue
        for site in singles.items():
            nearsites = testsingles[testsingles > (site[1] - longestoverlap)]
            nearsites = nearsites[nearsites < site[1] - shortestoverlap]
            if len(nearsites) == 0:
                break
            for nearsite in nearsites.items():
                pairs.loc[len(pairs.index)] = [row['Name'],row['Start'],row['End'],site[0],site[1],
                                            testrow['Name'],testrow['Start'],testrow['End'],nearsite[0],nearsite[1]]
    return(pairs)

#make all the possible pairs of BACs with overlap defined by two enzymes that linearize each BAC
def makePairs(cpustouse=1,longestoverlap=200,shortestoverlap=20):
    #Prep folders
    cwd = os.getcwd()
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')
    clonesPairs = os.path.join(cwd,'pairs')
    
    
    os.makedirs(clonesPairs, exist_ok=True)

    allchrs = []
    libraries = os.listdir(clonesMaps)
    libpaths = [os.path.join(clonesMaps, x) for x in libraries]
    for libpath in libpaths:
        for chr in os.listdir(libpath):
            if chr not in allchrs:
                allchrs.append(chr)
    allchrs.remove('index.csv')

    for chr in allchrs:
        #print(chr)
        chrset = [os.path.join(x, chr) for x in libpaths if os.path.isfile(os.path.join(x, chr))]
        locmaps = pd.DataFrame()
        for chrpath in chrset:
            clones = pd.read_csv(chrpath, sep = '\t')
            locmaps = pd.concat([locmaps, clones])
        locmaps['Start'] = locmaps['Start'].astype(int)
        locmaps = locmaps.sort_values(by=['Start'])
        locmaps = locmaps.reset_index(drop=True)
        tozip = []
        for i in range(len(locmaps)-1):
            line = locmaps.iloc[i]
            test = locmaps.iloc[i:]
            test = test[test['Start'] < int(line['End'])]
            if len(test) > 1:
                tozip.append(test)
        savepath = os.path.join(clonesPairs,chr)
        pairs = pd.DataFrame(columns = ['Name1','Start1','End1','Enzyme1','Site1',
                                        'Name2','Star2t','End2','Enzyme2','Site2'])
        zipped = zip(tozip, repeat(shortestoverlap), repeat(longestoverlap))
        p = Pool(cpustouse)
        for result in p.imap_unordered(findPairs, zipped):
            pairs = pd.concat([pairs,result],ignore_index=True)
        p.close()
        p.join()
        pairs.to_csv(savepath, sep='\t')

#get all maps from a clone
def getMaps(name):
    #Prep folders
    cwd = os.getcwd()
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    sequencedMaps = os.path.join(clonesMapsBase,'sequenced')
    placedMaps = os.path.join(clonesMapsBase,'placed')

    #find clone in sequenced maps
    found = False
    findLib = name[:name.find('-')]
    libPath = os.path.join(sequencedMaps,findLib)
    if os.path.isdir(libPath) == True:
        seqind = pd.read_csv(os.path.join(libPath,'index.csv'), low_memory=False)
        cloneLine = seqind[seqind['Name'] == name]
        if len(cloneLine) > 0:
            found = True
    if found == False:
        libPath = os.path.join(placedMaps,findLib)
        if os.path.isdir(libPath) == False:
            raise NameError('clone library not found')
        plaind = pd.read_csv(os.path.join(libPath,'index.csv'), low_memory=False)
        cloneLine = plaind[plaind['Name'] == name]
        if len(cloneLine) > 0:
            found = True
    if found == False:
        raise NameError('clone not found')
    if len(cloneLine) > 1:
        cloneLine = cloneLine.iloc[0]
    chrom = str(cloneLine['Chrom'].item())
    chrompath = os.path.join(libPath,chrom+'.csv')
    mapsList = pd.read_csv(chrompath, sep='\t', low_memory=False)
    rmaps = mapsList[mapsList['Name'] == name]
    trmaps = rmaps.applymap(lambda x: x[1:-1].split(',') if str(x)[0] == '[' else x)
    return(trmaps)

#return restriction map for single enzyme
def getRestrictionMap(name,enzyme):
    maps = getMaps(name)
    nenzyme, r = getRightIsoschizomer(enzyme)
    return(maps[nenzyme].item())

#get insert sequence of a BAC from name
def getSequenceFromName(name):
    row, local = getRow(name)
    sequence = getSequence(row,local)
    return(sequence)

#get sequence between start and end on a chromosome
def getSequenceFromLoc(chrom,start,end):
    #Prep folders
    cwd = os.getcwd()
    clonesSequences = os.path.join(cwd,'sequences')
    
    #Sequence files
    accessions = {}
    for n in range(1,10):
        accessions[str(n)] = 'NC_00000'+str(n)+'.'
    for n in range(10,23):
        accessions[str(n)] = 'NC_0000'+str(n)+'.'
    accessions['X'] = 'NC_000023.'
    accessions['Y'] = 'NC_000024.'
    accession = accessions[str(chrom)]
    accession = [x for x in os.listdir(clonesSequences) if accession in x]
    if len(accession) == 0:
        raise NameError('clone accession not found')
    accession = accession[-1]
    accPath = os.path.join(clonesSequences,accession)
    record_iter = SeqIO.parse(open(accPath), "fasta")
    record = list(record_iter)[0]
    seq = record.seq[start:end]
    return(seq)

#get all the maps between a start and end on a chromosome
def getMapsFromLoc(chrom,start,end, inclusive=True):
    #Prep folders
    cwd = os.getcwd()
    clonesMapsBase = os.path.join(cwd,'maps')
    placedMaps = os.path.join(clonesMapsBase,'placed')
    
    #Get sequenced maps
    maps = pd.DataFrame()
    for lib in os.listdir(placedMaps):
        mapPath = os.path.join(os.path.join(placedMaps,lib),str(chrom)+'.csv')
        if os.path.isfile(mapPath) == False:
            continue
        mapsList = pd.read_csv(mapPath, sep='\t')
        if inclusive == True:
            newmaps = mapsList[mapsList['Start'].astype(int) < int(end)] 
            newmaps = newmaps[newmaps['End'].astype(int) > int(start)]
        else:
            newmaps = mapsList[mapsList['Start'].astype(int) < int(start)]
            newmaps = newmaps[mapsList['End'].astype(int) > int(end)]
        maps = pd.concat([maps, newmaps], ignore_index = True)
    return(maps)

#Draw a circular or linear map of an enzyme
def drawMap(name,enzyme,circular = False):
    maps = getMaps(name)
    nenzyme, r = getRightIsoschizomer(enzyme)
    rmap = maps[nenzyme].item()
    if type(rmap) == str:
        raise NameError('too many cuts')
    nums = [int(x) for x in rmap]
    totallength = np.abs(int(maps['End']) - int(maps['Start']))
    enz = eval('rst.'+enzyme)
    figure, axes = plt.subplots()
    axes.set_aspect( 1 )
    if circular == False:
        plt.plot([0.005, .995],[0.5,0.5],color='dimgray',linestyle='-',linewidth=6)
        plt.text(-.05,0.55,'0')
        plt.text(1.03,.55,str(totallength))
        plt.axis('off')
        plt.plot([0, 0],[0.495,0.575],color='dimgray',linestyle='-',linewidth=3, dash_capstyle='round')
        plt.plot([1, 1],[0.495,0.575],color='dimgray',linestyle='-',linewidth=3, dash_capstyle='round')
        for v in range(0,totallength,2000):
            piv = v/totallength
            if v%10000 == 0:
                plt.plot([piv, piv], [0.5, .55], color='dimgray', linestyle='-', linewidth=1, dash_capstyle='round')
            else:
                plt.plot([piv, piv], [0.5, .525], color='dimgray', linestyle='-', linewidth=1, dash_capstyle='round')
        dropv = 0
        lastv = -1
        cwi = .02
        for v in nums:
            piv = v/totallength
            if (piv - (len(str(v))*cwi)) > lastv:
                dropv = 0
            else:
                dropv += 0.03
            if enz.is_blunt() == True:
                plt.plot([piv, piv], [0.45-dropv, .54], color='red', linestyle='-', linewidth=1, zorder=5)
            elif enz.is_3overhang:
                plt.plot([piv, piv], [0.45-dropv, .5], color='red', linestyle='-', linewidth=1, zorder=5)
                plt.plot([piv, piv+0.005], [.5, .5], color='red', linestyle='-', linewidth=1, zorder=5)
                plt.plot([piv+0.005, piv+0.005], [.5,.54], color='red', linestyle='-', linewidth=1, zorder=5)
            else:
                plt.plot([piv, piv], [0.45-dropv, .5], color='red', linestyle='-', linewidth=1, zorder=5)
                plt.plot([piv, piv-0.005], [.5, .5], color='red', linestyle='-', linewidth=1, zorder=5)
                plt.plot([piv-0.005, piv-0.005], [.5,.54], color='red', linestyle='-', linewidth=1, zorder=5)
            plt.text(piv+0.005,(0.45-dropv),str(v), zorder=10)
            rect = patches.Rectangle((piv+0.005,(0.45-dropv)), (len(str(v))*cwi), 0.03, facecolor='white', zorder=7)
            axes.add_patch(rect)
            lastv = piv
    else:
        cc = plt.Circle(( 0 , 0 ), 1 , fill= False) 
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

#get all the possible pairs of a specific bac by name
def findPairsFromName(name, longestoverlap, shortestoverlap):
    cwd = os.getcwd()
    clonesSequences = os.path.join(cwd,'sequences')
    clonesMapsBase = os.path.join(cwd,'maps')
    clonesMaps = os.path.join(clonesMapsBase,'placed')

    library = name[:name.find('-')]
    libpath = os.path.join(clonesMaps, library)
    indf = pd.read_csv(os.path.join(libpath,'index.csv'))
    placedCloneLine = indf[indf['Name'] == name]
    if len(placedCloneLine) > 0:
        placedCloneLine = placedCloneLine.iloc[0]
    else:
        raise NameError('clone not found')
    #ind = placedCloneLine.iloc[0]
    #accession = placedCloneLine['Accession']
    chrom = placedCloneLine['Chrom']
    chrompath = os.path.join(libpath,chrom+'.csv')
    locmaps = pd.read_csv(chrompath, sep = '\t', low_memory=False)
    locmaps['Start'] = locmaps['Start'].astype(int)
    locmaps['End'] = locmaps['End'].astype(int)
    locmaps = locmaps.sort_values(by=['Start'])
    locmaps = locmaps.reset_index(drop=True)
    line = locmaps[locmaps['Name'] == name]
    if len(line) > 0:
        if len(line) > 1:
            line = line.iloc[0]
        ind = np.where(locmaps['Name']==name)
    else:
        raise NameError('clone not found')
    slocmaps = pd.concat([line,locmaps[locmaps['Start'] > line['Start'] & locmaps['Start'] < line['End']]])
    send = [slocmaps,shortestoverlap,longestoverlap]
    pairs = findPairs(send)
    return(pairs)

def findOverlappingBACs(name):
    #Prep folders
    cwd = os.getcwd()
    clonesDetails = os.path.join(cwd,'details')
    clonesSequences = os.path.join(cwd,'sequences')
    clonesAccessions = os.path.join(clonesDetails,'reordered')
    row, local = getRow(name)
    acc = row['seqid'].item()
    start = row['start'].item()
    end = row['end'].item()
    bacsfile = os.path.join(clonesAccessions,acc)
    bacset = pd.read_csv(bacsfile, sep='\t')
    subset = bacset[bacset['start'] > start]
    subset = subset[subset['start'] < end]
    subset['overlaplength'] = (end - subset['start'])
    return(subset)
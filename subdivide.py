import sys
import gzip
import time

# Be exclusive in what we copy
chr1Range=chr2Range=chr3Range=chr4Range=chr5Range='0-0'
batch = ''
# Uses Batch files to specify which sequences we will copy a section of.
# More sensible than file by file, since expected use is to perform this script on
# a number of files, a population of interest, to compare them in specific places.

startTime = time.time()
index = 1
while ( index <= len(sys.argv)-1 ):
    arg = sys.argv[index]

#------------- chromosome range evaluation -------------#
    if (arg == 'chr' and index+2 <= len(sys.argv)-1): # avoid IndexOutOfBounds
        try:
            c = int(sys.argv[index+1])
            assert c > 0 and c <= 5
            r = sys.argv[index+2]
            # test validity of range
            test = r.split('-')
            testMin = int(test[0])
            testMax = int(test[1])
            if c == 1:
                chr1Range = r
            elif c == 2:
                chr2Range = r
            elif c == 3:
                chr3Range = r
            elif c == 4:
                chr4Range = r
            elif c == 5:
                chr5Range = r
            index = index+2
        except ValueError:
            print('Error: \'chr\' must be followed by the desired chromosome number [1-5], then by a range #-###.')
            sys.exit()

#------------- Consider input to be a methylome sequence file we want to compare -------------#
    else:
        try:
            with open('./Batch/' + arg, 'r') as test:
                batch = arg
            # saves given filename for later. This remembers only the last valid file.
        except FileNotFoundError:
            print('Given filename \'' + arg + '\' was not found under Batch directory.')
            sys.exit()

    index = index+1 # end of loop

#------------- describes the arguments of this operation -------------#
description = ''
if chr1Range != '0-0':
    description += ('_chr1_'+chr1Range)
if chr2Range != '0-0':
    description += ('_chr2_'+chr2Range)
if chr3Range != '0-0':
    description += ('_chr3_'+chr3Range)
if chr4Range != '0-0':
    description += ('_chr4_'+chr4Range)
if chr5Range != '0-0':
    description += ('_chr5_'+chr5Range)

#------------- Make new Batch to hold names of all newly created Data files  -------------#
newBatch = open('./Batch/subseq_from_'+batch+'_on'+description,'w')

#------------- for every given Data file -------------#
with open('./Batch/'+batch, 'r') as src:
    for line in src:
        seqFile = line.strip('\n').split(',')[0]

#------------- result sequence filename -------------#
        name = seqFile.split('.')[0]+'_SUBSEQ'+description+'.tsv.gz'

        subSeq = gzip.open('./Data/'+name, 'w')
        with gzip.open('./Data/'+seqFile, 'r') as data:
            for line in data:
                cols = line.decode('utf-8').split('\t')

                # check validity of chromosome and range
                if cols[0] == '1' and int( cols[1] ) >= int( chr1Range.split('-')[0] ) and int( cols[1] ) <= int( chr1Range.split('-')[1] ):
                    subSeq.write(line)
                elif cols[0] == '2' and int( cols[1] ) >= int( chr2Range.split('-')[0] ) and int( cols[1] ) <= int( chr2Range.split('-')[1] ):
                    subSeq.write(line)
                elif cols[0] == '3' and int( cols[1] ) >= int( chr3Range.split('-')[0] ) and int( cols[1] ) <= int( chr3Range.split('-')[1] ):
                    subSeq.write(line)
                elif cols[0] == '4' and int( cols[1] ) >= int( chr4Range.split('-')[0] ) and int( cols[1] ) <= int( chr4Range.split('-')[1] ):
                    subSeq.write(line)
                elif cols[0] == '5' and int( cols[1] ) >= int( chr5Range.split('-')[0] ) and int( cols[1] ) <= int( chr5Range.split('-')[1] ):
                    subSeq.write(line)

        # write newly created file to new Batch file
        newBatch.write(name + '\n')
        subSeq.close()
        print(str(round(time.time() - startTime, 2)) + 'seconds for '+ seqFile)

newBatch.close()

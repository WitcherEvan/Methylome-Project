import os
import sys
import gzip
import re

# Command line arguments:
# None required, no ordering of the three enforced.
# Keep 'chr' and 'range' paired with their # as indicated.
#
# chr [1-5]     tag is chr, followed by an integer 1-5 for chromosome number
# range #-#     range over which we will analyze.
#               Numbers must be positive integers seperated by '-'
# filename      a file containing the filenames of sequences we are interested in.
#               Current plan: have these filename files in seperate subdirectory

# The different columns of Methylome csv data:
# [0] : chromosome number, 1 through 5
# [1] : position in genome
# [2] : strand. plus/minus value, unclear use
# [3] : mc_class, the variety of Methylated Cytosine
# [4] : methylated_bases, number found to be methylated_bases
# [5] : total_bases, total number observed
# [6] : methylation_call, 1 indicates this position considered Methylated

# Important for use: Data is the folder/directory name under which the
# sequence files are kept. Change here or on system for consistency.
# This read.py script should be in the same directory as your Data directory.
dirName = 'Data/'

# default values. Change based on goals and convenience
batch = 'default.txt'
minPos = 0
maxPos = sys.maxint
chromo = 1

# Parsing input section:
# The current setup would keep taking inputs, and save the last valid one.
# Fix? Depends on whether we want option to put many filenames directly into command line

if (len(sys.argv) > 1): # then we have parameters
    index = 1
    while ( index <= len(sys.argv)-1 ):
        arg = sys.argv[index]

        if (arg == 'chr' and index+1 <= len(sys.argv)-1): # avoid IndexOutOfBounds
            try:
                num = int(sys.argv[index+1])
                # this line depends on our knowledge that there are 5 chromosomes
                if num > 0 and num <= 5:
                    chromo = num
                else:
                    print('Positive integer please. Using default chr.')
            except ValueError:
                print('Error: chr must be followed by the chromosome number, a positive integer.')

        elif (arg == 'range' and index+1 <= len(sys.argv)-1): # avoid IndexOutOfBounds
            try:
                rn = sys.argv[index+1]
                rmin = int(rn.split('-')[0])
                rmax = int(rn.split('-')[1])
                if rmin >= 0 and rmax > 0:
                    minPos = rmin
                    maxPos = rmax
                else:
                    print('Positive integers please. Using default range.')
            except ValueError:
                print('Error: range must be followed by min-max position on the genome.')

        else:
            try:
                with open('./Batch/' + arg, 'r') as test:
                    batch = arg
            except FileNotFoundError:
                print('Given filename \'' + arg + '\' was not found under Batch directory.')

        index = index+1


# start reading sequences from designated batch-file
with open('./Batch/' + batch, 'r') as seqList:
    for filename in seqList:
        path = dirName + filename.strip('\n')

        # Every file
        with gzip.open(path, 'r') as myzip:
            next(myzip) # skip title headers
            methylCount = 0
            n = 0

            # Every line of csv
            for byteForm in myzip:
                row = byteForm.decode('utf-8').strip('\n')
                columns = row.split('\t')

                if (int(columns[1]) <= minPos): # if current position less that target,
                    continue # not at start of target range

                if (int(columns[1]) <= maxPos): # we are in desired range. Do things

                    # here we do the processing, add info to dictionaries, etc
                    n = n+1
                    methylCount = methylCount + int(columns[6])

                else: # we have completed target range. end of code section
                    print(path + '\t' + str(methylCount) + ' methlyated out of ' + str(n))

                    break # ends reading of myzip file, move to next filesOfInterest

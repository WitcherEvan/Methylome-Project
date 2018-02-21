import os
import sys
import gzip
import re

# The different columns of Methylome csv data:
# [0] : chromosome number, 1 through 5
# [1] : position in genome
# [2] : strand. plus/minus value, unclear use
# [3] : mc_class, the variety of Methylated Cytosine
# [4] : methylated_bases, number found to be methylated_bases
# [5] : total_bases, total number observed
# [6] : methylation_call, 1 indicates this position considered Methylated

# Important for use: Data is the folder/directory name under which the
# sequence files are kept. Change here or on system for consistenty.
# This script should be in the same directory as your Data directory.
dirName = 'Data/'

# default values
filesOfInterest = os.listdir(dirName)
minPos = 0
maxPos = 100000
chromo = 1

# Parsing input section:


# start reading sequences
for el in filesOfInterest:
    path = dirName + el

    with gzip.open(path, 'r') as myzip:
        next(myzip) # skip title headers
        n = 0

        for byteForm in myzip:
            row = byteForm.decode('utf-8').strip('\n')
            columns = row.split('\t')

            if (int(columns[1]) >= minPos):
                continue # not at start of target range

            if (int(columns[1]) > maxPos):
                break # we have completed target range

            # here we do the processing, add info to dictionaries, etc
                # print(row)
                n = n+1

        print('lines processed: '+str(n)) # for our information, maybe useful

import os
import sys
import gzip
import re
import time

start = time.time()
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
# This summarize.py script should be in the same directory as your Data directory.
dirName = 'Data/'

# default values. Change based on goals and convenience
batch = 'default.txt'
minPos = 0
maxPos = 31000000 # 31 million covers the longerst chromosome.
chromo = 'all' # will consider all chromosomes by default
defaultRange = 1

# Parsing input section:
# The current setup would keep taking inputs, and save the last valid one.

if (len(sys.argv) > 1): # then we have parameters
    # Sort out the given parameters.
    index = 1
    while ( index <= len(sys.argv)-1 ):
        arg = sys.argv[index]

#------------- chromosome number evaluation -------------#
        if (arg == 'chr' and index+1 <= len(sys.argv)-1): # avoid IndexOutOfBounds
            try:
                num = int(sys.argv[index+1])
                # this line depends on our knowledge that there are 5 chromosomes
                if num > 0 and num <= 5:
                    chromo = num
                else:
                    print('Please use valid chromosome number, 1 to 5. Using default value.')
                # Increment index past the value given, since the lack of valueError suggests the input value
                # was intended as a chromosome number, but was invalid.
                index = index+1
            except ValueError:
                print('Error: \'chr\' must be followed by the desired chromosome number, a positive integer.')
                sys.exit()

#------------- Base pair position evaluation -------------#
        elif (arg == 'range' and index+1 <= len(sys.argv)-1): # avoid IndexOutOfBounds
            try:
                defaultRange = 0
                rn = sys.argv[index+1]
                rmin = int(rn.split('-')[0])
                rmax = int(rn.split('-')[1])
                if rmin >= 0 and rmax > 0:
                    minPos = rmin
                    maxPos = rmax
                else:
                    print('Positive integers please. Using default range.')
                index = index+1
            except ValueError:
                print('Error: \'range\' must be followed by dash-seperated min-max base pair positions.')
                sys.exit()

#------------- Consider input to be the batch-file, specifying the methylomes we want to work on -------------#
        else:
            try:
                with open('./Batch/' + arg, 'r') as test:
                    batch = arg
                # saves given filename for later. This remembers only the last valid file.
            except FileNotFoundError:
                print('Given filename \'' + arg + '\' was not found under Batch directory.')
                sys.exit()

# Increment index at end of loop
        index = index+1

# ------ Prep the result file ------
if defaultRange == 1:
    description = './Results/Summary_chr'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
else:
    description = './Results/Summary_chr'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'
try:
    result = open(description, 'w')
    result.write('sequence,overall,CG_context,CHH_context,CHG_context\n')
    # add others if necessary

except FileNotFoundError: # When dir Results does not exist
    os.mkdir('./Results')
    result = open(description, 'w')
    result.write('ecotype,chromosome,start-pos,end-pos,overall,CG_context,CHH_context,CHG_context\n')

# start reading sequences from designated batch-file
with open('./Batch/' + batch, 'r') as seqList:
    for filename in seqList:
        path = dirName + filename.strip('\n')
        subStart = time.time()
        # Reset for every new file
        CG_methyl=CG_total=CHH_methyl=CHH_total=CHG_methyl=CHG_total=linecount=all_methyls = 0

        # Seperate branch of opertions for default behaviour. Doesnt feel great,
        # quasi-duplicating my code this way, but it helps keep both functionalities
        # within this one summarize.py script
        if chromo == 'all':
            try:
                with gzip.open(path, 'r') as myzip:
                    next(myzip) # skip title headers

                    # Every line of csv
                    for byteForm in myzip:
                        row = byteForm.decode('utf-8').strip('\n')
                        columns = row.split('\t')

                        # Target range not reached, continue
                        if int(columns[1]) < minPos or int(columns[1]) > maxPos:
                            continue
                        # reached, do stuff
                        all_methyls += int(columns[6])
                        linecount+=1

                        match = re.search(r'CG', columns[3])
                        if match: #columns[3] == 'CG':
                            CG_methyl = CG_methyl + int(columns[6])
                            CG_total += 1

                        match = re.search(r'C[A,T,C]G', columns[3])
                        if match:
                            CHG_methyl = CHG_methyl + int(columns[6])
                            CHG_total += 1

                        match = re.search(r'C[A,T,C][A,T,C]', columns[3])
                        if match:
                            CHH_methyl = CHH_methyl + int(columns[6])
                            CHH_total += 1
                        # End of line-loop. Gets next line in myzip file

                    # on finishing all lines of a file,
                    overall = str(round(all_methyls / linecount, 6))
                    # CG_context = str(round((CG_methyl / CG_total), 6))
                    # CHH_context = str(round((CHH_methyl / CHH_total), 6))
                    # CHG_context = str(round((CHG_methyl / CHG_total), 6))
                    CG_context = str(round((CG_methyl / all_methyls), 6))
                    CHH_context = str(round((CHH_methyl / all_methyls), 6))
                    CHG_context = str(round((CHG_methyl / all_methyls), 6))

                    # write results to csv file to save answers
                    result.write(path+','+str(chromo)+','+str(minPos)+','+str(maxPos)+','+ overall +','+ CG_context +','+ CHH_context +','+ CHG_context +'\n')
                    # timer per file
                    print(str(round(time.time() - subStart, 2)) + 'seconds for '+ path)

            except FileNotFoundError:
                print('File ' + path + " was not found.")

        else: # This branch is where a single chromosome is specified.
            try:
                with gzip.open(path, 'r') as myzip:
                    next(myzip) # skip title headers

                    # Every line of csv
                    for byteForm in myzip:
                        row = byteForm.decode('utf-8').strip('\n')
                        columns = row.split('\t')

                        if int(columns[0]) < chromo or int(columns[1]) < minPos: # if current chr/position less that target,
                                continue # with next line in myzip

                        # Here we've passed the checks for finding our target chromosome.
                        if int(columns[0]) == chromo and int(columns[1]) <= maxPos: # Do things
                            all_methyls += int(columns[6])
                            linecount+=1

                            match = re.search(r'CG', columns[3])
                            if match: #columns[3] == 'CG':
                                CG_methyl = CG_methyl + int(columns[6])
                                CG_total += 1

                            match = re.search(r'C[A,T,C]G', columns[3])
                            if match:
                                CHG_methyl = CHG_methyl + int(columns[6])
                                CHG_total += 1

                            match = re.search(r'C[A,T,C][A,T,C]', columns[3])
                            if match:
                                CHH_methyl = CHH_methyl + int(columns[4])
                                CHH_total += 1

                            continue # to next line in myzip file

                        else: # we have completed target chr-range on a sequence
                            overall = str(round(all_methyls / linecount, 6))
                            CG_context = str(round((CG_methyl / CG_total), 6))
                            CHH_context = str(round((CHH_methyl / CHH_total), 6))
                            CHG_context = str(round((CHG_methyl / CHG_total), 6))

                            # write results to csv file to save answers
                            result.write(path+','+str(chromo)+','+str(minPos)+','+str(maxPos)+','+ overall +','+ CG_context +','+ CHH_context +','+ CHG_context +'\n')
                            # timer per file
                            print(str(round(time.time() - subStart, 2)) + 'seconds for '+ path)
                            break # from line-reading from current file

            except FileNotFoundError:
                print('File ' + path + " was not found.")


        # end of loop for a file

result.close()
print(str(round(time.time() - start,2)) + ' seconds to completion.')

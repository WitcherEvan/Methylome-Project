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
# This segSites.py script should be in the same directory as your Data directory.
dirName = 'Data/'

# default values. Change based on goals and convenience
minPos = 0
maxPos = 31000000 # 31 million covers the longerst chromosome.
chromo = -1 # will consider all chromosomes by default

# Parsing input section. Saves the last valid inputs for a category.
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

#------------- Consider input to be a methylome sequence file we want to compare -------------#
        else:
            try:
                with open('./Batch/' + arg, 'r') as test:
                    batch = arg
                # saves given filename for later. This remembers only the last valid file.
            except FileNotFoundError:
                print('Given filename \'' + arg + '\' was not found under Batch directory.')
                sys.exit()

# --------- Increment to read next command-line argument. ---------
        index = index+1


# -------- Confirm validity of given sequences --------
confirmed = []
with open('./Batch/' + batch, 'r') as seqList:
    for el in seqList:
        try:
            el = el.strip('\n')
            f = gzip.open(dirName+el, 'r')
            next(f) # verify this is safe for later, that it is not empty.
            confirmed.append(el)
        except FileNotFoundError:
            print('File ' + dirName + el + " was not found. Excluding.")

if len(confirmed) < 2:
    print('Need minimum 2 valid methylomes to compare segregating sites.')
    sys.exit()


# -------- Prep the result file --------
description = './Results/comparison_chr'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+str(confirmed)+'.txt'
try:
    result = open(description, 'w')
except FileNotFoundError: # When dir Results does not exist
    os.mkdir('./Results')
    result = open(description, 'w')


# The branch for looking at whole methylomes
if chromo == -1:

    # loop to do all n*(n-1)/2 comparisons given n methylomes
    index1 = 0
    while index1 < len(confirmed)-1:
        seq1 = confirmed[index1]
        index2 = index1 +1

        while index2 < len(confirmed):
            subStart = time.time() # to time each sequence to sequence comparison.
            seq2 = confirmed[index2]
            f1 = gzip.open(dirName+seq1, 'r')
            next(f1) # skip headers
            f2 = gzip.open(dirName+seq2, 'r')
            next(f2) # skip headers

            total_m = total_s = CG_m = CG_s = CHH_m = CHH_s = CHG_m = CHG_s = contextMismatch = 0

            # -------- Now we are in the comparison section --------
            try: # This 'try' block allows us to leave a loop on reaching the end of a file.
                line1 = next(f1).decode('utf-8').strip('\n').split('\t') # get next lines
                line2 = next(f2).decode('utf-8').strip('\n').split('\t')

                while 1: # We leave this otherwise infinite loop by hitting an IndexError at the end of a file.
                    if line1[0] == line2[0]: # chromosomes are the same
                        if  line1[1] == line2[1]: # Positions are matched.
                            check = int(line2[6]) + int(line1[6])

# Note: do we want to count a match only when context is the same?
# Here I count them as methy-matches when they have same context.

                            # matching
                            if check == 2: # line1[6] == 1 and line2[6] == 1:
                                if line1[3] == line2[3]:

                                    match = re.search(r'CG', line1[3])
                                    if match:
                                        CG_m = CG_m + 1
                                        total_m = total_m +1

                                    match = re.search(r'C[A,T,C]G', line1[3])
                                    if match:
                                        CHG_m = CHG_m + 1
                                        total_m = total_m +1

                                    match = re.search(r'C[A,T,C][A,T,C]', line1[3])
                                    if match:
                                        CHH_m = CHH_m +1
                                        total_m = total_m +1

                                else:
                                    contextMismatch = contextMismatch +1

                            # segregation
                            elif check == 1:
                                if line1[3] == line2[3]:

                                    match = re.search(r'CG', line1[3])
                                    if match:
                                        CG_s = CG_s + 1
                                        total_s = total_s +1

                                    match = re.search(r'C[A,T,C]G', line1[3])
                                    if match:
                                        CHG_s = CHG_s + 1
                                        total_s = total_s +1

                                    match = re.search(r'C[A,T,C][A,T,C]', line1[3])
                                    if match:
                                        CHH_s = CHH_s +1
                                        total_s = total_s +1
                                else:
                                    contextMismatch = contextMismatch +1

                            # End of position comparison, get next lines.
                            line1 = next(f1).decode('utf-8').strip('\n').split('\t')
                            line2 = next(f2).decode('utf-8').strip('\n').split('\t')

                        # Adjust positions.
                        elif line1[1] > line2[1]:
                            line2 = next(f2).decode('utf-8').strip('\n').split('\t')
                        elif line1[1] < line2[1]:
                            line1 = next(f1).decode('utf-8').strip('\n').split('\t')

                    else: # dealing with unusual files, since the chr don't match.
                        while line1[0] > line2[0]:
                            line2 = next(f2).decode('utf-8').strip('\n').split('\t')
                        while line1[0] < line2[0]:
                            line1 = next(f1).decode('utf-8').strip('\n').split('\t')

            except StopIteration:
                print('End of comparison '+seq1+' & '+seq2)
                print(str(round(time.time() - subStart, 2)) + 'seconds for a comparison.')

            # Write results to file. So many things, feels better to write as a paragraph.
                result.write('Comparing '+str(seq1)+' to '+str(seq2)+ ':\n')
                result.write('All Shared Methylations: '+str(total_m)+', Percent All Shared Methylations: '+str(total_m/(total_m+total_s))+'.\n')
                result.write('All Segregations: '+str(total_s)+', Percent All Segregations: '+str(total_s/(total_s+total_m))+'\n')
                result.write('Mismatched contexts (1 or 2 methylations at a position where mc_class differs): '+str(contextMismatch))

                result.write('CG context:\nCG Shared Methylations: '+str(CG_m)+', Percent Methylation: '+str(CG_m/(CG_m+CG_s))+'.\n')
                result.write('CG Segregations: '+str(CG_s)+', Percent Segregation: '+str(CG_s/(CG_s+CG_m))+'.\n')

                result.write('CHG context:\nCHG Shared Methylations: '+str(CHG_m)+', Percent Methylation: '+str(CHG_m/(CHG_m+CHG_s))+'.\n')
                result.write('CHG Segregations: '+str(CHG_s)+', Percent Segregation: '+str(CHG_s/(CHG_s+CHG_m))+'.\n')

                result.write('CHH context:\nCHH Shared Methylations: '+str(CHH_m)+', Percent Methylation: '+str(CHH_m/(CHH_m+CHH_s))+'.\n')
                result.write('CHH Segregations: '+str(CHH_s)+', Percent Segregation: '+str(CHH_s/(CHH_s+CHH_m))+'.\n\n')

            f1.close()
            f2.close()
            index2 = index2 +1 # end of inner

        index1 = index1 +1 # end of outer

else: # This branch is where a single chromosome is specified.

    # loop to do all n*(n-1)/2 comparisons given n methylomes
    index1 = 0
    while index1 < len(confirmed)-1:
        seq1 = confirmed[index1]
        index2 = index1 +1

        while index2 < len(confirmed):
            subStart = time.time()
            seq2 = confirmed[index2]
            f1 = gzip.open(dirName+seq1, 'r')
            next(f1) # skip headers
            f2 = gzip.open(dirName+seq2, 'r')
            next(f2) # skip headers

            total_m = total_s = CG_m = CG_s = CHH_m = CHH_s = CHG_m = CHG_s = 0

            # -------- Now we are in the comparison part --------
            try:
                line1 = next(f1).decode('utf-8').strip('\n').split('\t') # get next lines
                line2 = next(f2).decode('utf-8').strip('\n').split('\t')

                while 1:
                    # get the desired chromosome. Not assuming both sequences
                    # will reach it at the same time, so we increment seperately.
                    while chromo != line1[0]:
                        line1 = next(f1).decode('utf-8').strip('\n').split('\t') # get next lines
                    while chromo != line2[0]:
                        line2 = next(f2).decode('utf-8').strip('\n').split('\t')

                    while chromo == line1[0] and chromo == line2[0]:
                        if  line1[1] == line2[1]: # positions match. Compare:
                            check = int(line2[6]) + int(line1[6])

                            # matching
                            if check == 2: # line1[6] == 1 and line2[6] == 1:
                                if line1[3] == line2[3]:

                                    match = re.search(r'CG', line1[3])
                                    if match:
                                        CG_m = CG_m + 1
                                        total_m = total_m +1

                                    match = re.search(r'C[A,T,C]G', line1[3])
                                    if match:
                                        CHG_m = CHG_m + 1
                                        total_m = total_m +1

                                    match = re.search(r'C[A,T,C][A,T,C]', line1[3])
                                    if match:
                                        CHH_m = CHH_m +1
                                        total_m = total_m +1

                                else:
                                    contextMismatch = contextMismatch +1

                            # segregation
                            elif check == 1:
                                if line1[3] == line2[3]:

                                    match = re.search(r'CG', line1[3])
                                    if match:
                                        CG_s = CG_s + 1
                                        total_s = total_s +1

                                    match = re.search(r'C[A,T,C]G', line1[3])
                                    if match:
                                        CHG_s = CHG_s + 1
                                        total_s = total_s +1

                                    match = re.search(r'C[A,T,C][A,T,C]', line1[3])
                                    if match:
                                        CHH_s = CHH_s +1
                                        total_s = total_s +1
                                else:
                                    contextMismatch = contextMismatch +1

                            # end of section matched, move forward both.
                            line1 = next(f1).decode('utf-8').strip('\n').split('\t') # get next lines
                            line2 = next(f2).decode('utf-8').strip('\n').split('\t')

                        # adjust positions
                        elif line1[1] > line2[1]:
                            line2 = next(f2).decode('utf-8').strip('\n').split('\t')
                        elif line1[1] < line2[1]:
                            line1 = next(f1).decode('utf-8').strip('\n').split('\t')

                    # the chromosome of interest is complete.
                    print('~~End of comparison '+seq1+' & '+seq2)
                    print(str(round(time.time() - subStart, 2)) + 'seconds for a comparison.')
                    break # from `while 1` loop

            except IndexError:
                print('End of comparison '+seq1+' & '+seq2)
                print(str(round(time.time() - subStart, 2)) + 'seconds for a comparison.')

            # Write results to file. So many things, feels better to write as a paragraph.
            result.write('Comparing '+str(seq1)+' to '+str(seq2)+ ':\n')
            result.write('All Shared Methylations: '+str(total_m)+', Percent All Shared Methylations: '+str(total_m/(total_m+total_s))+'.\n')
            result.write('All Segregations: '+str(total_s)+', Percent All Segregations: '+str(total_s/(total_s+total_m))+'\n')
            result.write('Mismatched contexts (1 or 2 methylations at a position where mc_class differs): '+str(contextMismatch))

            result.write('CG context:\nCG Shared Methylations: '+str(CG_m)+', Percent Methylation: '+str(CG_m/(CG_m+CG_s))+'.\n')
            result.write('CG Segregations: '+str(CG_s)+', Percent Segregation: '+str(CG_s/(CG_s+CG_m))+'.\n')

            result.write('CHG context:\nCHG Shared Methylations: '+str(CHG_m)+', Percent Methylation: '+str(CHG_m/(CHG_m+CHG_s))+'.\n')
            result.write('CHG Segregations: '+str(CHG_s)+', Percent Segregation: '+str(CHG_s/(CHG_s+CHG_m))+'.\n')

            result.write('CHH context:\nCHH Shared Methylations: '+str(CHH_m)+', Percent Methylation: '+str(CHH_m/(CHH_m+CHH_s))+'.\n')
            result.write('CHH Segregations: '+str(CHH_s)+', Percent Segregation: '+str(CHH_s/(CHH_s+CHH_m))+'.\n\n')

            f1.close()
            f2.close()
            index2 = index2 +1 # end of inner

        index1 = index1 +1 # end of outer

result.close()
print(str(round(time.time() - start,2)) + ' seconds to completion.')

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
# HRR           for a Human Readable Report

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
chromo = 'all' # will consider all chromosomes by default
hrr = 0
defaultRange = 1
minBases = ''
methylRatio = ''

# Parsing input section. Saves the last valid inputs for a category.
if (len(sys.argv) > 1): # then we have parameters
    # Sort out the given parameters.
    index = 1
    while ( index <= len(sys.argv)-1 ):
        arg = sys.argv[index]

#------------- Human readable report -------------#
        if arg == 'HRR':
            hrr = 1

#------------- chromosome number evaluation -------------#
        elif (arg == 'chr' and index+1 <= len(sys.argv)-1): # avoid IndexOutOfBounds
            try:
                num = int(sys.argv[index+1])
                # this line depends on our knowledge that there are 5 chromosomes
                if num > 0 and num <= 5:
                    chromo = num
                else:
                    print('Please use valid chromosome number, 1 to 5, or omit \'chr\' to consider all chromosomes.')
                    sys.exit()
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

#------------- Minimum bases requirement -------------#
        elif (arg == 'minB' and index+1 <= len(sys.argv)-1): # avoid IndexOutOfBounds
            try:
                minBases = int(sys.argv[index+1])
                index = index+1
            except ValueError:
                print('Error: \'minB\' must be followed by a number, to represent minimum number of bases for a position to be considered methylated.')
                sys.exit()

#------------- Methylation ratio requirement -------------#
        elif (arg == 'ratio' and index+1 <= len(sys.argv)-1): # avoid IndexOutOfBounds
            try:
                methylRatio = float(sys.argv[index+1])
                index = index+1
                if methylRatio > 1:
                    print('Ratio must be less than 1.0')
                    sys.exit()
            except ValueError:
                print('Error: \'ratio\' must be followed by a decimal number representing minimum methylation ratio for a position to be considered methylated.')
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
try:
    with open('./Batch/' + batch, 'r') as seqList:
        for el in seqList:
            try:
                el = el.strip('\n').split(',')[0]
                f = gzip.open(dirName+el, 'r')
                next(f) # verify this is safe for later, that it is not empty.
                confirmed.append(el)
                f.close()
            except FileNotFoundError:
                print('Data file ' + el + " was not found. Excluding.")
except NameError:
    print('Batch file was not given.')
    sys.exit()

if len(confirmed) < 2:
    print('Need minimum 2 valid methylomes to compare segregating sites.')
    sys.exit()


# -------- Prep the result file. As many names as there are options... --------
if defaultRange == 1:
    if minBases == '':
        if methylRatio == '':
            description = './Results/Segregating_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
        else:
            description = './Results/Segregating_methylation_ratio_'+str(methylRatio)+'_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
    elif methylRatio == '':
        description = './Results/Segregating_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
    else:
        description = './Results/Segregating_methylation_ratio_'+str(methylRatio)+'_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
else:
    if minBases == '':
        if methylRatio == '':
            description = './Results/Segregating_chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'
        else:
            description = './Results/Segregating_methylation_ratio_'+str(methylRatio)+'chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'
    elif methylRatio == '':
        description = './Results/Segregating_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'
    else:
        description = './Results/Segregating_methylation_ratio_'+str(methylRatio)+'_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'

try:
    result = open(description, 'w')
except FileNotFoundError: # When dir Results does not exist
    os.mkdir('./Results')
    result = open(description, 'w')
# HEADERS
result.write('total-seg-sites,total-perc-seg,CG-seg,CG-perc,CHG-seg,CHG-perc,CHH-seg,CHH-perc,ecotype1,ecotype2,chromosome,start-pos,end-pos,minBase,methylRatio\n')

if hrr == 1:
    if defaultRange == 1:
        report = open('./Results/SegReport_chr_'+str(chromo)+'_range_full_from_'+batch+'.txt', 'w')
    else:
        report = open('./Results/SegReport_chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch+'.txt', 'w')

# The branch for looking at whole methylomes
if chromo == 'all':

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

# keep same 0,1,2 code for 'check' but new determination for what is considered methylated
                            if methylRatio == '':
                                # branch where no ratio is given
                                if minBases == '':
                                    # no arguments given, so use given methylation_call value
                                    check = int(line2[6]) + int(line1[6])
                                elif int(line2[5]) >= minBases and int(line1[5]) >= minBases:
                                    # minBases given. Both must meet minBases to qualify
                                    check = int(line2[6]) + int(line1[6])
                                else: # considered unmethylated for not meeting minBases
                                    check = 0

                            # branch where methylRatio is given, but minBases is not
                            elif minBases == '':
                                # Only methylRatio given. Both must qualify to be considered
                                if int(line1[4]) / int(line1[5]) >= methylRatio:
                                    mc1 = 1
                                else:
                                    mc1 = 0
                                if int(line2[4]) / int(line2[5]) >= methylRatio:
                                    mc2 = 1
                                else:
                                    mc2 = 0
                                check = mc1 + mc2

                            # branch with methylRatio and minBases
                            elif int(line2[5]) >= minBases and int(line1[5]) >= minBases:
                                # since both have enough bases, we can evaluate based on methylRatio
                                if int(line1[4]) / int(line1[5]) >= methylRatio:
                                    mc1 = 1
                                else:
                                    mc1 = 0
                                if int(line2[4]) / int(line2[5]) >= methylRatio:
                                    mc2 = 1
                                else:
                                    mc2 = 0
                                check = mc1 + mc2
                            else: # considered unmethylated for not meeting minBases
                                check = 0

                            # matching
                            if check == 2: # line1[6] == 1 and line2[6] == 1:
                                if line1[3] == line2[3]: # same mc_class

                                    match = re.search(r'CG', line1[3])
                                    if match:
                                        CG_m += 1
                                        total_m += 1

                                    match = re.search(r'C[A,T,C]G', line1[3])
                                    if match:
                                        CHG_m += 1
                                        total_m += 1

                                    match = re.search(r'C[A,T,C][A,T,C]', line1[3])
                                    if match:
                                        CHH_m += 1
                                        total_m += 1

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
                        elif int(line1[1]) > int(line2[1]):
                            line2 = next(f2).decode('utf-8').strip('\n').split('\t')
                        elif int(line1[1]) < int(line2[1]):
                            line1 = next(f1).decode('utf-8').strip('\n').split('\t')

                    else: # dealing with unusual files, since the chr don't match.
                        while line1[0] > line2[0]:
                            line2 = next(f2).decode('utf-8').strip('\n').split('\t')
                        while line1[0] < line2[0]:
                            line1 = next(f1).decode('utf-8').strip('\n').split('\t')

            except StopIteration:
                print('End of comparison '+seq1+' & '+seq2)
                print(str(round(time.time() - subStart, 2)) + 'seconds for a comparison.')

            # Write results to file according to HEADERS: total-seg-sites,total-perc-seg,CG-seg,CG-perc,CHG-seg,CHG-perc,CHH-seg,CHH-perc,ecotype1,ecotype2,chromosome,start-pos,end-pos,minBase,methylRatio
            # the requested metric for each context divides a mc_class segregation over all positions with resolved methyl status
            # CHG_s/(CHG_s+CHG_m) is the format to give the percent of segregation exclusive to a particular mc_class
            result.write(str(total_s)+','+str(total_s/(total_s+total_m))+','+str(CG_s)+','+str(CG_s/(total_s+total_m))+',')
            result.write(str(CHG_s)+','+str(CHG_s/(total_s+total_m))+','+str(CHH_s)+','+str(CHH_s/(total_s+total_m))+',')
            result.write(seq1+','+seq2+','+str(chromo)+','+str(minPos)+','+str(maxPos)+','+str(minBases)+','+str(methylRatio)+'\n')

            # Optional textfile report
            if hrr == 1:
                if minBases == '':
                    mb = 'not applicable'
                else:
                    mb = str(minBases)
                if methylRatio == '':
                    mr = 'not applicable'
                else:
                    mr = str(methylRatio)
                report.write('Comparing '+str(seq1)+' to '+str(seq2)+ 'where chromosome of interest is '+str(chromo)+', range is '+str(minPos)+'-'+str(maxPos)+'.\n')
                report.write('Custom minimum required observed bases is '+mb+', custom ratio of observed methylation is '+mr+'.\n')
                report.write('Total Segregations: '+str(total_s)+', Total Percent Segregation: '+str(total_s/(total_s+total_m))+'\n')
                report.write('Mismatched contexts (1 or 2 methylations at a position where mc_class differs): '+str(contextMismatch)+'\n')
                report.write('CG Segregations: '+str(CG_s)+', Percent Segregation within CG-class: '+str(CG_s/(CG_s+CG_m))+', Percent CG Segregation overall: '+str(CG_s/(total_s+total_m))+'.\n')
                report.write('CHG Segregations: '+str(CHG_s)+', Percent Segregation within CHG-class: '+str(CHG_s/(CHG_s+CHG_m))+', Percent CHG Segregation overall: '+str(CHG_s/(total_s+total_m))+'.\n')
                report.write('CHH Segregations: '+str(CHH_s)+', Percent Segregation within CHH-class: '+str(CHH_s/(CHH_s+CHH_m))+', Percent CHH Segregation overall: '+str(CHH_s/(total_s+total_m))+'\n\n')

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

            total_m = total_s = CG_m = CG_s = CHH_m = CHH_s = CHG_m = CHG_s = contextMismatch = 0

            # -------- Now we are in the comparison part --------
            line1 = next(f1).decode('utf-8').strip('\n').split('\t') # get next lines
            line2 = next(f2).decode('utf-8').strip('\n').split('\t')

            # Set files to start at designated chr and minPos:
            while chromo > int(line1[0]):
                line1 = next(f1).decode('utf-8').strip('\n').split('\t') # get next lines
            while int(line1[1]) < minPos: # get to start position
                line1 = next(f1).decode('utf-8').strip('\n').split('\t')
            while chromo > int(line2[0]):
                line2 = next(f2).decode('utf-8').strip('\n').split('\t')
            while int(line2[1]) < minPos:
                line2 = next(f2).decode('utf-8').strip('\n').split('\t')

            try:
                while chromo == int(line1[0]) and chromo == int(line2[0]):
                    if  line1[1] == line2[1]: # positions match. Compare:

                    # leave when we go past taget range
                        if int(line1[1]) > maxPos:
                            break

# keep same 0,1,2 code for 'check' but new determination for what is considered methylated
                        if methylRatio == '':
                            # branch where no ratio is given
                            if minBases == '':
                                # no arguments given, so use given methylation_call value
                                check = int(line2[6]) + int(line1[6])
                            elif int(line2[5]) >= minBases and int(line1[5]) >= minBases:
                                # minBases given. Both must meet minBases to qualify
                                check = int(line2[6]) + int(line1[6])
                            else: # considered unmethylated for not meeting minBases
                                check = 0

                        # branch where methylRatio is given
                        elif minBases == '':
                            # Only methylRatio given. Both must qualify to be considered
                            if int(line1[4]) / int(line1[5]) >= methylRatio:
                                mc1 = 1
                            else:
                                mc1 = 0
                            if int(line2[4]) / int(line2[5]) >= methylRatio:
                                mc2 = 1
                            else:
                                mc2 = 0
                            check = mc1 + mc2
                        elif int(line2[5]) >= minBases and int(line1[5]) >= minBases:
                            # since both have enough bases, we can evaluate based on methylRatio
                            if int(line1[4]) / int(line1[5]) >= methylRatio:
                                mc1 = 1
                            else:
                                mc1 = 0
                            if int(line2[4]) / int(line2[5]) >= methylRatio:
                                mc2 = 1
                            else:
                                mc2 = 0
                            check = mc1 + mc2
                        else: # considered unmethylated for not meeting minBases
                            check = 0

                        if check == 2: # matching
                            if line1[3] == line2[3]:

                                match = re.search(r'CG', line1[3])
                                if match:
                                    CG_m = CG_m + 1
                                    total_m += 1

                                match = re.search(r'C[A,T,C]G', line1[3])
                                if match:
                                    CHG_m += 1
                                    total_m += 1

                                match = re.search(r'C[A,T,C][A,T,C]', line1[3])
                                if match:
                                    CHH_m += 1
                                    total_m += 1

                            else:
                                contextMismatch = contextMismatch +1

                        elif check == 1: # segregation
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
                    elif int(line1[1]) > int(line2[1]):
                        line2 = next(f2).decode('utf-8').strip('\n').split('\t')
                    elif int(line1[1]) < int(line2[1]):
                        line1 = next(f1).decode('utf-8').strip('\n').split('\t')

                # the chromosome of interest is complete.
                print('End of comparison '+seq1+' & '+seq2)
                print(str(round(time.time() - subStart, 2)) + 'seconds for a comparison.')

            except StopIteration: # if we stop early which we shouldnt
                print('~End of comparison '+seq1+' & '+seq2)
                print(str(round(time.time() - subStart, 2)) + 'seconds for a comparison.')

                # Write results to file according to HEADERS: total-seg-sites,total-perc-seg,CG-seg,CG-perc,CHG-seg,CHG-perc,CHH-seg,CHH-perc,ecotype1,ecotype2,chromosome,start-pos,end-pos,minBase,methylRatio
                result.write(str(total_s)+','+str(total_s/(total_s+total_m))+','+str(CG_s)+','+str(CG_s/(total_s+total_m))+',')
                result.write(str(CHG_s)+','+str(CHG_s/(total_s+total_m))+','+str(CHH_s)+','+str(CHH_s/(total_s+total_m))+',')
                result.write(seq1+','+seq2+','+str(chromo)+','+str(minPos)+','+str(maxPos)+','+str(minBases)+','+str(methylRatio)+'\n')

                # Optional textfile report
                if hrr == 1:
                    if minBases == '':
                        mb = 'not applicable'
                    else:
                        mb = str(minBases)
                    if methylRatio == '':
                        mr = 'not applicable'
                    else:
                        mr = str(methylRatio)
                    report.write('Comparing '+str(seq1)+' to '+str(seq2)+ 'where chromosome of interest is '+str(chromo)+', range is '+str(minPos)+'-'+str(maxPos)+'.\n')
                    report.write('Custom minimum required observed bases is '+mb+', custom ratio of observed methylation is '+mr+'.\n')
                    report.write('Total Segregations: '+str(total_s)+', Total Percent Segregation: '+str(total_s/(total_s+total_m))+'\n')
                    report.write('Mismatched contexts (1 or 2 methylations at a position where mc_class differs): '+str(contextMismatch)+'\n')
                    report.write('CG Segregations: '+str(CG_s)+', Percent Segregation within CG-class: '+str(CG_s/(CG_s+CG_m))+', Percent CG Segregation overall: '+str(CG_s/(total_s+total_m))+'.\n')
                    report.write('CHG Segregations: '+str(CHG_s)+', Percent Segregation within CHG-class: '+str(CHG_s/(CHG_s+CHG_m))+', Percent CHG Segregation overall: '+str(CHG_s/(total_s+total_m))+'.\n')
                    report.write('CHH Segregations: '+str(CHH_s)+', Percent Segregation within CHH-class: '+str(CHH_s/(CHH_s+CHH_m))+', Percent CHH Segregation overall: '+str(CHH_s/(total_s+total_m))+'\n\n')

            f1.close()
            f2.close()
            index2 = index2 +1 # end of inner

        index1 = index1 +1 # end of outer

result.close()
report.close()
print(str(round(time.time() - start,2)) + ' seconds to completion.')

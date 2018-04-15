import os
import sys
import gzip
import re
import time

def eval(line, minBases, methylRatio):
    if methylRatio == '':
        # branch where no ratio is given
        if minBases == '':
            # no arguments given, so use given methylation_call value
            return int(line[6])
        elif int(line[5]) >= minBases:
            # minBases requirement met.
            return int(line[6])
        else: # considered unknown for not meeting minBases
            return -1
    # branch where methylRatio is given
    elif minBases == '':
        if int(line[4]) / int(line[5]) >= methylRatio:
            return 1
        else:
            return 0
    elif int(line[5]) >= minBases:
        if int(line[4]) / int(line[5]) >= methylRatio:
            return 1
        else:
            return 0
    else: # considered unknown for not meeting minBases req.
        return -1

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
minPos = 0
maxPos = 31000000 # 31 million covers the longerst chromosome.
chromo = 'all' # will consider all chromosomes by default
defaultRange = 1
minBases = ''
methylRatio = ''
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

# -------- Prep the result file. As many names as there are options... --------
if defaultRange == 1:
    if minBases == '':
        if methylRatio == '':
            description = './Results/Summary_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
        else:
            description = './Results/Summary_methylation_ratio_'+str(methylRatio)+'_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
    elif methylRatio == '':
        description = './Results/Summary_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
    else:
        description = './Results/Summary_methylation_ratio_'+str(methylRatio)+'_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_full_from_'+batch.split('.')[0]+'.csv'
else:
    if minBases == '':
        if methylRatio == '':
            description = './Results/Summary_chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'
        else:
            description = './Results/Summary_methylation_ratio_'+str(methylRatio)+'chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'
    elif methylRatio == '':
        description = './Results/Summary_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'
    else:
        description = './Results/Summary_methylation_ratio_'+str(methylRatio)+'_min_bases_'+str(minBases)+'_chr_'+str(chromo)+'_range_'+str(minPos)+'-'+str(maxPos)+'_from_'+batch.split('.')[0]+'.csv'

try:
    result = open(description, 'w')
except FileNotFoundError: # When dir Results does not exist
    os.mkdir('./Results')
    result = open(description, 'w')
result.write('overall-methylation,CG_context,CHH_context,CHG_context,undetermined-status,ecotype,chromosome,start-pos,end-pos,minBase,methylRatio\n')

# start reading sequences from designated batch-file
with open('./Batch/' + batch, 'r') as seqList:
    for line in seqList:
        fname = line.strip('\n').split(',')[0]
        subStart = time.time()
        # Reset for every new file
        CG_methyl=CG_total=CHH_methyl=CHH_total=CHG_methyl=CHG_total=linecount=all_methyls=unk = 0

        # Seperate branch of opertions for default behaviour. Doesnt feel great,
        # quasi-duplicating my code this way, but it helps keep both functionalities
        # within this one summarize.py script

        try:
            with gzip.open(dirName+fname, 'r') as f:
                next(f)# skip title headers

                for byteLine in f:
                    line = byteLine.decode('utf-8').strip('\n').split('\t')
                    # Check all the reasons NOT to look at this base
                    if chromo != 'all':
                        if chromo == int(line[0]):
                            if int(line[1]) < minPos:
                                line = next(f).decode('utf-8').strip('\n').split('\t')
                                continue
                            elif int(line[1]) > maxPos:
                                break
                        elif chromo > int(line[0]):
                            line = next(f).decode('utf-8').strip('\n').split('\t')
                            continue

                        else: # past target chromosome
                            break

                    # in target region. Evaluate
                    res = eval(line, minBases, methylRatio)

                    if res == -1:
                        unk += 1 # keep track of undetermined positioins
                    linecount += 1 # to calc real % methyls, % unknown, % non-methyl

                    match = re.search(r'CG', line[3])
                    if match:
                        #CG_total += 1
                        if res == 1:
                            CG_methyl += 1
                            all_methyls += res

                    match = re.search(r'C[A,T,C]G', line[3])
                    if match:
                        #CHG_total += 1
                        if res == 1:
                            CHG_methyl += 1
                            all_methyls += res

                    match = re.search(r'C[A,T,C][A,T,C]', line[3])
                    if match:
                        #CHH_total += 1
                        if res == 1:
                            CHH_methyl += 1
                            all_methyls += res

                    # End of line-loop. Gets next line in myzip file

            # on finishing all lines of a file,
            m_overall = str(round(all_methyls / linecount, 6)) #
            unk_count = str(unk)
            # CG_context = str(round((CG_methyl / CG_total), 6))
            # CHH_context = str(round((CHH_methyl / CHH_total), 6))
            # CHG_context = str(round((CHG_methyl / CHG_total), 6))
            CG_context = str(round((CG_methyl / linecount), 6))
            CHH_context = str(round((CHH_methyl / linecount), 6))
            CHG_context = str(round((CHG_methyl / linecount), 6))

            # write results to csv file to save answers
            result.write(m_overall+','+CG_context+','+CHH_context+','+CHG_context+','+unk_count+','+fname+','+str(chromo)+','+str(minPos)+','+str(maxPos)+','+str(minBases)+','+str(methylRatio)+'\n')
            # timer per file
            print(str(round(time.time() - subStart, 2)) + ' seconds for '+ fname)

        except FileNotFoundError:
            print('File ' + fname + " was not found.")

result.close()
print(str(round(time.time() - start,2)) + ' seconds to completion.')

# The objective here is to better organise our methylome sequences.
# With the previously constructed masterList, we move all sequence files not on
# this masterList to a parrallel directory, Unmatched. This way we aren't deleting
# any data and keep it out of the way for future potential use.

# STEP 1
with open('masterList.csv', 'r') as src:
    next(src) # because headers
    confirmed = []
    for line in src:
        try:
            # python list indexing lets us grab the last element this way, which is filename
            f = line.split(',')[-1]
            confirmed.append(f.strip('\n'))
        except IndexError:
            print ('-------------IndexError---------'+line)

# STEP 2
with open('move.bash', 'w') as dest:
    # longlist refers to the output of `ls -1` on Data/
    # The longlist.txt given has all the filenames that I got from downloading the Data, including duplicates.
    # The original data had duplications indicated with (1) in name. They can be removed.
    with open('BuildingBlocks/sequenceFiles.txt', 'r') as src:
        dest.write('#!/bin/bash\n')
        for filename in src:
            if filename.strip('\n') not in confirmed:
                dest.write('mv \"Data/' + filename.strip('\n') + '\" Unmatched\n')
# The script is made, now just run it from command line: bash move.bash

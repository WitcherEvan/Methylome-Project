import re

with open('file-per-line.txt', 'r') as files:
    # Create dictionary for GSM# -> filename
    numToFile = {}
    for filename in files:
        numToFile[filename.split('_')[0]] = filename.strip('\n')

with open('query.csv', 'r') as info:
    next(info) # skip headers

    # Create dictionaries for eid -> row and name -> row
    eidRow = {}
    nameRow = {}
    for row in info:
        eidRow[row.split(',')[0]] = row.strip('\n')
        nameRow[row.split(',')[1]] = row.strip('\n')

final = open('associated.csv', 'w') # where we write completed data
final.write('tg_ecotypeid,name,CS_number,country,latitude,longitude,collector,seq_by, filename') # headers

with open('gsm-title.csv', 'r') as numTitle:
    remainder = []
    for row in numTitle:
        # looking for tg_ecotypeId in filename. Most have this: 704/926
        match = re.search(r'\([0-9]*\)', row)

        if match: # search() returned something: ecotypeId was available
            # find data associated with this eid. May not exist.
            eid = match.group(0).strip('()')
            dataRow = eidRow.get(eid)

            if dataRow != None:
                # Since eID relates to an entry we have gps data on,
                # we can do the assciations, add new csv row
                gsm = row.split(',')[0]
                filename = numToFile.get(gsm)
                if filename != None:
                    final.write(dataRow + ',' + filename + '\n')
                    continue
                # go to next entry, doesnt put this in 'remainder' pile

        # failing above conditions brings us here.
        # We store this entry to handle later
        remainder.append(row.strip('\n'))

# if 'name' from dict is in a remainder entry, then we have a match by Name
names = nameRow.keys()
exluded = []

for entry in remainder:
    messyName = entry.split(',')[1]
    cleanName = messyName.split('_(')[0]

    if cleanName in names:
        dataRow = nameRow.get(cleanName)
        gsm = entry.split(',')[0]
        filename = numToFile.get(gsm)

        if (filename != None):
            final.write(dataRow + ',' + filename + '\n')
            continue # when successful

    # Store files that did not pair to any gps data.
    # Good for verifying that we haven't missed any.
    exluded.append(entry + '\n')

final.close()
with open('excluded.txt', 'w') as ex:
    ex.write(''.join(exluded))

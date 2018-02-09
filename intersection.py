with open('compact.csv', 'r') as gpsfiles:
    # We make the unique names keys to their entire csv row.
    rowInfo= {}
    for row in gpsfiles:
        cols = row.split(',')
        key = cols[0]
        rowInfo[key] = row.strip('\n')

with open('file-per-line.txt', 'r') as allNames:
    with open('confirmed-sequences3', 'w') as dest:
        # We match every fileName to the appropriate row.
        # We write the new csv file to contain this connection
        # between gps coordinates and related sequence file.
        for longname in allNames:
            # in progress
            #ifname)

            # for key in rowInfo.keys():
            #     if (key in el):
            #         dest.write(rowInfo.get(key) + ',' + el + '\n')
            #         break

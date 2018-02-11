with open('all927.csv','r') as src:
    with open('condensed.csv', 'w') as dest:
        #next(src) # skip header descriptions
        for row in src:
            cols = row.split(',')
            dest.write(cols[0]  + ',' + cols[1] + '\n')

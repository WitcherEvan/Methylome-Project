with open('query.csv','r') as src:
    with open('hamsterCustomCoord.txt', 'w') as dest:
        next(src)
        for row in src:
            cols = row.split(',')
            dest.write(cols[4]  + '\t' + cols[5] + '\t' + 'circle1' + '\t' + 'green' '\n')

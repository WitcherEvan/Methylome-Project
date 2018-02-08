src = open('query.csv','r')
dest = open('hamsterCustomCoord.txt', 'w')
for row in src:
    cols = row.split(',')
    dest.write(cols[4]  + '\t' + cols[5] + '\t' + 'circle1' + '\t' + 'green' '\n')

src.close()
dest.close()

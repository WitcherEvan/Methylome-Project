with open('query.csv','r') as src:
    names = []
    count = 0
    for row in src:
        cols = row.split(',')
        names.append(cols[1])
        count = count + 1
while (count > 1):
    aName = names.pop()
    count = count - 1
    if (aName in names):
        print(aName)

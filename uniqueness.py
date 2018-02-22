# A debugging tool I modify to look for duplications in csv's

with open('masterList.csv','r') as src:
    items = []
    count = 0
    for row in src:
        cols = row.split(',')
        items.append(cols[len(cols)-1])
        count = count + 1
while count > 1:
    aItem = items.pop()
    count = count - 1
    if aItem in items:
        print(aItem)
print('Unique if this is the only line!')

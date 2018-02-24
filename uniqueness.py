# A debugging tool I modify for many things

with open('masterList.csv','r') as src:
    items = []
    count = 0
    with open('checkL2', 'w') as dest:
        for row in src:
            dest.write(row.split(',')[-1])
            #items.append(row.split(',')[-1])
            count = count + 1
for
# # uniqueness checking
# while count > 1:
#     aItem = items.pop()
#     count = count - 1
#     if aItem in items:
#         print(aItem)
# print('Unique if this is the only line!')

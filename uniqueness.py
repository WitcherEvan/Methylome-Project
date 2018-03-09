# A debugging tool I modify for many things

# with open('masterList.csv','r') as src:
#     items = []
#     count = 0
#     with open('checkL2', 'w') as dest:
#         for row in src:
#             dest.write(row.split(',')[-1])
#             #items.append(row.split(',')[-1])
#             count = count + 1

# # To prove that base-pair position depends on chromosome.
# # Each chromosome goes from 1 to ~30M. Not one continuous thing
# import gzip
# with gzip.open('./Data/GSM2099121_allc_10015.tsv.gz', 'r') as f:
#     next(f)
#     try:
#         count = 0
#         last = 0
#         for l in f:
#             count = count + 1
#
#             col = l.decode('utf-8').split('\t')
#             cur = int(col[1])
#             if cur < last:
#                 print(col)
#                 print('Checking for bpp rest on new chromosome')
#                 break
#             last = cur
#     except ValueError:
#         print(l)
#         print (count)

# # uniqueness checking
# while count > 1:
#     aItem = items.pop()
#     count = count - 1
#     if aItem in items:
#         print(aItem)
# print('Unique if this is the only line!')

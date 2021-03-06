print('\n**************************************************************************************\n')
print('Quick reference for python scripts. For more deatils, see README.md')
print('\n**************************************************************************************\n')
print('summarize.py arguments:')
print('--- chr [1-5] -------- optional, paired. Specifies a chromosome to summarize.')
print('--- range #-# -------- optional, paired. Specifies a range on chromosome to summarize.')
print('--- minB # ----------- optional, minimum observations for valid methylation status.')
print('--- ratio # ---------- optional, minimum ratio for methylated status.')
print('--- Batch_file ------- required. Name of file in Batch.')
print('Output: Result csv file.')
print('\n**************************************************************************************\n')
print('segSites.py arguments:')
print('--- chr [1-5] -------- optional, paired. Specifies a chromosome to summarize.')
print('--- range #-# -------- optional, paired. Specifies a range on chromosome to summarize.')
print('--- minB # ----------- optional, minimum observations for valid methylation status.')
print('--- ratio # ---------- optional, minimum ratio for methylated status.')
print('--- Batch_file ------- required. Name of file in Batch.')
print('--- HRR -------------- optional. Produces a Human Readable Report in Results.')
print('Output: Result csv file, optional Result txt file.')
print('\n**************************************************************************************\n')
print('subdivide.py arguments:')
print('--- Batch_file ------- required. Name of file in Batch.')
print('--- chr [1-5] #-# ---- minimum one required. Indentifies the range on a chromosome')
print('---------------------- which is saved to resulting Data file.')
print('---------------------- Adds one chr at a time. Repeat for each desired chromosome.')
print('Output: Data file, Batch file.')
print('\n**************************************************************************************\n')
print('createPop.py arguments:')
print('--- kmRadius --------- optional. Cutoff range from central ecotype.')
print('--- popSize ---------- optional. Max number of ecotypes in simulated popualtion.')
print('--- ecotypeID # ------ Pick central ecotype by tg_ecotypeid.')
print('--- ecotypeName name - Pick central ecotype by name.')
print('--- GSMnumber # ------ Pick central ecotype by GSM part of filename, include \'GSM\'.')
print('Output: Batch file with distances.')
print('\n**************************************************************************************\n')

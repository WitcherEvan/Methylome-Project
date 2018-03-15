# Methylome-Project
Developping tools for the analysis of large methylome files.

intersection.py does inital work of matching gps coordinates to their corresponding methylome, where such matches exist. It requires use of GPSdata.csv, gsm-title.csv, and sequenceFiles.txt. It produces masterList.csv, which contains data matched to its methylome sequence. As a secondary output, it creates excluded.txt to hold the names of unmatched methylome files.
In the event that we have new data, methylome sequences and/or gps coordinates, update sequenceFiles.txt (output of ls -1 in Data directory) and/or GPSdata.csv (by adding rows), then run intersection.py.
Currently, 792 of 922 methylome sequences were paired to their gps-data.

Methylome sequence files were not named uniformly. Pairing data with sequences required using both 'name' and 'tg_ecotypeid', from the gsm-title file (from ncbi).
The included csv files have the dashes and spaces changed to underscores, to match the format of filenames (on Ubuntu/Linux, that is). You can see in all927.csv, directly from ncbi site, that the format was different, and needed standardizing.

For housekeeping purposes, there exists culling.py. After running intersection.py and getting masterList.csv, running culling.py will create the Bash script move.bash. Running this script sends all excluded methylomes from Data to Unmatched.

read.py has command-line language to make specific queries about sequences. Computes fraction of methylation overall and for each cytosine context. It also outputs to command-line the time taken on completing each methylome.
--- Command line arguments ---
None are required, and no ordering enforced.
Keep 'chr' and 'range' paired with their value as indicated.
-- chr [1-5]    tag is chr, followed by an integer 1-5 for chromosome number
-- range #-#    range over which we will analyze. Numbers must be positive integers seperated by '-'.
-- filename     a file containing the filenames of sequences we are interested in. These files should be in Batch directory.

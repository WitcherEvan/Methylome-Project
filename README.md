# Methylome-Project
791 of 922 sequence files were paired to their gps-data.
Adds the filenames of sequences we have to their gps coordinates and other data, should they exist.
There are more gps-data entries than sequences, so not every gps-data entry is useful.
Equally, downloaded sequences may not have any associated gps-data.
Sequence files were not named uniformly. Pairing data and sequences required using both 'name' and 'tg_ecotypeid', from the gsm-title file (from ncbi).
Use the included csv files, because they have the dashes and spaces changed to underscores, to match the format of filenames (on Ubuntu/Linux, that is). You can see in all927.csv, from ncbi site, that the format was different, and needed standardizing.

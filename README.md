# Methylome-Project
Developping tools for the analysis of large methylome files.

**intersection.py**
intersection.py does inital work of matching gps coordinates to their corresponding methylome, where such matches exist. It requires use of GPSdata.csv, gsm-title.csv, and sequenceFiles.txt. It produces masterList.csv, which contains data matched to its methylome sequence. As a secondary output, it creates excluded.txt to hold the names of unmatched methylome files.
In the event that we have new data, methylome sequences and/or gps coordinates, update sequenceFiles.txt (output of ls -1 in Data directory) and/or GPSdata.csv (by adding rows), then run intersection.py.
Notes:
Currently, 792 of 922 methylome sequences were paired to their gps-data.
Three methylomes paired to respective entries in the given GPSdata.csv are missing GPS coordinates. By filename, they are GSM2099277_allc_7236.tsv.gz and GSM1085274_mC_calls_Litva.tsv.gz which share tg_ecotypeid, and GSM2099312_allc_7427.tsv.gz.

***About methylome Data filenames***
Methylome sequence files were not named uniformly. Pairing data with sequences required using both 'name' and 'tg_ecotypeid', from the gsm-title file (from ncbi).
The included csv files have the dashes and spaces changed to underscores, to match the format of methylome filenames (on Ubuntu/Linux, that is). You can see in all927.csv, directly from ncbi site, that the format was different, and needed standardizing. We conform to the formatting of the methylome sequence files so that more files can be included without modifying them.

**culling.py and move.bash**
For housekeeping purposes, there exists culling.py. After running intersection.py and getting masterList.csv, running culling.py will create the Bash script move.bash. Running this script sends all excluded methylomes from Data to Unmatched.

**summarize.py**
summarize.py computes fraction of methylation overall and for each cytosine context. Uses command-line language to summarize specific regions of methlyomes. Outputs to command-line the time to complete each methylome.
--- Command line arguments ---
None are required. No ordering enforced. Keep 'chr' and 'range' paired with their value as indicated.
-- chr [1-5]    tag is chr, followed by an integer 1-5 for chromosome number
-- range #-#    range over which we will analyze. Numbers must be positive integers seperated by '-'.
-- minB         Specify a minimum number of observations at a position necessary for methylation status to be included.
-- ratio        Specify a minimum ratio of methylations over observations that must be met for a position to be considered methylated.
-- filename     the file containing multiple filenames of sequences we are interested in. This 'filename' should be in Batch directory.

**segSites.py**
segSites.py compares two methylomes to count segregating sites (where the same positions have same mc_class but have different methylation status). Sites with same methylation status but different mc_class are not currently considered a valid segregating site, but are counted seperately in case it can give some helpful information.
Functions just like summarize.py for command-line arguments, with one additional option. As long as the 'filename' found in Batch has 2 or more methylomes, it compares all of them to each other, to produce number and percent of Segregating Sites, across all contexts. 1 - % segregation = % shared methylation.
--- Command line arguments ---
None are required, no ordering enforced. Keep 'chr' and 'range' paired with their value as indicated.
-- chr [1-5]    chr is followed by an integer 1-5 for chromosome number of interest. Default consideres all chromosomes.
-- range #-#    range over which we will analyze. Numbers must be positive integers seperated by '-'. Default range covers everything.
-- minB         Specify a minimum number of observations at a position necessary for methylation status to be included.
-- ratio        Specify a minimum ratio of methylations over observations that must be met for a position to be considered methylated.
-- filename     the file containing the filenames of sequences we are interested in. These files should be in Batch directory.
-- HRR          Short for Human Readable Report, this argument will write a textfile describing the result of each comparison.

**createPop.py**
createPop.py creates a population of ecotypes based on closest proximity to a central ecotype we designate. Requires masterList.csv, currently in the 'root' directory, same as createPop's script.
We specify a radius around the central ecotype, then give a maximum number of ecotypes to include in the population.
As output, a new Batch file is made, where each line has a second word, seperated by a comma ',' giving the distance (in kilometers) from a line's methylome to the central ecotype.
--- Command line arguments ---
-- kmRadius #    Max distance from central. distance < kmRadius is a necessary condition to be included in result.
-- popSize #     Max number of best (closest) ecotypes to be included in result.
--> One of ecotypeID, ecotypeName, or GSMnumber is required to identify central ecotype. If multiple are given in command line, the last will take precedence.  Follow with appropriate value, ex:
    ecotypeID 8424
    ecotypeName Kas_2
    GSMnumber GSM2099347

***About Batch files***
There is currently no required naming or extension (.txt .csv) for Batch files. It is the contents of Batch files that have restrictions:
One filename per line as the first word on said line. That filename is the complete name of a methylome sequence file held in Data.
Batch files made via createPop will have a comma-seperated distance for a second word per line. This is to keep information about the simulated pop.
So long as the filename is seperated by a comma, anything else can follow on that line. If that is a restriction we do not want, tweak the one line in segSites.py and summarize.py where a line of the given Batch file is parsed to open the specified sequence file from Data directory.

***About naming conventions***
The naming format used on files made in these scripts can be lengthly at times. This is so that filenames/titles are
--> unique, will only be overwritten by the same funciton and same arguments
--> directly identifing the parameters used to produce the result.
This depends on Batch files also being unique. Example:
Running a command with certain arguments including Batch file 'bf1' creates some Result file. Say 'bf1' is changed. Running the same command with same arguments produces a modified Result which overwrites the previous, because they would share the same name.
Pro: Avoids creating copies of the same Result, reducing clutter.
Con: Could loose Results if human-written Batch files change contents but keep same name.
Fix: Add a timestamp to Result file naming format. Makes names longer and doesn't provide new information, but means no Result will be overwritten.

**subdivide.py**
subdivide.py creates new methylome sequence files in Data by copying selected regions of given original files. Can be used with a subdivided Data file.
--- Command line arguments ---
-- filename    Takes the filename of a Batch file, containing the names of Data files to be copied and modified.
Slightly more complex format to designate desired regions:
-- chr [1-5] minPosition-maxPosition    chr is the key word, indicating that the next two argments are a chromosome number, then the range of interest on said chromosome. Repeat for every desired chromosome.
We use Batch files not only to keep with the style of other scripts, but also because it makes sense that we want to compare multiple methylomes over specific regions.
The 'chr' could be removed for brievity, but is included for clarity, punctuating what would otherwise be a mess of numbers.

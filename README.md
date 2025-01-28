PARERSv2 WIP as of 11-25-24

What you need:
- Python
- R
- BBTools suite
- Java (dependecy for BBMerge)
- MUSCLE alignment software (script runs from command line)
- AmpliconsRaw input file (provided)
- InputTemp input file (template provided)
- R1 and R2 sequencing files (folder with example files has been provided)
  - The files I have provided are truncated versions of the actual read files we got off the sequencer. Each contains 200 of the original ~1,000,000 sequences. I uploaded truncated versions of the files to speed up processing (for two files took 18 seconds on my work computer) and because even the zipped files far exceed the 25 mb limit GitHub will let me upload. I have uploaded the output folder I generated when I ran these two sets of data so you can compare your outputs. 

User guide in brief:
1. Download the dependencies listed above plus python and R packages listed in "environment_versions_TR_BD_JC_11-22-24"
    - This file plus additional background information on the pipeline can be found in the "PARERSv2_info_for_RSC" folder
    - Note: not all the packages in the spreadsheet (e.g. Rpy2) are used in this pipeline, but I included them so you could 100% recreate the run environment
2. Fill out your input temp file. To replicate the results I have included in the repo you don't need to change everything, but below is what you do/can.
  - Optional: output file name
  - Mandatory: output directory, path to BBmerge, path to MUSCLE, path to R, path to folder containing R scripts (just the folder, the script calls out the script names specifically)
3.  Define the paths to your R1 and R2 read files in parallel and assemble the files into separate R1 and R2 lists
4.  Make a list of your cell line names and define the control
5.  Define the path to your input file
6.  Define the path to your amplicons raw fasta file
7.  Since you are running the script in a linux environment, you may need to make adjustments to some of the portions of the script that call for the command line (ie BBMerge, MUSCLE, R scripts). The line locations of these are:
  - BBMerge = 271
  - MUSCLE = 629
  - R scripts = 2480, 2499, 2517, 2540, 2557

Support files (all in the info folder):
- Powerpoint with background information on RNA editing, the lab work done prior to data analysis, and the PARERS pipeline
- Notes from our meeting on 11-22
- Excel doc containing the programs and packages present on Tyler, Jason, and Brittney's computers
- Text document containing comments for what every line of code in the script is doing. This is extremely detailed and more about code function than intent, so  I don't expect it will be of useful for you, but if you are confused about what a particular chunk is doing this may provide some insight. You are also welcome to reach out to Tyler for clarification. 

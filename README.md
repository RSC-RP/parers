# Pipeline Analyzing RNA Editing RNA Sequencing (PARERS)
Authors: Rodshagen, Tyler; Tracy, Maxwell; Davidge, Brittney; Carnes, Jason; Morton, Glenn  
Maintainer: Rodshagen, Tyler  
[Stuart Lab](https://www.seattlechildrens.org/research/centers-programs/global-infectious-disease-research/research-areas-and-labs/stuart-lab/), [Center for Global Infectious Disease](https://www.seattlechildrens.org/research/centers-programs/global-infectious-disease-research/research-areas-and-labs/)  
[Research Scientific Computing](https://github.com/RSC-RP)  
[Seattle Childrens Research Institute](https://www.seattlechildrens.org/research/research-institute/  )

## Environmental Dependencies
- Python (3.12+ with pandas, python-docx, biopython, xlsxwriter)
- R (4.4 with tidyverse)
- [BBMap](https://sourceforge.net/projects/bbmap/) for BBMerge
- [MUSCLE](https://github.com/rcedgar/muscle) alignment software

You also need:
- [`AmpliconsRaw`](./test_data/AmpliconsRaw.fasta) input file (fasta)
- R1 and R2 sequencing files
  - There is test_data provided as examples
  - [`curated`](./test_data/curated/) contains synthetic examples around the A6 gene to demonstrate each type of editing event.
  - [`truncated`](./test_data/truncated/) contains 200 of an original ~1,000,000 sequences from MURF2 gene.
- `parers.cfg` - a configuration file that informs the pipeline on how to run the sample(s)

## Configure `parers.cfg`
It is recommended to use either of these as a template and modify for your samples:
- [`./test_data/curated/parers.cfg`](./test_data/curated/parers.cfg)
- [`./test_data/truncated/parers.cfg`](./test_data/truncated/parers.cfg)

**Note**: You may add as many pairs of R1, R2 sequences as you want as long as they are configured the same. i.e. Don't mix genes or primers.

## Apptainer
### Build
Provided is [`parers.def`](./parers.def) which is an [Apptainer](https://apptainer.org/docs/user/main/index.html) definition file.
``` bash
sudo apptainer build parers.sif parers.def
```
### Open interactive shell inside the container
Once built, launch the container and ensure the directories mentioned in `parers.cfg` are in the [bind paths](https://apptainer.org/docs/user/main/bind_paths_and_mounts.html) for the container.
``` bash
apptainer shell --bind /path/to/include parers.sif
```
### Run the pipeline on your configured sample in the container
The PARERS pipeline can be run by going to the directory with your `parers.cfg` and then running `python3 /parers/parers`.
``` bash
cd /path/to/parers.cfg
python3 /parers/parers
```

## Conda/Mamba
### Setup
Also provided is a [`parers.yml`](./env/parers.yml) so you can [create the parers conda or mamba environment from yml file](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

In order to set up the script to run in your environment, you need to tell it where somethings are:
- `bbmerge.sh` script
- `muscle` binary
- `Rscript` binary
- `R_for_cmd` directory

Activate your parers conda environment that was created from the [`parers.yml`](./env/parers.yml) in the step.
Now you will be able to use `whereis` to find the full path to `bbmerge.sh` script and each of the `muscle` and `Rscript` binaries.
``` bash
whereis bbmerge.sh
whereis muscle
whereis Rscript
```
You also need to note the full path of the `R_for_cmd` directory.

You will need to edit the following variables in [`parers.py`](./parers.py) with the appropriate paths:
``` python
path_to_bbmerge = "/path/to/bbmerge.sh"
path_to_muscle = "/path/to/muscle"
path_to_r = "/path/to/Rscript"
path_to_r_scripts = "/path/to/R_for_cmd"
```
### Run the pipeline on your configured sample in the conda/mamba environment
The PARERS pipeline can be run by going to the directory with your `parers.cfg` and then running `python3 /parers/parers`.
``` bash
cd /path/to/parers.cfg
python /path/to/parers.py
```

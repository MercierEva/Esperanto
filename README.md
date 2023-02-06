# Esperanto

## Description 

An app to analyse Sequencing Data from ONT technology.

This app allows to browse or drag and drop some fastq files. The user may choose few parameters to filter data. Then, all that remains is to launch the program. 
The workflow in the background is building with snakemake.
The result is a sequence Fasta consensus with or without the variant calling and a statistical report. 

## Install with Conda

### Requirements

#### Linux machine Ubuntu 20.04 or Ubuntu 22.04
- conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html>
- libgtk-3-dev 
- libreoffice

#### Install
`sudo apt-get install libgtk-3-dev libreoffice`
 
`conda create --name esperanto --file conda_env_for_esperanto.yml`

`conda activate esperanto`

#### Running App
`cd Esperanto`

`python3.8 u_app.py` 

or 

`python3.10 u_app.py`


## Or Install with Docker (But each docker build leads to rebuild each conda env)

### Requirements 

- Docker <https://docs.docker.com/engine/install/>

### Install and Running App
`cd Esperanto && mkdir data`

data directory will be the directory where you may put the fastq.gz files in order to access them from the Docker container.

`docker build -t esperanto . `

`docker run -it --env="DISPLAY" --net=host -v ${PWD}/workflow:/app/workflow esperanto:latest //from a VM`


## User Guide

- Upload one or more of fastq.gz files
- Choose a name for the final repository of the results (default 'WorkingSpace_YYYYMMDD')
- Set min and max length for the step of data filtration
- Set minimum of quantity of reads after the filtration step (Q17-->Q7 progressive to exceed this threshold parameter)
- Set type of sequencing reads
- Set number of threads (minimum 4 is recommended)
- Set Reverse Primer for the orientation of sequences 

## Go to Results

The results are located inside the **workflow/** directory in the directory with the name you specified. (for example **workflow/WorkingSpace_YYYYMMDD/**) 

Main files are :

- All_Consensus_fastas.fasta : Display of consensus sequences after the variant calling, so with degenerated nucleotides.
- All_fastas.fasta : Display of sequences with major nucleotides.
- StatisticReport.tsv : All data from the analysis with fasta sequence.





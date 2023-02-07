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
- libwebkitgtk-1.0-0
- libffi6
- libcanberra-gtk-module
- dbus-x11

#### Install
`sudo add-apt-repository universe`

`sudo apt update`

`sudo apt install -y libwebkitgtk-1.0-0 libffi6 libgtk-3-dev libreoffice libcanberra-gtk-module dbus-x11`

`conda create --name esperanto --file conda_env_for_esperanto.yml`

`conda activate esperanto`

#### Running App
`cd Esperanto`

`python3.8 u_app.py` 

or 

`python3.10 u_app.py`


## Or Install with Docker 

### Requirements 

- Docker <https://docs.docker.com/engine/install/>

### Install and Running App
`cd Esperanto`

`docker build -t esperanto . `

`docker run -it --env="DISPLAY" --net=host -v ${PWD}/:/app/ esperanto:latest //from a VM`


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

Understanding each field of StatisticReport:

- QTreshold_1 : Threshold of quality of filtration of reads \
- Depth_1 : Number of reads filtered \
- Identity_Percent_2 : Percentage identity threshold by quality to form the main cluster in order to create a consensus with vsearch \
- Depth_2: Number of reads inside the main cluster in order to create a consensus sequence \ 
- Mean_read_depth: The average of the depth of the multi-alignment assembling by samtools \
- Breath of coverage: The coverage of the multi-alignment assembling by samtools \
- Length: The length of final consensus \
- Sequence fasta: Sequence consensus (without degenerated nucleotides) 

 



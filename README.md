# Esperanto

## Description 

An app to analyze Sequencing Data from ONT technology.

This app allows browsing or dragging and dropping FASTQ files. The user may choose a few parameters to filter the data. Then, the program is launched. 
The workflow in the background is built using **Snakemake**.
The result is a sequence Fasta consensus, with or without the variant calling, along with a statistical report.

## Install with Virtualenv (Python Environment)

### Requirements

#### Linux machine (Ubuntu 20.04, Ubuntu 22.04, Fedora)
- Python 3.8 or later
- **Virtualenv** (for creating a Python virtual environment)
- **Dependencies**:
  - Ubuntu 20.04 / 22.04 or Fedora:
    - `libgtk-3-dev`
    - `libwebkit2gtk-4.0-dev`
    - `libffi6` (Ubuntu only)
    - `libcanberra-gtk-module`
    - `dbus-x11`
    - `build-essential`
    - **Other development libraries** (see the detailed dependencies)

### Install Steps

1. **Install system dependencies (Ubuntu 20.04 / 22.04)**:
    ```bash
    sudo add-apt-repository universe
    sudo apt update
    sudo apt install -y libwebkitgtk-1.0-0 libffi6 libgtk-3-dev libreoffice libcanberra-gtk-module dbus-x11
    ```

    **For Fedora**:
    ```bash
    sudo dnf install -y python3-devel freeglut-devel mesa-libGL-devel \
        mesa-libGLU-devel gstreamer1-plugins-base-devel gtk3-devel \
        libjpeg-devel libnotify-devel libpng-devel SDL2-devel libSM-devel \
        libtiff-devel webkit2gtk3-devel libXtst-devel
    ```

2. **Setup Python Virtual Environment**:
   1. Install Python dependencies:
      ```bash
      python3 -m pip install --upgrade pip
      python3 -m pip install virtualenv
      ```

   2. Create a new virtual environment:
      ```bash
      python3 -m venv env_wxpython
      ```

   3. Activate the virtual environment:
      ```bash
      source env_wxpython/bin/activate
      ```

3. **Install wxPython**:
   Inside the activated environment, install **wxPython**:
   ```bash
   pip install wheel
   pip install wxPython PyYAML pandas snakemake pypubsub


### Or run `bash setup_virtualenv.sh`

### Launch GUI on a remote VM 
`sudo nano /etc/ssh/sshd_config`

- set the field X11forwarding to 'yes':

X11Forwarding yes

- save 

- Install xauth on ssh server
`sudo apt update && sudo apt install --assume-yes xauth `# Ubuntu and other Debian-based distribution
`sudo dnf install --assumeyes xorg-x11-xauth ` # CentOS and other Red Hat based distributions

- Launch SSH connection 
`ssh -X user@IP`

- Set in the end of ~/.bashrc 
export DISPLAY=:0

#### Running App
`cd Esperanto`

`python u_app.py` 


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

 



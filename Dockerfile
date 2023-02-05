FROM python:3.8.1-buster


RUN apt-get update
RUN apt-get install -y build-essential && \
	apt-get install -y wget 

RUN apt-get install -y libgtk-3-dev libreoffice
RUN pip install attrdict
RUN pip install -f https://extras.wxpython.org/wxPython4/extras/linux/gtk3/debian-9 wxPython
RUN pip install PyYAML snakemake pandas 

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
	/bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda install -c conda-forge mamba

COPY . /app/
WORKDIR /app

VOLUME ${PWD}/data /data/
CMD ["python3.8", "u_app.py"])

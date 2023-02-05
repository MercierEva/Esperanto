FROM python:3.8.1-buster


RUN apt-get update
RUN apt-get install -y build-essential libgtk-3-dev libreoffice
RUN pip install attrdict
RUN pip install -f https://extras.wxpython.org/wxPython4/extras/linux/gtk3/debian-9 wxPython
RUN pip install PyYAML snakemake pandas 

COPY . /app/
WORKDIR /app

VOLUME ${PWD}/data /data/
CMD ["python3.8", "u_app.py"])

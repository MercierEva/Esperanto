# -*- coding: utf-8 -*-
import gzip
import time
import pandas as pd
import subprocess
import numpy as np

def filesCheck(wildcards):
    file = checkpoints.filtration.get(**wildcards).output[0]
    print(file)
    print(count_seq(file))
    if count_seq(file) > int(config["params"]["coverage"]):           
        print("pass")
        pass_list=[]
        pass_list = expand(config["folder"]+"06_stats/report_{sample}_final.tsv", sample= wildcards.sample ) + expand(config["folder"]+"05_porechop/{sample}_ont.fasta", sample =wildcards.sample ) 
        return pass_list
    else:      
        return expand(config["folder"]+"06_stats/Temp/{sample}_state.temp", sample=wildcards.sample)


configfile: "config_wf.yaml"
include: "rule_workflow_otu.smk"


onsuccess:
    print("Le workflow est fini; vous pouvez récupérer les fasta ci-dessous")
    shell("rm .snakemake/log/*")

onerror:
    print("Le workflow ne s'est pas déroulé comme prévu...")


rule all:
    input:
        expand(config["folder"]+"06_stats/Temp/{sample}_finished.temp", sample=config["samples"])

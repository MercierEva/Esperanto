import os

checkpoint filtration:
    input :  
        ech = lambda wildcards: config["samples"][wildcards.sample]
    output : 
        config["folder"]+"/01_nanofilt/PASS/{sample}_filt.fastq.gz",
    params : 
        min_length = config["params"]["filtration"]["min_length"],
        max_length = config["params"]["filtration"]["max_length"],
        rd = config["params"]["filtration"]["readtype"], 
        cov = config["params"]["coverage"], 
        folder = config["folder"]
    conda : 
        "envs/nanofilt.yaml"
    message : 
        "The filtration between {params.min_length} and {params.max_length} on a variable quality score according to the samples is launched."
    shell: """
        bash scripts/sh_scripts/run_filtration.sh {input.ech} {params.min_length} {params.max_length} {params.rd} {params.cov} {wildcards.sample} {output[0]} {params.folder} || true
    """

rule alt_setting :
    input:
        config["folder"]+"01_nanofilt/PASS/{sample}_filt.fastq.gz"
    output:
        config["folder"]+"06_stats/Temp/{sample}_state.temp"
    shell: """
        echo "not pass filtration step" > {output}
    """

rule fastq_to_fasta :
    input:
        config["folder"]+"/01_nanofilt/PASS/{sample}_filt.fastq.gz"
    output:
        config["folder"]+"/01_nanofilt/PASS/{sample}_filt.fasta"
    conda: 
        "envs/pandas.yaml"
    message : "Converting fastq.gz to fasta"
    shell: """
        zcat {input} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr '\\t' '\\n' > {output}
    """

rule cluster_to_consensus :
    input: 
        config["folder"]+"/01_nanofilt/PASS/{sample}_filt.fasta",     
    output: 
        config["folder"]+"/02_vsearch/PASS/consensus_{sample}.fasta",
    conda :
        "envs/vsearch.yaml"
    threads : config["params"]["threading"]
    params : 
        cluster_dir = config["folder"]+"/02_vsearch/{sample}_cluster_", 
        cov = config["params"]["coverage"], 
        folder = config["folder"]
    message : 
        "Cluster formation according to a rate of similarity between reads:"
        " - Exclusion of reads that are too different"
        " - Consensus building"  
    shell : """
        bash scripts/sh_scripts/run_vsearch_otu.sh {input[0]} {threads} {wildcards.sample} {params.cluster_dir} {params.cov} {output[0]} {params.folder} || true
        """

rule rename : 
    input : 
        config["folder"]+"/02_vsearch/PASS/consensus_{sample}.fasta"
    output :
        config["folder"]+"/02_vsearch/PASS/sequence_{sample}_consensus.fasta"
    message : 
        "The consensus header of the selected (dominant) cluster is renamed."
    shell: """
        awk \'/^>/{{ split($0,a,";seqs=");  print \">{wildcards.sample}_consensus_vsearch_with_\"a[2]; next }}{{ print $0 }}\' {input} > {output} || true
        """

rule orientation_to_forward : 
    input : 
        config["folder"]+"/02_vsearch/PASS/sequence_{sample}_consensus.fasta"
    output : 
        config["folder"]+"/02_vsearch/PASS/sequence_{sample}_consensus_forward.fasta"
    message : 
        "Reorientation in Sense if the consensus of the strand is in Anti-sense."
    params : 
        amorce = config["params"]["amorce_Reverse"], 
        folder=config["folder"]
    shell :
        "python scripts/py_scripts/match_primer_deg.py {input} {params.folder} {wildcards.sample} {output} {params.amorce}"

rule correct_and_polish :
    input : 
        config["folder"]+"/02_vsearch/PASS/sequence_{sample}_consensus_forward.fasta",
    threads : config["params"]["threading"]
    conda : 
        "envs/medaka.yaml"
    params : 
        cluster_dir = config["folder"]+"/02_vsearch/{sample}_cluster_", 
        model = config["params"]["model"],
        folder = config["folder"]
    output :
        config["folder"]+"/03_medaka/{sample}_medaka_RN.fasta",
        report=temp(config["folder"]+"/05_stats/report_{sample}_complementary2.tsv")
    message : 
        "Correction/Polishing by Medaka, be careful with the choice of the model..."
    shell: """ 
        bash scripts/sh_scripts/run_medaka_otu.sh {params.cluster_dir} {input} {threads} {params.model} {output[0]} {wildcards.sample} {output.report} {params.folder}
    """

rule delete_adapters : 
    input : 
        config["folder"]+"/03_medaka/{sample}_medaka_RN.fasta",
        report=config["folder"]+"/05_stats/report_{sample}_complementary2.tsv"
    output : 
        config["folder"]+"/04_porechop/{sample}_ont.fasta",
        report=config["folder"]+"/05_stats/report_{sample}_final.tsv" 
    message : 
        "Cutting of the adapters."
    conda : 
        "envs/porechop.yaml"
    shell :"""
        porechop -i {input[0]} -o {output[0]} || true  
        python scripts/py_scripts/report_stats_final_otu.py {output[0]} {input.report} {output.report} || true
    """


rule aggregated: 
    input:
        filesCheck
    output: 
        config["folder"]+"05_stats/Temp/{sample}_finished.temp"
    run:
        if len(input) > 1 :
            with open(config["folder"]+"All_fastas.fasta", "a+") as filefinal:
                with open(input[1], "r") as file :
                    filefinal.write(file.read())
                       
            with open(config["folder"] + "ReportStatistics.tsv", "a+") as filefinal:
                with open(input[0], "r") as file :
                    print(filefinal.tell())
                    if filefinal.tell() != 0 :
                        filefinal.write(file.readlines()[1])
                    else:
                        filefinal.write(file.read())

            df = pd.read_csv(config["folder"]+"ReportStatistics.tsv", sep='\t')      
            pivot_df = pd.pivot_table(df, index='sequencies', columns='sample', values='Depth_2', aggfunc=np.sum, fill_value=0)
            samples_list = list(pivot_df.columns)
            new_df = pivot_df.sort_values(by=samples_list, ascending=len(samples_list) * [False])

            new_df.to_csv(config["folder"]+"Table_OTU.tsv", sep='\t', encoding='utf-8')
                
            with open(output[0], "w") as filefinal:
                filefinal.write("finished")
             
        else:
            with open(output[0], "w") as filefinal:
                filefinal.write("finished")
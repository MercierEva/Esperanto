checkpoint filtration :
    input:
        ech=lambda wildcards: config["samples"][wildcards.sample] 
    output : 
        files=config["folder"]+"01_nanofilt/PASS/{sample}_filt.fastq.gz"
    params : 
        min_length = config["params"]["filtration"]["min_length"],
        max_length = config["params"]["filtration"]["max_length"],
        rd = config["params"]["filtration"]["readtype"], 
        content = config["params"]["filtration"]["min_content"],
        folder = config["folder"]
    threads : config["params"]["threading"]
    conda: "envs/nanofilt.yaml"  
    message : 
        "The filtration between {params.min_length} and {params.max_length} on a variable quality score according to the samples is launched."
    shell: """
        python workflow/scripts/py_scripts/run_filtration.py {input.ech} {params.min_length} {params.max_length} {params.rd} {wildcards.sample} {output.files} {params.folder} {params.content}
    """

rule alt_setting :
    input:
        config["folder"]+"01_nanofilt/PASS/{sample}_filt.fastq.gz"
    output:
        config["folder"]+"07_stats/Temp/{sample}_state.temp"
    shell: """
        echo "not pass filtration step" > {output}
    """

rule fastq_to_fasta :
    input:
        config["folder"]+"01_nanofilt/PASS/{sample}_filt.fastq.gz"
    output:
        config["folder"]+"01_nanofilt/PASS/{sample}_filt.fasta"
    message : "Converting fastq.gz to fasta"
    shell: """
        zcat {input[0]} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr '\\t' '\\n' > {output[0]}	
    """

rule cluster_to_consensus :
    input: 
        fasta=config["folder"]+"01_nanofilt/PASS/{sample}_filt.fasta"
    output: 
        config["folder"]+"02_vsearch/PASS/consensus_{sample}.fasta"
    conda :
        "envs/vsearch.yaml"
    threads : config["params"]["threading"]
    params : 
        cluster_dir = config["folder"]+"02_vsearch/{sample}_cluster_", 
        folder = config["folder"]
    message : 
        "Cluster formation according to a rate of similarity between reads:"
        " - Exclusion of reads that are too different"
        " - Consensus building"   
    shell : """
        bash workflow/scripts/sh_scripts/run_vsearch.sh {input.fasta} {threads} {wildcards.sample} {params.cluster_dir} {output[0]} {params.folder}
        """


rule rename : 
    input : 
        config["folder"]+"02_vsearch/PASS/consensus_{sample}.fasta"
    output :
        config["folder"]+"02_vsearch/PASS/sequence_{sample}_consensus.fasta"
    message : 
        "The consensus header of the selected (dominant) cluster is renamed."
    shell: r"""
        awk '/^>/{{print ">{wildcards.sample}_consensus_vsearch";next}}{{print $1}}' {input[0]} | head -2 > {output[0]} || true
        """


rule orientation_to_forward : 
    input : 
        config["folder"]+"02_vsearch/PASS/sequence_{sample}_consensus.fasta"
    output : 
        config["folder"]+"02_vsearch/PASS/sequence_{sample}_consensus_forward.fasta"
    message : 
        "Reorientation in Sense if the consensus of the strand is in Anti-sense."
    params : 
        amorce = config["params"]["amorce_Reverse"]
    shell :"""
        cat {input[0]} | while read L
        do 
            if [[ $L =~ ^'>' ]] 
            then 
                echo $L > {output[0]} 
            else if [[ $L =~ \"{params.amorce}\" ]] 
            then 
                echo $L | tr ATCG TAGC | rev >>  {output[0]} 
            else 
                echo $L >>  {output}
                fi
            fi 
        done || true
    """

rule delete_adapters : 
    input : 
        config["folder"]+"02_vsearch/PASS/sequence_{sample}_consensus_forward.fasta"
    output :
    	config["folder"]+"03_porechop/{sample}_ont.fasta"
    message : 
        "Cutting of the adapters."
    conda : 
        "envs/porechop.yaml"
    shell : """
        porechop -i {input[0]} -o {output[0]} 
    """

rule multialignment : 
    input:
        ref=config["folder"]+"03_porechop/{sample}_ont.fasta"
    output:
        config["folder"]+"04_multialignment/{sample}_MSA.sam"
    params:
        cluster_dir = config["folder"]+"02_vsearch/{sample}_cluster_"
    threads : config["params"]["threading"]    
    conda: "envs/ngmlr.yaml"
    message: "Multialignment"
    shell:"""
        bash workflow/scripts/sh_scripts/run_ngmlr.sh {threads} {input.ref} {params.cluster_dir} {output}
    """

rule samtools view :
    input : 
        samfile=config["folder"]+"04_multialignment/{sample}_MSA.sam"
    output : 
        config["folder"]+ "04_multialignment/{sample}_MSA.bam"
    conda: "envs/samtools.yaml"
    params : folder = config["folder"]
    message:"sam_to_bam"
    shell:"""
        samtools view -S -b {input.samfile} > {output}
       """

rule samtools_sort :
    input: config["folder"]+ "04_multialignment/{sample}_MSA.bam"
    output : config["folder"]+ "04_multialignment/{sample}_MSA.sorted.bam"
    conda: "envs/samtools.yaml"
    params : folder = config["folder"]
    message:"sorting bamfile"    
    shell:"""
        samtools sort {input} -o {output}
        """

rule samtools_mpileup :
    input: 
        bam = config["folder"]+ "04_multialignment/{sample}_MSA.sorted.bam"
    output : config["folder"]+"04_multialignment/{sample}.pileup"
    conda: "envs/samtools.yaml"
    params : folder = config["folder"]
    message:"bulding pileup file"
    shell:"""
        samtools index {input.bam}
        samtools mpileup -f "{params.folder}03_porechop/{wildcards.sample}_ont.fasta" {input.bam} > {output[0]}
        """
    

rule variant_calling:
    input : 
        config["folder"]+"04_multialignment/{sample}.pileup"
    output :
        config["folder"]+"05_varscan/{sample}_CNS.tsv"
    conda : "envs/varscan.yaml"
    params: 
        folder = config["folder"]
    message : "Variant calling"
    shell:"""
        bash workflow/scripts/sh_scripts/run_varscan.sh {input} {params.folder}05_varscan/{wildcards.sample}_SNP.tsv {params.folder}05_varscan/{wildcards.sample}_INDEL.tsv {output} {wildcards.sample} {params.folder} 
    """

rule consensus : 
    input :
        config["folder"]+"05_varscan/{sample}_CNS.tsv"
    output:
        config["folder"]+"06_seq_inform/{sample}_SCAN.fa"
    message :
        "Sequence reconstruction with integrated variants."
    shell: r"""
        awk -F"\\t" 'BEGIN{{print ">{wildcards.sample}"}}
        {{
            if(NR>2){{
                gsub(/[\*\+-].*\/[\+-]/, "", $4)
                if($4 ~ /\/\+/ ){{
                    printf $3$4}}
                else if($4 ~ /\/-/ ){{
                    printf $3
                    NR++}}
                else {{
                    printf $4}}
                }}
        }}END{{print "\\n"}}' {input} > {output} || true    
        """

rule report :
    input :
        config["folder"]+"06_seq_inform/{sample}_SCAN.fa"
    output :
        config["folder"]+"07_stats/report_{sample}_final.tsv"
    params:
        folder=config["folder"]
    conda: "envs/report.yaml"
    message : 
        "Creating ReportStatistics" 
    shell : """
        python workflow/scripts/py_scripts/report_stats.py {wildcards.sample} {params.folder} {output[0]}
    """

rule aggregated: 
    input:
        filesCheck
    output: 
        config["folder"]+"07_stats/Temp/{sample}_finished.temp"
    threads:
        config["params"]["threading"]
    run:
        if input[0].find("final") != -1:
            with open(config["folder"]+"All_fastas.fasta", "a+") as filefinal:
                with open(config["folder"]+"03_porechop/"+wildcards.sample+"_ont.fasta", "r") as file :
                    filefinal.write(file.read())
            
            with open(config["folder"]+ "All_Consensus_fastas.fasta", "a+") as filefinal:
                with open(config["folder"]+"06_seq_inform/"+wildcards.sample+"_SCAN.fa","r") as file:
                    filefinal.write(file.read())
                
            with open(config["folder"] + "StatisticReport.tsv", "a+") as filefinal:
                with open(input[0], "r") as file :
                    if filefinal.tell() == 0 :
                        filefinal.write(file.read())
                    else:
                        filefinal.write(file.readlines()[1])
                
            with open(output[0], "w") as filefinal:
                filefinal.write("finished")
              
        else :
            with open(output[0], "w") as filefinal:
                filefinal.write("finished")

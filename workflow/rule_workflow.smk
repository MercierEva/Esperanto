

def count_seq(input_file):
    proc = subprocess.Popen("gunzip -c " + input_file + " | wc -l ", shell=True, stdout=subprocess.PIPE) 
    out, err = proc.communicate()
    count= int(out.decode('utf-8'))
    return count 

checkpoint filtration :
    input:
        ech=lambda wildcards: config["samples"][wildcards.sample] 
    output : 
        files=config["folder"]+"01_nanofilt/PASS/{sample}_filt.fastq.gz"
    params : 
        min_length = config["params"]["filtration"]["min_length"],
        max_length = config["params"]["filtration"]["max_length"],
        rd = config["params"]["filtration"]["readtype"], 
        cov = config["params"]["coverage"], 
        folder = config["folder"]
    threads : config["params"]["threading"]
    conda: "envs/nanofilt.yaml"  
    message : 
        "The filtration between {params.min_length} and {params.max_length} on a variable quality score according to the samples is launched."
    shell: """
        bash scripts/sh_scripts/run_filtration.sh {input.ech} {params.min_length} {params.max_length} {params.rd} {params.cov} {wildcards.sample} {output.files} {params.folder} 
    """


rule set_of_fastq : 
    input:
        config["folder"]+"01_nanofilt/PASS/{sample}_filt.fastq.gz"
    output:
        report=temp(config["folder"]+"06_stats/report_{sample}_nanofilt.tsv"),
        files=config["folder"]+"02_filtlong/{sample}_filt_set.fastq.gz"
    params:
        cov = config["params"]["coverage"],
        folder=config["folder"]
    conda:
        "envs/filtlong.yaml"
    message : "Extracting best of reads."
    shell:"""
	python scripts/py_scripts/report_stats.py {wildcards.sample} {output.report} {input} {params.folder} 
	bash scripts/sh_scripts/run_filtlong.sh {input} {output.files} {params.cov}        
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
        config["folder"]+"02_filtlong/{sample}_filt_set.fastq.gz"
    output:
        config["folder"]+"02_filtlong/PASS/{sample}_filt_set.fasta"
    message : "Converting fastq.gz to fasta"
    shell: """
        zcat {input[0]} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr '\\t' '\\n' > {output[0]}	
    """

rule cluster_to_consensus :
    input: 
        fasta=config["folder"]+"02_filtlong/PASS/{sample}_filt_set.fasta"
    output: 
        config["folder"]+"03_vsearch/PASS/consensus_{sample}.fasta"
    conda :
        "envs/vsearch.yaml"
    threads : config["params"]["threading"]
    params : 
        cluster_dir = config["folder"]+"03_vsearch/{sample}_cluster_", 
        cov = config["params"]["coverage"], 
        folder = config["folder"]
    message : 
        "Cluster formation according to a rate of similarity between reads:"
        " - Exclusion of reads that are too different"
        " - Consensus building"   
    shell : """
        bash scripts/sh_scripts/run_vsearch.sh {input.fasta} {threads} {wildcards.sample} {params.cluster_dir} {params.cov} {output[0]} {params.folder}
        """

rule multialignment : 
    input:
        ref=config["folder"]+"05_porechop/{sample}_ont.fasta",
        reads=config["folder"]+"02_filtlong/{sample}_filt_set.fastq.gz"
    output:
        config["folder"]+"07_multialignment/{sample}_MSA.sam"
    threads : config["params"]["threading"]    
    conda: "envs/VC.yaml"
    message: "Multialignment"
    shell:"""
        ngmlr -x ont -t {threads} -r {input.ref} -q {input.reads} -o {output}  
    """

rule samtools :
    input : 
        ref=config["folder"]+"05_porechop/{sample}_ont.fasta",
        samfile=config["folder"]+"07_multialignment/{sample}_MSA.sam"
    output : 
        config["folder"]+"07_multialignment/{sample}.pileup"
    conda: "envs/VC.yaml"
    params : folder = config["folder"]
    message:"bulding pileup file"
    shell:"""
        samtools view -S -b {input.samfile} > {params.folder}07_multialignment/{wildcards.sample}_MSA.bam
        samtools sort {params.folder}07_multialignment/{wildcards.sample}_MSA.bam -o {params.folder}07_multialignment/{wildcards.sample}_MSA.sorted.bam
        samtools index {params.folder}07_multialignment/{wildcards.sample}_MSA.sorted.bam
        samtools mpileup -f {input.ref} {params.folder}07_multialignment/{wildcards.sample}_MSA.sorted.bam > {output}
    """

rule variant_calling:
    input : 
        config["folder"]+"07_multialignment/{sample}.pileup"
    output :
        SNP=config["folder"]+"08_varscan/{sample}_SNP.tsv",
        INDEL=config["folder"]+"08_varscan/{sample}_INDEL.tsv",
        CNS=config["folder"]+"08_varscan/{sample}_CNS.tsv"
    conda : "envs/VC.yaml"
    params: folder=config["folder"]
    message : "Variant calling"
    shell:"""
        bashrun_varscan.sh {input} {output[0]} {output[1]} {output[2]} {wildcards.sample} {params.folder}
    """


rule consensus : 
    input :
        config["folder"]+"08_varscan/{sample}_CNS.tsv"
    output:
        config["folder"]+"09_seq_inform/{sample}_SCAN.fa"
    message :
        "Sequence reconstruction with integrated variants. "
    shell : """    
        awk -F\"\\t\" \'BEGIN{{print \">{wildcards.sample}\"}}
        {{
            if(NR>2){{
                {{gsub(/[\*\+-].*\/[\+-]/, \"\", $4)}}
                if($4 ~ /\/\+/ ){{
                    printf $3$4 }}
                else if($4 ~ /\/-/ ){{
                    printf $3
                    NR++ }}
                else {{
                    printf $4}}
                }}
            }}END{{ print \"\\n\" }}\' {input} > {output} || true    
    """

rule rename : 
    input : 
        config["folder"]+"03_vsearch/PASS/consensus_{sample}.fasta"
    output :
        config["folder"]+"03_vsearch/PASS/sequence_{sample}_consensus.fasta"
    message : 
        "The consensus header of the selected (dominant) cluster is renamed."
    shell: """
        awk \'/^>/{{ print \">{wildcards.sample}_consensus_vsearch\" ; next }}{{ print $1 }}\' {input[0]} | head -2 > {output[0]} || true
        """

rule orientation_to_forward : 
    input : 
        config["folder"]+"03_vsearch/PASS/sequence_{sample}_consensus.fasta"
    output : 
        config["folder"]+"03_vsearch/PASS/sequence_{sample}_consensus_forward.fasta"
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

rule correct_and_polish :
    input : 
        config["folder"]+"03_vsearch/PASS/sequence_{sample}_consensus_forward.fasta"
    threads : config["params"]["threading"]
    conda : 
        "envs/medaka.yaml"
    params : 
        cluster_dir = config["folder"]+"03_vsearch/{sample}_cluster_", 
        model = config["params"]["model"], 
    output :
        directory(config["folder"]+"04_medaka/{sample}_Medaka")
    message : 
        "Correction/Polishing by Medaka, be careful with the choice of the model..."
    shell: """ 
        bash scripts/sh_scripts/run_medaka.sh {params.cluster_dir} {input[0]} {threads} {output[0]} {params.model} || true
	"""

rule rename2 :
    input :
        rules.correct_and_polish.output,
	    rules.set_of_fastq.output.report
    output :
        config["folder"]+"04_medaka/{sample}_medaka_RN.fasta",
        temp(config["folder"]+"06_stats/report_{sample}_complementary2.tsv")
    conda:
        "envs/medaka.yaml"
    params:
        folder=config["folder"]
    message : 
        "The final fasta file is renamed." 
    shell : """
        awk \'/^>/{{ print \">{wildcards.sample}_medaka\" ; next }}{{print $1}}\' {input[0]}/consensus.fasta > {output[0]} || true
        python scripts/py_scripts/report_stats_complementary2.py {wildcards.sample} {input[0]} calls_to_draft.bam consensus.fasta {input[1]} {output[1]} {params.folder} || true 
    """

rule delete_adapters : 
    input : 
        config["folder"]+"04_medaka/{sample}_medaka_RN.fasta",
        config["folder"]+"06_stats/report_{sample}_complementary2.tsv"
    output :
        config["folder"]+"06_stats/report_{sample}_final.tsv",
	config["folder"]+"05_porechop/{sample}_ont.fasta"
    message : 
        "Cutting of the adapters."
    conda : 
        "envs/porechop.yaml"
    shell : """
        porechop -i {input[0]} -o {output[1]} 
        python scripts/py_scripts/report_stats_final.py {wildcards.sample} {input[0]} {input[1]} {output[0]}
    """


rule aggregated: 
    input:
        filesCheck
    output: 
        config["folder"]+"06_stats/Temp/{sample}_finished.temp"
    run:
        if len(input) > 1 :
            with open(config["folder"]+"All_fastas.fasta", "a+") as filefinal:
                with open(input[1], "r") as file :
                    filefinal.write(file.read())
            
            with open(config["folder"]+ "All_Consensus_fastas.fasta", "a+") as filefinal:
                with open(input[2],"r") as file:
                    filefinal.write(file.read())
                
            with open(config["folder"] + "ReportStatistics.tsv", "a+") as filefinal:
                with open(input[0], "r") as file :
                    print(filefinal.tell())
                    if filefinal.tell() != 0 :
                        filefinal.write(file.readlines()[1])
                    else:
                        filefinal.write(file.read())
                
            with open(output[0], "w") as filefinal:
                filefinal.write("finished")

                
        else:
            with open(output[0], "w") as filefinal:
                filefinal.write("finished")

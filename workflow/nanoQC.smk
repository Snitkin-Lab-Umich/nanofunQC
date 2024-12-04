# Author: Ali Pirani and Dhatri Badri 
configfile: "config/confi.yaml"

import pandas as pd
import os
import numpy as np

samples_df = pd.read_csv(config["samples"])
BARCODE = list(samples_df['barcode_id'])
PREFIX = config["prefix"]

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

# Snakefile to generate summary report
include: "nanoQC_summary.smk"

rule all:
    input:
        trimmed = expand("results/{prefix}/filtlong/{barcode}/{barcode}.trimmed.fastq.gz", barcode=BARCODE, prefix=PREFIX),
	    nanoplot = expand("results/{prefix}/nanoplot/{barcode}/{barcode}_preqcNanoPlot-report.html", barcode=BARCODE, prefix=PREFIX),
        flye_assembly = expand("results/{prefix}/flye/{barcode}/{barcode}_flye.fasta", barcode=BARCODE, prefix=PREFIX),
        flye_circ_assembly = expand("results/{prefix}/flye/{barcode}/{barcode}_flye_circ.fasta", barcode=BARCODE, prefix=PREFIX),
        medaka_out = expand("results/{prefix}/medaka/{barcode}/{barcode}_medaka.fasta", barcode=BARCODE, prefix=PREFIX),
        prokka_out = expand("results/{prefix}/prokka/{barcode}/{barcode}_medaka.gff", barcode=BARCODE, prefix=PREFIX),
        quast_out = expand("results/{prefix}/quast/{barcode}/{barcode}_flye/report.txt", barcode=BARCODE, prefix=PREFIX),
        busco_out = expand("results/{prefix}/busco/{barcode}/{barcode}.medaka/busco_medaka.txt", barcode=BARCODE, prefix=PREFIX),
        mlst_out = expand("results/{prefix}/mlst/{barcode}/report.tsv", barcode=BARCODE, prefix=PREFIX),
        skani_out = expand("results/{prefix}/skani/{barcode}/{barcode}_skani_output.txt", barcode=BARCODE, prefix=PREFIX),
        nanoplot_report = expand("results/{prefix}/{prefix}_report/{prefix}_nanoplot_results.csv", barcode=BARCODE, prefix=PREFIX),
        summary_report = expand("results/{prefix}/{prefix}_report/{prefix}_report.csv", prefix=PREFIX)
        
rule filtlong:
    input:
        longreads = config["long_reads"] + "/{barcode}"
    output:
        trimmed = "results/{prefix}/filtlong/{barcode}/{barcode}.trimmed.fastq.gz"
    #log:
    #    "logs/{prefix}/filtlong/{barcode}/{barcode}.log"   
    singularity:
        "docker://staphb/filtlong:0.2.1"
    #envmodules:
    #    "Bioinformatics",
    #    "filtlong"
    shell:
        "filtlong --min_length 1000 --keep_percent 95 {input.longreads}/*.fastq.gz  --mean_q_weight 10 --target_bases 500000000 | gzip > {output.trimmed}"
        
rule nanoplot:
    input:
        longreads = config["long_reads"] + "/{barcode}",
        trimmed = lambda wildcards: expand(f"results/{wildcards.prefix}/filtlong/{wildcards.barcode}/{wildcards.barcode}.trimmed.fastq.gz"),
    output:
        nanoplot_preqc = "results/{prefix}/nanoplot/{barcode}/{barcode}_preqcNanoPlot-report.html"
    #log:
    #    "logs/{prefix}/nanoplot/{barcode}/{barcode}.log"
    params:
        outdir="results/{prefix}/nanoplot/{barcode}",
        prefix="{barcode}",
    singularity:
        "docker://staphb/nanoplot:1.42.0"
    #envmodules:
    #    "Bioinformatics",
    #    "nanoplot"
    shell:
        """
        cat {input.longreads}/*.fastq.gz > /tmp/{params.prefix}.gz && 
        NanoPlot -o {params.outdir} -p {params.prefix}_preqc --tsv_stats --info_in_report --N50 --title {params.prefix}_preqc --fastq /tmp/{params.prefix}.gz && 
        NanoPlot -o {params.outdir} -p {params.prefix}_postqc --tsv_stats --info_in_report --N50 --title {params.prefix}_postqc --fastq {input.trimmed} && rm /tmp/{params.prefix}.gz 
        """
        
rule flye:
    input:
        trimmed = lambda wildcards: f"results/{wildcards.prefix}/filtlong/{wildcards.barcode}/{wildcards.barcode}.trimmed.fastq.gz"
    output:
        assembly = "results/{prefix}/flye/{barcode}/{barcode}_flye.fasta",
    params:
        assembly_dir = "results/{prefix}/flye/{barcode}/",
        size = config["genome_size"],
        threads = config["threads"],
        #flye_options = config["flye_options"],
        prefix = "{barcode}",
    #log:
    #    "logs/{prefix}/flye/{barcode}/{barcode}_flye.log" # Flye has its own log
    singularity:
        "docker://staphb/flye:2.9.4"
    #envmodules:
    #    "Bioinformatics",
    #    "flye"
    shell:   
        """
        flye --nano-hq {input.trimmed} -g {params.size} -o {params.assembly_dir} -t {params.threads} --debug && 
        cp {params.assembly_dir}/assembly.fasta {params.assembly_dir}/{params.prefix}_flye.fasta 
        """
     
rule flye_add_circ:
    input:
        flye_assembly = "results/{prefix}/flye/{barcode}/{barcode}_flye.fasta",
    output:
        assembly = "results/{prefix}/flye/{barcode}/{barcode}_flye_circ.fasta",
    params:
        assembly_dir = "results/{prefix}/flye/{barcode}/",
        #size = config["genome_size"],
        #threads = config["threads"],
        #flye_options = config["flye_options"],
        prefix = "{barcode}",
    #log:
    #    "logs/{prefix}/flye/{barcode}/{barcode}_flye_add_circ.log",
    run:
        shell("cp {params.assembly_dir}/{params.prefix}_flye.fasta {params.assembly_dir}/{params.prefix}_flye_circ.fasta")
        assembly_info = pd.read_csv(f"{params.assembly_dir}/assembly_info.txt", sep='\t', header=0)
        assembly_info["circular"] = np.where(assembly_info["circ."] == "Y", "true", "false")
        flye_assembly_circ = output.assembly
        for index, row in assembly_info.iterrows():
            circular = f"{row['#seq_name']};circular={row['circular']}"
            shell(f"sed -i 's/\\<{row['#seq_name']}\\>/{circular}/g' {flye_assembly_circ}")

rule medaka:
    input:
        trimmed = lambda wildcards: expand(f"results/{wildcards.prefix}/filtlong/{wildcards.barcode}/{wildcards.barcode}.trimmed.fastq.gz"),
        flye_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/flye/{wildcards.barcode}/{wildcards.barcode}_flye_circ.fasta"),
    output:
        medaka_out = f"results/{{prefix}}/medaka/{{barcode}}/{{barcode}}_medaka.fasta",
    params:
        medaka_out_dir = "results/{prefix}/medaka/{barcode}",
        threads = config["threads"],
        prefix = f"{{barcode}}",
    #log:
    #    "logs/{prefix}/medaka/{barcode}/{barcode}.log"
    singularity:
        "docker://staphb/medaka:1.2.0"
    #envmodules:
    #    "Bioinformatics",
    #    "medaka",
    #    "bcftools"
    shell:
        """
        medaka_consensus -i {input.trimmed} -d {input.flye_assembly} -o {params.medaka_out_dir} -t {params.threads} -m r941_min_high_g303 &&
        cp {params.medaka_out_dir}/consensus.fasta {params.medaka_out_dir}/{params.prefix}_medaka.fasta 
        """ 
        
rule prokka:
    input:
        medaka = f"results/{{prefix}}/medaka/{{barcode}}/{{barcode}}_medaka.fasta"
    output:
        medaka_annotation = f"results/{{prefix}}/prokka/{{barcode}}/{{barcode}}_medaka.gff",
    params:
        threads = config["ncores"],
        prefix = "{barcode}",
        #options = config["prokka_options"],
        prokka_dir = directory("results/{prefix}/prokka/{barcode}"),
    log:
        "logs/{prefix}/prokka/{barcode}/{barcode}.log"
    singularity:
        "docker://staphb/prokka:1.14.6"
    #envmodules:
    #    "Bioinformatics",
    #    "prokka"
    shell:
        "prokka --force --kingdom Bacteria --rfam --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_medaka {input.medaka} &>{log}"

rule quast:
    input:
        flye_assembly = f"results/{{prefix}}/flye/{{barcode}}/{{barcode}}_flye_circ.fasta",
        medaka_out = f"results/{{prefix}}/medaka/{{barcode}}/{{barcode}}_medaka.fasta",
    output:
        quast_out_flye = f"results/{{prefix}}/quast/{{barcode}}/{{barcode}}_flye/report.txt",
        quast_out_medaka = f"results/{{prefix}}/quast/{{barcode}}/{{barcode}}_medaka/report.txt",
    params:
        threads = config["ncores"],
        quast_dir = directory("results/{prefix}/quast/{barcode}/{barcode}"),
    log:
        "logs/{prefix}/quast/{barcode}/{barcode}.log"
    singularity:
        "docker://staphb/quast:5.2.0"
    #envmodules:
    #    "Bioinformatics",
    #    "quast"
    shell:
        """
        quast.py {input.flye_assembly} -o {params.quast_dir}_flye &&
        quast.py {input.medaka_out} -o {params.quast_dir}_medaka &>{log}
        """   
  
rule busco:
    input:
        medaka_assembly = f"results/{{prefix}}/medaka/{{barcode}}/{{barcode}}_medaka.fasta",
        flye_assembly = f"results/{{prefix}}/flye/{{barcode}}/{{barcode}}_flye_circ.fasta",
    output:
        busco_flye_medaka_out = f"results/{{prefix}}/busco/{{barcode}}/{{barcode}}.medaka/busco_medaka.txt", 
        busco_flye_assembly_out = f"results/{{prefix}}/busco/{{barcode}}/{{barcode}}.flye_assembly/busco_flye_assembly.txt",     
    params:
        busco_outpath = f"results/{{prefix}}/busco/{{barcode}}/{{barcode}}",
        medaka_busco_out = f"short_summary.specific.bacteria_odb10.{{barcode}}.medaka.txt",
        flye_assembly_busco_out = f"short_summary.specific.bacteria_odb10.{{barcode}}.flye_assembly.txt",
        threads = config["ncores"],
    #log:
    #    "logs/{prefix}/busco/{barcode}/{barcode}.log" # BUSCO has it own logs folder
    singularity:
        "docker://staphb/busco:5.7.1-prok-bacteria_odb10_2024-01-08"
    #envmodules:
    #    "Bioinformatics",
    #    "busco"
    shell:
        """ 
        busco -f -i {input.medaka_assembly} -m genome -l bacteria_odb10 -o {params.busco_outpath}.medaka && 
        cp {params.busco_outpath}.medaka/{params.medaka_busco_out} {params.busco_outpath}.medaka/busco_medaka.txt &&    
        
        busco -f -i {input.flye_assembly} -m genome -l bacteria_odb10 -o {params.busco_outpath}.flye_assembly && 
        cp {params.busco_outpath}.flye_assembly/{params.flye_assembly_busco_out} {params.busco_outpath}.flye_assembly/busco_flye_assembly.txt 
        """
     
rule mlst:
    input:
        medaka_out = lambda wildcards: expand(f"results/{wildcards.prefix}/medaka/{wildcards.barcode}/{wildcards.barcode}_medaka.fasta")
    output:
        mlst_report = f"results/{{prefix}}/mlst/{{barcode}}/report.tsv",
    log:
        "logs/{prefix}/mlst/{barcode}/{barcode}.log"
    singularity:
        "docker://staphb/mlst:2.23.0-2024-03"
    #envmodules:
    #    "Bioinformatics",
    #    "mlst"
    shell:
        "mlst {input.medaka_out} > {output.mlst_report} 2>{log}"
        
rule skani:
    input:
        medaka_out = lambda wildcards: expand(f"results/{wildcards.prefix}/medaka/{wildcards.barcode}/{wildcards.barcode}_medaka.fasta")
    output:
        skani_output = f"results/{{prefix}}/skani/{{barcode}}/{{barcode}}_skani_output.txt"
    params:
        skani_ani_db = config["skani_db"],
        threads = 4
    log:
        "logs/{prefix}/skani/{barcode}/{barcode}.log"
    singularity:
        "docker://staphb/skani:0.2.1"
    shell:
        "skani search {input.medaka_out} -d {params.skani_ani_db} -o {output.skani_output} -t {params.threads} 2>{log}"
       


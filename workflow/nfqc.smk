# Author: Ali Pirani, Dhatri Badri, and Joseph Hale 
configfile: "config/config.yaml"

import pandas as pd
import os
import numpy as np

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])
PREFIX = config["prefix"]

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")


rule all:
    input:
        summary_output = expand("results/{prefix}/multiqc/{prefix}_final_qc_summary.tsv",prefix=PREFIX)

        
rule filtlong:
    input:
        longreads = config["long_reads"] + "/{sample}.fastq.gz"
    output:
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz"
    singularity:
        "docker://staphb/filtlong:0.2.1"
    resources:
        mem_mb = 5000,
        runtime = 20
    shell:
        "filtlong --min_length 1000 --keep_percent 95 {input.longreads} --mean_q_weight 10 --target_bases 500000000 | gzip > {output.trimmed}"
        
rule nanoplot:
    input:
        longreads = config["long_reads"] + "/{sample}.fastq.gz",
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz",
    output:
        nanoplot_preqc = "results/{prefix}/nanoplot/{sample}/{sample}_preqcNanoPlot-report.html",
        nanoplot_postqc = "results/{prefix}/nanoplot/{sample}/{sample}_postqcNanoStats.txt",
    params:
        outdir="results/{prefix}/nanoplot/{sample}",
    singularity:
        "docker://staphb/nanoplot:1.42.0"
    resources:
        mem_mb = 5000,
        runtime = 30
    shell:
        """ 
        NanoPlot -o {params.outdir} -p {wildcards.sample}_preqc --tsv_stats --info_in_report --N50 --title {wildcards.sample}_preqc --fastq {input.longreads} && 
        NanoPlot -o {params.outdir} -p {wildcards.sample}_postqc --tsv_stats --info_in_report --N50 --title {wildcards.sample}_postqc --fastq {input.trimmed} 
        """
        
rule flye_assembly:
    input:
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz"
    output:
        assembly = "results/{prefix}/flye/{sample}/{sample}_flye.fasta",
    params:
        assembly_dir = "results/{prefix}/flye/{sample}/",
        size = config["genome_size"],
    singularity:
        "docker://staphb/flye:2.9.4"
    threads: 8
    resources:
        mem_mb = 8000,
        runtime = 30
    shell:   
        """
        flye --nano-hq {input.trimmed} -g {params.size} -o {params.assembly_dir} -t {threads} --debug && 
        cp {params.assembly_dir}/assembly.fasta {params.assembly_dir}/{wildcards.sample}_flye.fasta 
        """
     

rule medaka:
    input:
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz",
        flye_assembly = "results/{prefix}/flye/{sample}/{sample}_flye.fasta",
    output:
        medaka_out = "results/{prefix}/medaka/{sample}/{sample}_medaka.fasta",
    params:
        medaka_out_dir = "results/{prefix}/medaka/{sample}",
    singularity:
        "docker://staphb/medaka:2.0.1"
    threads: 6
    resources:
        mem_mb = 15000,
        runtime = lambda wildcards, attempt: 20 + ((attempt-1)*60) 
    shell:
        """
        medaka_consensus -i {input.trimmed} -d {input.flye_assembly} -o {params.medaka_out_dir} -t {threads} &&
        cp {params.medaka_out_dir}/consensus.fasta {params.medaka_out_dir}/{wildcards.sample}_medaka.fasta 
        """ 
        

rule quast:
    input:
        flye_assembly = "results/{prefix}/flye/{sample}/{sample}_flye.fasta",
        medaka_out = "results/{prefix}/medaka/{sample}/{sample}_medaka.fasta",
    output:
        quast_out_flye = "results/{prefix}/quast/{sample}/{sample}_flye/report.txt",
        quast_out_medaka = "results/{prefix}/quast/{sample}/{sample}_medaka/report.txt",
    params:
        quast_dir = "results/{prefix}/quast/{sample}/{sample}",
    log:
        "logs/{prefix}/quast/{sample}/{sample}.log"
    singularity:
        "docker://staphb/quast:5.2.0"
    resources:
        mem_mb = 1000,
        runtime = 20
    shell:
        """
        quast.py {input.flye_assembly} -o {params.quast_dir}_flye &&
        quast.py {input.medaka_out} -o {params.quast_dir}_medaka &>{log}
        """   
  
# copied from funQCD

rule auriclass:
    input:
        medaka_out = "results/{prefix}/medaka/{sample}/{sample}_medaka.fasta",
    output:
        auriclass_report = "results/{prefix}/auriclass/{sample}/{sample}_report.tsv",
    resources:
        mem_mb = 5000,
        runtime = 30
    singularity:
        "docker://quay.io/biocontainers/auriclass:0.5.4--pyhdfd78af_0"
    shell:
        "auriclass --name {wildcards.sample} -o {output.auriclass_report} {input.medaka_out}"


rule funannotate_sort:
    input:
        medaka_out = "results/{prefix}/medaka/{sample}/{sample}_medaka.fasta",
    output:
        sorted_assembly = "results/{prefix}/funannotate/{sample}/{sample}_sorted.fa"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        #sample = "{sample}"
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    resources:
        mem_mb = 5000,
        runtime=60
    shell:
        """
        funannotate clean -i {input.medaka_out} -o {params.out_dir}{wildcards.sample}_cleaned.fa
        funannotate sort -i {params.out_dir}{wildcards.sample}_cleaned.fa -o {params.out_dir}{wildcards.sample}_sorted.fa --minlen 0
        """

rule repeatmasker:
    input:
        sorted_assembly = "results/{prefix}/funannotate/{sample}/{sample}_sorted.fa"
    output:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa"
    params:
        out_dir = "results/{prefix}/repeatmasker/{sample}/",
        repeat_lib = config["funqcd_lib"] + "repeat_libraries/fungi_b8441/b8441_fungi_repeatlib.fa",
    threads: 8
    resources:
        mem_mb = 5000,
        runtime = 120,
    singularity:
        "docker://dfam/tetools:1.89.2"
    shell:
        """
        RepeatMasker -xsmall -dir {params.out_dir} -lib {params.repeat_lib} \
        results/{wildcards.prefix}/funannotate/{wildcards.sample}/{wildcards.sample}_sorted.fa -pa {threads}
        mv {params.out_dir}/{wildcards.sample}_sorted.fa.masked {params.out_dir}/{wildcards.sample}_masked.fa
        """

# this requires the RNA-seq data from Teresa, and assumes that the files are in a specific format
# specifically, that read 1 contains 'R1' in the file name, that all files are in fastq.gz format and in RF order for stranded RNA-seq, and that they can be ordered alphabetically
# This is also currently hard-coded to use 30G of memory
rule funannotate_train:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa"
    output:
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        funannotate_training_pasa_gff = "results/{prefix}/funannotate/{sample}/training/funannotate_train.pasa.gff3",
        funannotate_training_stringtie = "results/{prefix}/funannotate/{sample}/training/funannotate_train.stringtie.gtf",
        funannotate_training_transc_align = "results/{prefix}/funannotate/{sample}/training/funannotate_train.transcripts.gff3",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        rna_data_r1 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R1_' in x and 'fastq.gz' in x])),
        rna_data_r2 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R2_' in x and 'fastq.gz' in x])),
        mem_g = "30G",
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 780
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate train --input {input.masked_assembly} --out {params.out_dir} \
        --left {params.rna_data_r1} --right {params.rna_data_r2} --stranded RF \
        --jaccard_clip --species "Candida auris" --isolate {wildcards.sample} --cpus {threads} --memory {params.mem_g}
        """

# This should automatically detect the four training files generated previously, even without explicit input
# All steps should run with 'pasa' or 'rna-bam' under Training-Methods. Nothing should run with 'busco'.
rule funannotate_predict:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa",
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
    output:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        genemark_path = config["funqcd_lib"] + "genemark/gmes_linux_64_4/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 500
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate predict --input {input.masked_assembly} --out {params.out_dir} \
        --species {wildcards.sample} --force \
        --busco_seed_species candida_albicans --busco_db saccharomycetes_odb10 --cpus {threads} \
        --GENEMARK_PATH {params.genemark_path}
        """

rule funannotate_update:
    input:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 600
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate update --input {params.out_dir} --cpus {threads}
        """

# this step removes several large files generated by funannotate_train
# these files are needed for funannotate_update, but are too large for permanent storage
# funannotate_update can still be run with the remaining files present in the normalize/ directory
rule cleanup:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
    output:
        cleanup_check = "results/{prefix}/funannotate/{sample}/training/{sample}_cleanup_complete.txt",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
    threads: 1
    resources:
        mem_mb = 5000,
        runtime = 60,
    shell:
        """
        rm -f {params.out_dir}training/left.fq.gz
        rm -f {params.out_dir}training/right.fq.gz
        rm -r -f {params.out_dir}training/trimmomatic/
        rm -r -f {params.out_dir}training/trinity_gg/
        rm -f {params.out_dir}training/normalize/*CV10000.fq
        touch {output.cleanup_check}
        """

rule interproscan:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
    output:
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml"
    singularity:
        "docker://interpro/interproscan:5.71-102.0"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/interproscan/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 360
    shell:
        """
        bash /opt/interproscan/interproscan.sh --input {input.funannotate_update_proteins} --output-dir {params.out_dir} \
        --disable-precalc --cpu {threads}
        """


rule eggnog:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
    output:
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations"
    singularity:
        "docker://nanozoo/eggnog-mapper:2.1.9--4f2b6c0"
    params:
        eggnog_data_dir = config["funqcd_lib"] + "eggnog_data/",
        out_dir = "results/{prefix}/funannotate/{sample}/eggnog/",
    threads: 8
    resources:
        mem_mb = 10000,
        runtime = 300,
    shell:
        """
        emapper.py -i {input.funannotate_update_proteins} --itype proteins --data_dir {params.eggnog_data_dir} -m diamond \
        --output {wildcards.sample} --output_dir {params.out_dir} --cpu {threads} --override
        """

rule funannotate_annotate:
    input:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml",
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations",
        cleanup_check = "results/{prefix}/funannotate/{sample}/training/{sample}_cleanup_complete.txt",
    output:
        funannotate_annotate_proteins = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa",
        funannotate_annotate_assembly = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.scaffolds.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        funannotate_update_dir = "results/{prefix}/funannotate/{sample}/update_results/",
    threads: 8
    resources:
        mem_mb = 3000,
        runtime = 80
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate annotate -i {params.funannotate_update_dir} -o {params.out_dir} --cpus {threads} \
        --iprscan {input.interproscan_out} --eggnog {input.eggnog_out} --busco_db saccharomycetes_odb10
        """

# The line 'rm -rf RM_*' removes the directories that RepeatMasker generates in the working directory
rule busco_final:
    input:
        funannotate_annotate_proteins = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa", prefix = PREFIX, sample = SAMPLE),
        funannotate_annotate_nucleotides = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.scaffolds.fa", prefix = PREFIX, sample = SAMPLE),       
    output:
        busco_out_p = "results/{prefix}/busco/busco_output_prot/batch_summary.txt",
        busco_out_n = "results/{prefix}/busco/busco_output_nucl/batch_summary.txt",
    params:
        busco_db = config["funqcd_lib"] + "busco/"
    threads: 8
    resources:
        mem_mb = 20000,
        runtime = 600
    singularity:
        "docker://ezlabgva/busco:v5.7.0_cv1"
    shell:
        """
        mkdir -p results/{wildcards.prefix}/busco/input/prot/
        mkdir -p results/{wildcards.prefix}/busco/input/nucl/
        cp results/{wildcards.prefix}/funannotate/*/annotate_results/*.proteins.fa results/{wildcards.prefix}/busco/input/prot
        cp results/{wildcards.prefix}/funannotate/*/annotate_results/*.scaffolds.fa results/{wildcards.prefix}/busco/input/nucl
        busco -f --in results/{wildcards.prefix}/busco/input/prot --mode protein --lineage_dataset saccharomycetes_odb10 --out_path results/{wildcards.prefix}/busco/ -c {threads} --out busco_output_prot --offline --download_path {params.busco_db}
        busco -f --in results/{wildcards.prefix}/busco/input/nucl --mode genome --lineage_dataset saccharomycetes_odb10 --out_path results/{wildcards.prefix}/busco/ -c {threads} --out busco_output_nucl --offline --download_path {params.busco_db}
        rm -rf RM_*
        """

rule multiqc:
    input:
        quast_out_flye = expand("results/{prefix}/quast/{sample}/{sample}_flye/report.txt",prefix=PREFIX, sample=SAMPLE),
        quast_out_medaka = expand("results/{prefix}/quast/{sample}/{sample}_medaka/report.txt",prefix=PREFIX, sample=SAMPLE),
        nanoplot_preqc = expand("results/{prefix}/nanoplot/{sample}/{sample}_preqcNanoPlot-report.html",prefix=PREFIX, sample=SAMPLE),
        busco_out_p = "results/{prefix}/busco/busco_output_prot/batch_summary.txt",
        busco_out_n = "results/{prefix}/busco/busco_output_nucl/batch_summary.txt",
    output:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
    params:
        outdir = "results/{prefix}/multiqc",
        prefix = "{prefix}",
        quast_dir = "results/{prefix}/quast/",
        nanoplot_dir = "results/{prefix}/nanoplot/",
        busco_n_dir = "results/{prefix}/busco/busco_output_nucl/",
        busco_p_dir = "results/{prefix}/busco/busco_output_prot/",
    resources:
        mem_mb = 1000,
        runtime = 120
    threads: 1
    singularity:
        "docker://multiqc/multiqc:v1.25.1"
    shell:
        """
        multiqc -f --outdir {params.outdir} -n {params.prefix}_QC_report -i {params.prefix}_QC_report \
        {params.quast_dir} {params.busco_n_dir} {params.busco_p_dir} {params.nanoplot_dir}
        """

rule qc_report_final:
    input:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
        auriclass_report = expand("results/{prefix}/auriclass/{sample}/{sample}_report.tsv", sample = SAMPLE, prefix = PREFIX), 
        nanostat_report = expand("results/{prefix}/nanoplot/{sample}/{sample}_postqcNanoStats.txt", sample = SAMPLE, prefix = PREFIX),
        flye_assembly = expand("results/{prefix}/flye/{sample}/{sample}_flye.fasta", sample = SAMPLE, prefix = PREFIX),
    output:
        summary_output = "results/{prefix}/multiqc/{prefix}_final_qc_summary.tsv"
    params:
        multiqc_dir = "results/{prefix}/multiqc/{prefix}_QC_report_data/",
        auriclass_dir = "results/{prefix}/auriclass/", 
        nanostat_dir = "results/{prefix}/nanoplot/",
        flye_dir = "results/{prefix}/flye/",
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 30,
    script:
        "QC_summary.py"
    

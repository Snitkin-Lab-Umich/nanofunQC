
import pandas as pd
import os
import json
from functools import reduce
import numpy as np

def make_summary_report(input_path, output_path, report, type = 'tsv'):
    input_path,output_path = [x+'/' if x[-1] != '/' else x for x in [input_path,output_path]]
    # this is the path to the folder with each individual auriclass report
    input_dict = {}
    for input_dir in [x for x in os.listdir(input_path) if os.path.isdir(input_path + x)]:
        for input_file in os.listdir(input_path + input_dir):
            if '.' + type in input_file:
                if type == 'tsv':
                    fdata = pd.read_csv(input_path + input_dir + '/' + input_file, sep = '\t', header = 0)
                    if fdata['Sample'][0] in input_dict:
                        print('Error: overlapping data in input directories')
                    input_dict[fdata['Sample'][0]] = [fdata['QC_decision'][0], fdata['Clade'][0]]
                if type == 'json':
                    with open(input_path + input_dir + '/' + input_file, 'r') as fh:
                        fdata = json.load(fh)
                    sample = input_file.split('_coverage.json')[0]
                    if sample in input_dict:
                        print('Error: overlapping data in input directories')
                    fdata = [fdata['qc_stats']['read_total'], fdata['qc_stats']['total_bp'], fdata['qc_stats']['read_mean'], fdata['qc_stats']['coverage']]
                    input_dict[sample] = [str(x) for x in fdata]
                if type == 'txt' and 'postqcNanoStats.txt' in input_file:
                    sample = input_dir
                    with open(input_path + input_dir + '/' + input_file, 'r') as fh:
                        for line in fh:
                            line = line.strip().split('\t')
                            if line[0] == 'n50':
                                sample_read_n50 = line[1]
                            if line[0] == 'mean_qual':
                                sample_mean_read_qual = line[1]
                            if line[0] == 'mean_read_length':
                                sample_mean_read_length = line[1]
                        input_dict[sample] = [sample_read_n50,sample_mean_read_qual,sample_mean_read_length]
    with open(output_path + report + '_summary.tsv','w') as fhout:
        if report == 'auriclass':
            _ = fhout.write('Sample\tauriclass_QC_decision\tauriclass_clade\n')
        if report == 'raw_coverage':
            _ = fhout.write('Sample\ttotal_reads\ttotal_bp\tmean_read_length\tcoverage\n')
        if report == 'nanostat_postqc':
            _ = fhout.write('Sample\tsample_read_n50\tsample_mean_read_qual\tsample_mean_read_length\n')        
        for s in input_dict:
            _ = fhout.write('\t'.join([s] + input_dict[s]) + '\n')


def make_flye_coverage_report(input_path,output_path):
    # input path and output path should both be the flye directory in results/
    input_path,output_path = [x+'/' if x[-1] != '/' else x for x in [input_path,output_path]]
    # this is the path to the folder with each individual auriclass report
    input_dict = {}
    # input_dir is the sample name, such as 'UM_Caur_nano_1'
    for input_dir in [x for x in os.listdir(input_path) if os.path.isdir(input_path + x)]:
        for input_file in os.listdir(input_path + input_dir):
            if input_file == 'assembly_info.txt':
                fdata = pd.read_csv(input_path + input_dir + '/' + input_file, sep = '\t', header = 0)
                percent_of_genome = fdata['length']/sum(fdata['length'])
                average_coverage = sum(fdata['cov.'] * percent_of_genome)
                input_dict[input_dir] = [str(round(average_coverage,2))]
    with open(output_path + 'flye_coverage_summary.tsv','w') as fhout:
        _ = fhout.write('Sample\taverage_coverage\n')
        for s in input_dict:
            _ = fhout.write('\t'.join([s] + input_dict[s]) + '\n')
    

def final_qc_summary(multiqc_path,auriclass_path,nanostat_path,flye_coverage_path):

    if multiqc_path[-1] != '/':
        multiqc_path+='/'
    #fastqc_path = multiqc_path + 'multiqc_fastqc.txt'
    quast_path = multiqc_path + 'multiqc_quast.txt'
    busco_path = multiqc_path + 'multiqc_busco.txt'

    # load quast report
    data_quast = pd.read_csv(quast_path, sep = '\t', header = 0)
    # data_quast['Sample'].replace('_contigs_l1000','',regex=True, inplace=True)
    data_quast = data_quast.loc[['_medaka' in x for x in data_quast['Sample']],:]
    # this removes the rows with _flye
    data_quast['Sample'] = data_quast['Sample'].replace('_medaka','',regex=True)
    data_quast = data_quast[['Sample','N50','Total length','# contigs (>= 0 bp)']]

    # load busco report
    data_busco = pd.read_csv(busco_path, sep = '\t', header = 0)
    data_busco = data_busco[~(data_busco['Sample'] == 'short')]
    # removes the 'short' row
    data_busco['Sample'] = [x[-1] for x in data_busco['Sample'].str.split('odb10.')]
    # sample names should now be either {sample} or {sample}.proteins
    data_busco['Percent Complete BUSCOs'] = round((data_busco['complete'] / data_busco['total']) * 100,2)
    data_busco = data_busco[['Sample','Percent Complete BUSCOs']]
    # calculates and rounds the percentage of complete buscos, removes other columns
    data_busco_prot = data_busco.copy()[data_busco['Sample'].str.contains('.proteins')]
    data_busco_nucl = data_busco.copy()[~(data_busco['Sample'].str.contains('.proteins'))]
    # splits busco into prot and nucl results
    # data_busco_prot['Sample'].replace('.proteins','',regex=True, inplace=True)
    data_busco_prot['Sample'] = data_busco_prot['Sample'].replace('.proteins','',regex=True)
    data_busco_prot.rename(columns = {'Percent Complete BUSCOs':'Percent Complete Protein BUSCOs'}, inplace = True)
    data_busco_nucl.rename(columns = {'Percent Complete BUSCOs':'Percent Complete Nucleotide BUSCOs'}, inplace = True)
    # renames columns to distinguish nucl from prot
    data_busco_merge = pd.merge(left=data_busco_nucl,right=data_busco_prot,how='inner',on='Sample')
    # combine tables

    # load auriclass report
    data_auriclass = pd.read_csv(auriclass_path, sep = '\t', header = 0)

    # load nanostat report
    data_nanostat = pd.read_csv(nanostat_path, sep = '\t', header = 0)

    # load coverage report
    data_flye_coverage = pd.read_csv(flye_coverage_path, sep = '\t', header = 0)

    # merge everything on the sample column
    data_list = [data_nanostat, data_quast, data_busco_merge, data_auriclass,data_flye_coverage]
    data_merge = reduce(lambda left,right: pd.merge(left,right,on='Sample',how='inner'), data_list)
    
    # choose which tests to use for fastqc 
    # fastqc_tests = [
    #     'after_trim_basic_statistics','after_trim_per_base_sequence_quality','after_trim_per_tile_sequence_quality',
    #     'after_trim_per_sequence_quality_scores','after_trim_per_base_sequence_content','after_trim_per_sequence_gc_content',
    #     'after_trim_per_base_n_content','after_trim_sequence_length_distribution','after_trim_sequence_duplication_levels',
    #     'after_trim_overrepresented_sequences','after_trim_adapter_content'
    #     ]
    # data_merge['fastqc_tests_passed'] = [sum(data_merge.loc[x,fastqc_tests] == 'pass') for x in range(len(data_merge.index))]
    
    return(data_merge)


def qc_evaluate(input_df, n50_min, contig_number_max, contig_number_min, assembly_length_max, assembly_length_min, mean_read_length_min, mean_read_quality_min, busco_n_score_min, flye_coverage_min):
    qc_check = [
        input_df['N50'] < n50_min,
        input_df['# contigs (>= 0 bp)'] > contig_number_max,
        input_df['# contigs (>= 0 bp)'] < contig_number_min,
        input_df['Total length'] > assembly_length_max,
        input_df['Total length'] < assembly_length_min,
        input_df['sample_mean_read_length'] < mean_read_length_min,
        input_df['sample_mean_read_qual'] < mean_read_quality_min,
        #input_df['coverage'] < average_coverage_min,
        #input_df['fastqc_tests_passed'] < fastqc_tests_passed_min,
        input_df['Percent Complete Nucleotide BUSCOs'] < busco_n_score_min,
        input_df['average_coverage'] < flye_coverage_min,
    ]
    
    status = ['FAIL'] * len(qc_check)
    
    input_df['QC_EVALUTATION'] = np.select(qc_check, status, default='PASS')
    return(input_df)


def main(multiqc_path, auriclass_path, nanostat_path, flye_path, output_path, n50_min, contig_number_max, contig_number_min, assembly_length_max, assembly_length_min, mean_read_length_min, mean_read_quality_min, busco_n_score_min, flye_coverage_min):
    auriclass_path,nanostat_path = [x+'/' if x[-1] != '/' else x for x in [auriclass_path,nanostat_path]]
    # generate auriclass summary
    make_summary_report(input_path = auriclass_path, output_path = auriclass_path, report = 'auriclass', type = 'tsv')
    # generate nanostat summary
    make_summary_report(input_path = nanostat_path, output_path = nanostat_path, report = 'nanostat_postqc', type = 'txt')
    # generate flye coverage report
    make_flye_coverage_report(input_path = flye_path,output_path = flye_path)
    mult = multiqc_path
    auri = auriclass_path + 'auriclass_summary.tsv'
    nano = nanostat_path + 'nanostat_postqc_summary.tsv'
    fcov = flye_path + 'flye_coverage_summary.tsv'
    final = final_qc_summary(multiqc_path = mult, auriclass_path = auri, nanostat_path = nano, flye_coverage_path = fcov)
    final_eval = qc_evaluate(
        final,n50_min=n50_min,contig_number_max=contig_number_max,contig_number_min=contig_number_min,assembly_length_max=assembly_length_max,assembly_length_min=assembly_length_min,
        mean_read_length_min=mean_read_length_min, mean_read_quality_min=mean_read_quality_min,busco_n_score_min=busco_n_score_min,flye_coverage_min=flye_coverage_min
        )
    final_eval.to_csv(output_path, sep='\t',index=False)

if __name__ == "__main__":
    #inp = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/auriclass/'
    #outp = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/auriclass/'
    #inp = '/Users/josephhale/Desktop/qc/auriclass'
    #outp = '/Users/josephhale/Desktop/qc/auriclass'
    #make_summary_report(input_path = inp, output_path = outp, report = 'auriclass', type = 'tsv')
    #inp = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/raw_coverage/'
    #outp = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/raw_coverage/'
    #inp = '/Users/josephhale/Desktop/qc/raw_coverage'
    #outp = '/Users/josephhale/Desktop/qc/raw_coverage'
    #make_summary_report(input_path = inp, output_path = outp, report = 'raw_coverage', type = 'json')
    #mult = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/multiqc/2024-11-14_Chicago_strains_v2_QC_report_data/'
    #aur = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/auriclass/summary.tsv'
    #rcov = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/raw_coverage/summary.tsv'
    #outp = '/nfs/turbo/umms-esnitkin/Project_Cauris/Analysis/2024_Pipeline_testing/2024_12_16_nanofunsake/funQCD/results/2024-11-14_Chicago_strains_v2/multiqc/final_qc_summary.tsv'
    #mult = '/Users/josephhale/Desktop/qc/multiqc/2024-11-14_Chicago_strains_v2_QC_report_data/'
    #aur = '/Users/josephhale/Desktop/qc/auriclass/summary.tsv'
    #rcov = '/Users/josephhale/Desktop/qc/raw_coverage/summary.tsv'
    #outp = '/Users/josephhale/Desktop/qc/multiqc/final_qc_summary.tsv'
    #final_qc_summary(multiqc_path = mult,auriclass_path = aur,raw_coverage_path = rcov, output_path = outp)
    #mult = '/Users/josephhale/Desktop/qc/multiqc/2024-11-14_Chicago_strains_v2_QC_report_data/'
    #auri = '/Users/josephhale/Desktop/qc/auriclass/'
    #rcov = '/Users/josephhale/Desktop/qc/raw_coverage/'
    #outp = '/Users/josephhale/Desktop/qc/multiqc/final_qc_summary.tsv'
    #n50_min = 50000
    #contig_number_max = 500
    #contig_number_min = 10
    #assembly_length_max = 13000000
    #average_coverage_min = 20
    #fastqc_tests_passed_min = 7
    #busco_n_score_min = 97
    #main(mult, auri, rcov, outp, n50_min, contig_number_max, contig_number_min, assembly_length_max, average_coverage_min, fastqc_tests_passed_min, busco_n_score_min)
    mult = snakemake.params["multiqc_dir"]
    auri = snakemake.params["auriclass_dir"]
    nano = snakemake.params["nanostat_dir"]
    flye = snakemake.params["flye_dir"]
    outp = snakemake.output[0]
    n50_min = snakemake.config["min_n50"]
    contig_number_min = snakemake.config["min_contigs"]
    contig_number_max = snakemake.config["max_contigs"]
    assembly_length_min = snakemake.config["min_assembly_length"]
    assembly_length_max = snakemake.config["max_assembly_length"]
    #average_coverage_min = snakemake.config["min_avg_coverage"]
    #fastqc_tests_passed_min = snakemake.config["min_fastqc_tests_passed"]
    mean_read_length_min = snakemake.config["min_read_length"]
    mean_read_quality_min = snakemake.config["min_read_quality"]
    busco_n_score_min = snakemake.config["min_busco_nucl_score"]
    flye_coverage_min = snakemake.config["min_average_coverage"]
    main(
        multiqc_path=mult, auriclass_path=auri, nanostat_path=nano, flye_path=flye, output_path=outp, n50_min=n50_min, contig_number_max=contig_number_max, contig_number_min=contig_number_min, 
        assembly_length_max=assembly_length_max, assembly_length_min=assembly_length_min, mean_read_length_min=mean_read_length_min, mean_read_quality_min=mean_read_quality_min, 
        busco_n_score_min=busco_n_score_min, flye_coverage_min=flye_coverage_min
        )
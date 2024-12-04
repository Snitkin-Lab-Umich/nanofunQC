# Author: Dhatri Badri 
#configfile: "config/config.yaml"

import pandas as pd
import os
import numpy as np

samples_df = pd.read_csv(config["samples"])
BARCODE = list(samples_df['barcode_id'])
PREFIX = config["prefix"]

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

# Concatenate nanoplot outputs
def nanoplot_report(outdir, prefix):
    prefix = prefix.pop()
    outdir = "results/%s" % prefix
    report_dir = str(outdir) + "/%s_report" % prefix
    nanoplot_dir = os.path.join(outdir, "nanoplot")

    # Define the columns to extract and the header for the CSV
    metrics_to_extract = ["n50", "median_read_length", "mean_read_length", "number_of_reads"]
    header = ["sample_long_read"] + metrics_to_extract
    
    # Define the column renaming map
    column_rename = {
        'number_of_reads': 'number_of_reads_nanoplot',
        'n50': 'N50_nanoplot',
        'median_read_length': 'median_read_length_nanoplot',
        'mean_read_length': 'mean_read_length_nanoplot'
    }

    # List to hold all samples' data
    all_samples_data = []

    # Loop through each sample directory in the input directory
    for sample in os.listdir(nanoplot_dir):
        sample_path = os.path.join(nanoplot_dir, sample)

        # Define the path to the sample_postqcNanoStats.txt file
        stats_file = os.path.join(sample_path, f"{sample}_postqcNanoStats.txt")
        
        if os.path.exists(stats_file):
            # Initialize a dictionary to hold the metrics
            metrics = {}

            # Read the stats file and extract the required metrics
            with open(stats_file, "r") as file:
                for line in file:
                    line = line.strip()
                    if not line:  # Skip empty lines
                        continue
                    
                    # Split the line into key and value
                    key, value = line.split("\t")
                    
                    # Check if the key is in the list of metrics to extract
                    if key in metrics_to_extract:
                        metrics[key] = value
            
            # Add the current sample's data to the list of all samples' data
            sample_data = [sample] + [metrics.get(metric, '') for metric in metrics_to_extract]
            all_samples_data.append(sample_data)

    # Convert the list of all samples' data into a pandas DataFrame
    result_df = pd.DataFrame(all_samples_data, columns=header)

    # Rename columns in the DataFrame
    result_df.rename(columns=column_rename, inplace=True)

    # Create the output directory if it doesn't exist
    os.makedirs(report_dir, exist_ok=True)
    
    # Define the output CSV file path
    output_file = os.path.join(report_dir, "%s_nanoplot_results.csv" % prefix)

    # Write the DataFrame to a CSV file
    result_df.to_csv(output_file, index=False)

def quast_report(outdir, prefix):
    prefix = prefix.pop()
    outdir = "results/%s" % prefix
    report_dir = str(outdir) + "/%s_report" % prefix
    quast_dir = os.path.join(outdir, "quast")
    
    # Define the columns to extract and the header for the CSV
    metrics_to_extract = ["# contigs", "N50", "Total length"]
    header = ["sample_long_read"] + metrics_to_extract
    
    # Define the column renaming map
    column_rename = {
        '# contigs': 'number_of_contigs_medaka_quast',
        'N50': 'N50_medaka_quast',
        'Total length': 'total_length_medaka_quast'
    }

    # List to hold all samples' data
    all_samples_data = []

    # Loop through each sample directory in the input directory
    for sample in os.listdir(quast_dir):
        sample_path = os.path.join(quast_dir, sample)

        # Construct the path to the QUAST output directory (e.g., Lojek_G6C_1_medaka)
        quast_output_dir = os.path.join(sample_path, f"{sample}_medaka")
        
        # Define the path to the transposed.tsv file
        transposed_file = os.path.join(quast_output_dir, "transposed_report.tsv")
        
        # Debug: Check if the file path is correct
        #print(f"Checking for file: {transposed_file}")
        
        if os.path.exists(transposed_file):
            #print(f"File found: {transposed_file}")
            
            # Read the transposed.tsv file into a pandas DataFrame
            try:
                df = pd.read_csv(transposed_file, sep="\t")
                
                # Debug: Check if the DataFrame is empty
                if df.empty:
                    print(f"File {transposed_file} is empty.")
                else:
                    # Print the first row of the DataFrame for debugging
                    #print(df.iloc[0])
                    
                    # Initialize a list to hold the metrics for the current sample
                    sample_data = [sample]

                    # Extract the required metrics from the only row
                    row = df.iloc[0]
                    for metric in metrics_to_extract:
                        sample_data.append(row[metric])

                    # Add the current sample's data to the list of all samples' data
                    all_samples_data.append(sample_data)
            except Exception as e:
                print(f"Error reading {transposed_file}: {e}")
        else:
            print(f"File not found: {transposed_file}")

    # Create the output directory if it doesn't exist
    os.makedirs(report_dir, exist_ok=True)
    
    # Define the output CSV file path
    output_file = os.path.join(report_dir, "%s_quast_results.csv" % prefix)

    # Convert the list of all samples' data into a pandas DataFrame
    result_df = pd.DataFrame(all_samples_data, columns=header)

    # Rename columns in the DataFrame
    result_df.rename(columns=column_rename, inplace=True)
    
    # Write the DataFrame to a CSV file
    result_df.to_csv(output_file, index=False)

# Create skani report
def skani_report(outdir, prefix):
    prefix = prefix.pop()
    outdir = "results/%s" % prefix
    report_dir = str(outdir) + "/%s_report" % prefix
    skani_dir = os.path.join(outdir, "skani")
    #report_data_dir = report_dir + "/data"

    result_df = pd.DataFrame(columns=['sample_long_read', 'ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name']) # Create an empty DataFrame to store the results

    for root, dirs, files in os.walk(skani_dir):  # Iterate over the directory structure in the 'results/prefix' directory
        for sample_name in dirs:  # Iterate over each sample directory
            skani_file_path = os.path.join(root, sample_name, f'{sample_name}_skani_output.txt')

            if os.path.exists(skani_file_path):  # Check if the skani file exists
                skani_file = pd.read_csv(skani_file_path, sep='\t| ,', skipinitialspace=True, header=0, engine='python')  # Read the skani file
                first_row_df = skani_file[['ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name']].iloc[:1]  # Extract the first row from specific columns

                if first_row_df.empty:  # Check if the first row is empty
                    first_row_df = pd.DataFrame({'sample_long_read': [sample_name],  # Create a DataFrame with NaN values and add the sample name
                                                'ANI': ["NAs"],
                                                'Align_fraction_ref': ["NAs"],
                                                'Align_fraction_query': ["NAs"],
                                                'Ref_name': ["NAs"]})
                else:
                    first_row_df.loc[:, 'sample_long_read'] = sample_name  # Add a new column for the sample name

                first_row_df = first_row_df[['sample_long_read', 'ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name']]  # Reorder the columns
                result_df = pd.concat([result_df, first_row_df], ignore_index=True)  # Concatenate the result to the overall DataFrame
    
    result_file_path = os.path.join(report_dir, '%s_Skani_report_final.csv' % prefix) # Save the final result to CSV files in the output directory
    result_df.to_csv(result_file_path, index=False)

def summary(outdir, prefix):
    prefix = prefix.pop()
    outdir = outdir.pop()
    # Define the paths to the various files
    # Organize reports directory
    report_dir = os.path.join(outdir, "%s_report" % prefix)
    #report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix
    nanoplot_file = os.path.join(report_dir, "%s_nanoplot_results.csv" % prefix)
    quast_file = os.path.join(report_dir, "%s_quast_results.csv" % prefix)
    mlst_file = os.path.join(report_dir, "%s_mlst_results.csv" % prefix)
    skani_file = os.path.join(report_dir, "%s_Skani_report_final.csv" % prefix)
    flye_dir= os.path.join(outdir,"flye")
    
    # Load the nanoplot and quast results into pandas DataFrames
    nanoplot_df = pd.read_csv(nanoplot_file)
    quast_df = pd.read_csv(quast_file)

    # Merge the DataFrames on the 'sample_long_read' column
    merged_df = pd.merge(nanoplot_df, quast_df, on='sample_long_read')

    # Initialize the num_of_circ_contigs column with 0
    merged_df['num_of_circ_contigs_flye'] = 0

    # Loop through each sample and count the circular contigs
    for sample in merged_df['sample_long_read']:
        assembly_info_path = os.path.join(flye_dir, sample, "assembly_info.txt")

        if os.path.exists(assembly_info_path):
            # Count the number of 'Y' entries in the circ. column
            with open(assembly_info_path, "r") as file:
                circular_contigs_count = sum(1 for line in file if line.strip().split("\t")[3] == "Y")

            # Update the DataFrame with the count
            merged_df.loc[merged_df['sample_long_read'] == sample, 'num_of_circ_contigs_flye'] = circular_contigs_count

    # Load the MLST results file into a DataFrame
    mlst = pd.read_csv(mlst_file, sep='\t', header=0)
    #print(mlst.columns)  # Print the column names to debug

    # Apply the lambda function to the 'sample_long_read' column to extract the sample name
    mlst['sample_long_read'] = mlst['sample_long_read'].apply(lambda x: x.split('/')[-1].split('_medaka.fasta')[0])

    # Read final Skani output file
    skani_summary = pd.read_csv(skani_file, sep=',', skipinitialspace=True, header=0, engine='python')

    # Merge all dataframes
    merged_df_temp = pd.merge(merged_df, mlst, on='sample_long_read')
    merged_df_temp_2 = pd.merge(merged_df_temp, skani_summary, on='sample_long_read')

    # Calculate the coverage
    merged_df_temp_2['coverage'] = (merged_df_temp_2['number_of_reads_nanoplot'] * merged_df_temp_2['mean_read_length_nanoplot']) / merged_df_temp_2['total_length_medaka_quast']

    # Determine the QC check status
    merged_df_temp_2['qc_check'] = (merged_df_temp_2['coverage'] > 30) & (merged_df_temp_2['number_of_contigs_medaka_quast'] < 20)
    merged_df_temp_2['qc_check'] = merged_df_temp_2['qc_check'].apply(lambda x: 'PASS' if x else 'FAIL')

    # Define the output file path
    output_file = os.path.join(report_dir, "%s_report.csv" % prefix)

    # Save the final DataFrame to a CSV file
    merged_df_temp_2.to_csv(output_file, index=False)

rule concatenate_tool_outputs:
    input:
        mlst_out = expand("results/{prefix}/mlst/{barcode}/report.tsv", prefix=PREFIX, barcode=BARCODE),
        skani_out = expand("results/{prefix}/skani/{barcode}/{barcode}_skani_output.txt", prefix=PREFIX, barcode=BARCODE),
        quast_flye_out= expand("results/{prefix}/quast/{barcode}/{barcode}_flye/report.txt", prefix=PREFIX, barcode=BARCODE),
        quast_medaka_out = expand("results/{prefix}/quast/{barcode}/{barcode}_medaka/report.txt", prefix=PREFIX, barcode=BARCODE),
        nanoplot_out = expand("results/{prefix}/nanoplot/{barcode}/{barcode}_preqcNanoPlot-report.html", prefix=PREFIX, barcode=BARCODE),
    output:
        nanoplot_report = f"results/{{prefix}}/{{prefix}}_report/{{prefix}}_nanoplot_results.csv",
        quast_report = f"results/{{prefix}}/{{prefix}}_report/{{prefix}}_quast_results.csv",
        skani_report =f"results/{{prefix}}/{{prefix}}_report/{{prefix}}_Skani_report_final.csv",
        mlst_report = f"results/{{prefix}}/{{prefix}}_report/{{prefix}}_mlst_results.csv",
    params:
        prefix = "{prefix}",
        outdir = "results/{prefix}",
        mlst_dir = "results/{prefix}/mlst"
    run:
        nanoplot_report({params.outdir},{params.prefix})
        quast_report({params.outdir},{params.prefix})
        skani_report({params.outdir},{params.prefix})
        shell("""
        echo "sample_long_read\tScheme\tST" > {output.mlst_report} && cut -f1-3 {input.mlst_out} >> {output.mlst_report}
        """)

rule summary: 
    input:
        nanoplot_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_report/{wildcards.prefix}_nanoplot_results.csv"),
        quast_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_report/{wildcards.prefix}_quast_results.csv"),
        skani_report =lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_report/{wildcards.prefix}_Skani_report_final.csv"),
        mlst_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_report/{wildcards.prefix}_mlst_results.csv"),
    output:
        summary_report = f"results/{{prefix}}/{{prefix}}_report/{{prefix}}_report.csv"
    params:
        prefix = "{prefix}",
        outdir = "results/{prefix}",
    run:
        summary({params.outdir},{params.prefix})

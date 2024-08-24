# Overview
The script, biolog_proc.py, processes raw BIOLOG data files to identify the metabolites that support the growth of specfic strains. I implemented three quantitative metrics for this assessment: (1) endpoint OD; (2) area under the growth curve; and (3) specific growth rate. Each metric is evaluated based on two criteria: the average fold change of the metric (compared to negative control) must be >= `fc_cutoff`, and the p-value (based on a paried-sample t-test) must be < `pvalue_cutoff`.

Below is a guide on how to run the pipeline with various command-line arguments.

# Basic Usage
Download the folder and put your raw BIOLOG files into the folder `input_data_folder`. To run the script, use the following command:
python3 your_script_name.py input_folder_path [optional arguments]

## Positional Argument
input_folder_path (required): Specifies the path to the folder containing the input Biolog output files.

## Optional Arguments
--growth_model (optional): Specifies the growth curve model to be used for data fitting. The default model is 'LGX'.
To use the Gompertz model, run the following command:
python3 your_script_name.py /path/to/biolog/output/files --growth_model Gompertz

--min_r2 (optional): Specifies the minimum R² (coefficient of determination) required for the growth curve model fitting. The default value is 0.9.
To use a more stringent R² cutoff, run the following command:
python3 your_script_name.py /path/to/biolog/output/files --min_r2 0.95

--max_trials (optional): Sets the maximum number of trial attempts for initial guesses during the growth curve model fitting process. The default value is 50.
To use more trial attempts, run the following command:
python3 your_script_name.py /path/to/biolog/output/files --max_trials 100

--fc_cutoff (optional): Specifies the minimum mean fold change cutoff for considering growth. The default value is 1.0.
To use a more stringent cutoff for fold change, run the following command:
python3 your_script_name.py /path/to/biolog/output/files --fc_cutoff 1.5

--pvalue_cutoff (optional): Specifies the maximum p-value cutoff for growth determination. The default value is 0.05.
To use a more stringent pvalue cutoff, run the following command:
python3 your_script_name.py /path/to/biolog/output/files --pvalue_cutoff 0.01

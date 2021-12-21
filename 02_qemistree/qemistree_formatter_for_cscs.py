# %%

import pandas as pd
import sys
import shlex
import zipfile
import glob
import os
import subprocess


# %% First we externalize the variables we'll use in the rest of the script

job_id = '044e981ff0d84246ae5c91ef3db643a8'
gnps_job_path = '/path/of/your/choice/'
project_name = 'Qemistree_set'
base_filename = 'GNPS_output_' + project_name
filename_suffix = 'zip'
path_to_folder = os.path.join(gnps_job_path, base_filename)
path_to_file = os.path.join(gnps_job_path, base_filename + "." + filename_suffix)



# Actually here we also want to have acces to the GNPS data
# %% Downloading GNPS files


job_url_zip = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task="+job_id+"&view=download_cytoscape_data"

cmd = 'curl -d "" '+job_url_zip+' -o '+path_to_file
subprocess.call(shlex.split(cmd))

with zipfile.ZipFile(path_to_file, 'r') as zip_ref:
    zip_ref.extractall(path_to_folder)

# We finally remove the zip file
cmd = 'rm '+ path_to_file
subprocess.call(shlex.split(cmd))


# Here we get the peaklist from the GNPS folder

# the feature table path is generated from the base bath to the GNPS results folder

feature_table_path = os.path.join(path_to_folder,'quantification_table_reformatted','')

df_pl = pd.read_csv(feature_table_path + str(os.listdir(feature_table_path)[0]), sep=',')


# %% See here for formatting requirement in qiime metdata artifacts https://docs.qiime2.org/2020.11/tutorials/metadata/

df_pl.columns = df_pl.columns.str.replace(" Peak area", "")
df_pl = df_pl.sort_index(axis=1)
df_pl.rename(columns={'row ID': 'feature-id'}, inplace=True)
df_pl.reset_index(inplace=True)
df_pl.set_index('feature-id', inplace=True)
df_pl = df_pl[df_pl.columns.drop(list(df_pl.filter(regex='Unnamed: ')))]
# Here we filter only the columns containing the .mz pattern 
df_pl = df_pl[df_pl.columns[df_pl.columns.str.contains('.mz')]]



# We now export this one and check if it can be converted to biom table 

df_pl.to_csv(feature_table_path + 'feature_table_for_biom.tsv', sep = '\t', index = True)


# We now go to the feature_table_path folder and run the following command line (note that biom should be installed)

`biom convert -i feature_table_for_biom.tsv -o biom_feature_table.biom --table-type="OTU table" --to-hdf5`

# This should output a biom file.
# This file is the converted to a qiime FeatureTable[Frequency] using the following command line


`qiime tools import --type 'FeatureTable[Frequency]' --input-path biom_feature_table.biom --output-path feature_table.qza`


# We now launch the cscsc command and explicitly specify all options
# In your qiime activated env enter `qiime cscs cscs` to get all options
# here we use nohup and time options additionnaly.
# Launch the follwoing command from the root folder of the dowloaded GNPS job

# For the CSCS with normalisation and weighting

`nohup bash -c 'time qiime cscs cscs --p-css-edges networkedges_selfloop/3466497461974198a9ab8c9463d05b53..selfloop --i-features quantification_table_reformatted/feature_table.qza --p-cosine-threshold 0.7 --p-normalization --p-weighted --p-cpus 40 --o-distance-matrix qeemistree_set_cscs_distance_matrix_norm_weighted.qza'`

# For the CSCS with normalisation and no weighting

`nohup bash -c 'time qiime cscs cscs --p-css-edges networkedges_selfloop/3466497461974198a9ab8c9463d05b53..selfloop --i-features quantification_table_reformatted/feature_table.qza --p-cosine-threshold 0.7 --p-normalization --p-no-weighted --p-cpus 40 --o-distance-matrix qeemistree_set_cscs_distance_matrix_norm_unweighted.qza'`

# For the CSCS with no normalisation and weighting

`nohup bash -c 'time qiime cscs cscs --p-css-edges networkedges_selfloop/3466497461974198a9ab8c9463d05b53..selfloop --i-features quantification_table_reformatted/feature_table.qza --p-cosine-threshold 0.7 --p-no-normalization --p-weighted --p-cpus 40 --o-distance-matrix qeemistree_set_cscs_distance_matrix_nonorm_weighted.qza'`

# For the CSCS with no normalisation and no weighting

`nohup bash -c 'time qiime cscs cscs --p-css-edges networkedges_selfloop/3466497461974198a9ab8c9463d05b53..selfloop --i-features quantification_table_reformatted/feature_table.qza --p-cosine-threshold 0.7 --p-no-normalization --p-no-weighted --p-cpus 40 --o-distance-matrix qeemistree_set_cscs_distance_matrix_nonorm_unweighted.qza'`

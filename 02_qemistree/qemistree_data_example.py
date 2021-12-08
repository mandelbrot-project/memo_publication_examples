import memo_ms as memo
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import time
import scipy as sp
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from skbio.stats.distance import DistanceMatrix 
import plotly.graph_objects as go

#######################################################################################################

def conditions(df_meta):
    if ((df_meta['Proportion_Fecal_1']>0) & (df_meta['Proportion_Fecal_2']==0)& (df_meta['Proportion_Tomato']==0) & (df_meta['Proportion_NIST_1950_SRM']==0)):
        return 'Fecal_1'
    if ((df_meta['Proportion_Fecal_1']==0) & (df_meta['Proportion_Fecal_2']>0)& (df_meta['Proportion_Tomato']==0) & (df_meta['Proportion_NIST_1950_SRM']==0)):
        return 'Fecal_2'
    if ((df_meta['Proportion_Fecal_1']==0) & (df_meta['Proportion_Fecal_2']==0)& (df_meta['Proportion_Tomato']>0) & (df_meta['Proportion_NIST_1950_SRM']==0)):
        return 'Tomato'
    if ((df_meta['Proportion_Fecal_1']==0) & (df_meta['Proportion_Fecal_2']==0)& (df_meta['Proportion_Tomato']==0) & (df_meta['Proportion_NIST_1950_SRM']>0)):
        return 'Plasma'
    if ((df_meta['Proportion_Fecal_1']>0) & (df_meta['Proportion_Fecal_2']>0)& (df_meta['Proportion_Tomato']==0) & (df_meta['Proportion_NIST_1950_SRM']==0)):
        return 'Fecal_1 + Fecal_2'
    if ((df_meta['Proportion_Fecal_1']>0) & (df_meta['Proportion_Fecal_2']==0)& (df_meta['Proportion_Tomato']>0) & (df_meta['Proportion_NIST_1950_SRM']==0)):
        return 'Fecal_1 + Tomato'
    if ((df_meta['Proportion_Fecal_1']>0) & (df_meta['Proportion_Fecal_2']==0)& (df_meta['Proportion_Tomato']==0) & (df_meta['Proportion_NIST_1950_SRM']>0)):
        return 'Fecal_1 + Plasma'
    if ((df_meta['Proportion_Fecal_1']==0) & (df_meta['Proportion_Fecal_2']>0)& (df_meta['Proportion_Tomato']>0) & (df_meta['Proportion_NIST_1950_SRM']==0)):
        return 'Fecal_2 + Tomato'
    if ((df_meta['Proportion_Fecal_1']==0) & (df_meta['Proportion_Fecal_2']>0)& (df_meta['Proportion_Tomato']==0) & (df_meta['Proportion_NIST_1950_SRM']>0)):
        return 'Fecal_2 + Plasma'
    if ((df_meta['Proportion_Fecal_1']==0) & (df_meta['Proportion_Fecal_2']==0)& (df_meta['Proportion_Tomato']>0) & (df_meta['Proportion_NIST_1950_SRM']>0)):
        return 'Tomato + Plasma'
    if ((df_meta['Proportion_Fecal_1']>0) & (df_meta['Proportion_Fecal_2']>0)& (df_meta['Proportion_Tomato']>0) & (df_meta['Proportion_NIST_1950_SRM']>0)):
        return 'Fecal_1 + Fecal_2 + Tomato + Plasma' 
    else:
        return 'What is it? :)'

#######################################################################################################


''' Generate MEMO matrix from aligned feature table and spectra '''


start = time.time()
feat_table_qe = memo.FeatureTable(path="../01_input_data/02_qemistree/qemistree_memo_cscs_comparison/quantification_table-00000.csv", software='mzmine')
spectra_qe = memo.SpectraDocuments(path="../01_input_data/02_qemistree/qemistree_memo_cscs_comparison/qemistree_specs_ms.mgf", min_relative_intensity = 0.01,
            max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)

memo_qe = memo.MemoMatrix()

memo_qe.memo_from_aligned_samples(featuretable= feat_table_qe, spectradocuments= spectra_qe)
end = time.time()
print(f"computing memo from aligned took {end-start} seconds")

    # We use the filter function to remove blanks and 
memo_qe = memo_qe.filter(samples_pattern='blank|Qcmix')
feat_table_qe = feat_table_qe.filter(samples_pattern='blank|Qcmix')


''' Generate MEMO matrix from unaligned .mgf '''


    # Path to folder containning the mgf (QEC18)
path_to_file = "../01_input_data/02_qemistree/mgf_individual_files/Qemistree_QE_C18/"
    # Generating memo matrix
start = time.time()
memo_unaligned = memo.MemoMatrix()
memo_unaligned.memo_from_unaligned_samples(path_to_file, min_relative_intensity = 0.01,
            max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)
end = time.time()
print(f"computing memo from unaligned for C18 samples took {end-start} seconds")

memo_unaligned = memo_unaligned.filter(samples_pattern='blank|Qcmix')

    # Path to folder containning the mgf (QEC18RT)
path_to_file = "../01_input_data/02_qemistree/mgf_individual_files/Qemistree_QE_C18RT/"
    # Generating memo matrix
start = time.time()
memo_unaligned2 = memo.MemoMatrix()
memo_unaligned2.memo_from_unaligned_samples(path_to_file, min_relative_intensity = 0.01,
            max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)
end = time.time()
print(f"computing memo from unaligned for C18RT samples took {end-start} seconds")

memo_unaligned2= memo_unaligned2.filter(samples_pattern='blank|Qcmix')

memo_merged = memo_unaligned.merge_memo(memo_unaligned2, drop_not_in_common=False)
memo_merged.memo_matrix.index += '.mzML'
memo_merged.memo_matrix = memo_merged.memo_matrix.div(memo_merged.memo_matrix.sum(axis=1), axis=0)


''' Load and clean metadata '''


df_meta = pd.read_csv("../01_input_data/02_qemistree/1901_gradient_benchmarking_dataset_v4_sample_metadata.txt", sep='\t')
df_meta['Samplename'] = df_meta['Samplename'].str[:-6]
df_meta['Samplename'] = df_meta['Samplename'].str.replace('BLANK_', 'BLANK')
df_meta = df_meta[['Filename', 'Experiment', 'Samplename', 'Triplicate_number', 'Proportion_Fecal_1', 'Proportion_Fecal_2', 'Proportion_Tomato', 'Proportion_NIST_1950_SRM']]
df_meta['contains'] = df_meta.apply(conditions, axis=1)
df_meta['instrument'] = np.where(df_meta['Samplename'].str.contains('qTOF'), 'qTOF', 'QE')
df_meta['blank_qc'] = np.where(df_meta['Samplename'].str.contains('blank|qcmix', case = False), 'yes', 'no')
df_meta['Proportion_Fecal_1'] = df_meta['Proportion_Fecal_1'].apply(str)
df_meta['Proportion_Fecal_2'] = df_meta['Proportion_Fecal_2'].apply(str)
df_meta['Proportion_Tomato'] = df_meta['Proportion_Tomato'].apply(str)
df_meta['Proportion_NIST_1950_SRM'] = df_meta['Proportion_NIST_1950_SRM'].apply(str)


''' Generating/loading the distance matrix for pcoa '''


samples = feat_table_qe.feature_table.index.to_list()
df_meta = df_meta[df_meta['Filename'].isin(samples)]
df_meta.set_index('Filename', inplace=True)
df_meta = df_meta.loc[samples]

    # Bray-curtis
df_feature_trans = feat_table_qe.feature_table.loc[samples]
dm_braycurtis = sp.spatial.distance.pdist(df_feature_trans, 'braycurtis')

    # CSCS
dm_cscs_weighted = pd.read_csv("../01_input_data/02_qemistree/qemistree_memo_cscs_comparison/cscs_weighed_dm.tsv", sep='\t')
dm_cscs_weighted.set_index('Unnamed: 0', inplace=True)
dm_cscs_weighted = dm_cscs_weighted[samples]
dm_cscs_weighted = dm_cscs_weighted.loc[samples]

dm_cscs_unweighted = pd.read_csv("../01_input_data/02_qemistree/qemistree_memo_cscs_comparison/cscs_unweighed_dm.tsv", sep='\t')
dm_cscs_unweighted.set_index('Unnamed: 0', inplace=True)
dm_cscs_unweighted = dm_cscs_unweighted[samples]
dm_cscs_unweighted = dm_cscs_unweighted.loc[samples]

    #Qemistree
dm_qemistree_weighted = pd.read_csv("../01_input_data/02_qemistree/qemistree_memo_cscs_comparison/weighted_unifrac_dm.csv", sep=',')
dm_qemistree_weighted = dm_qemistree_weighted.drop('Unnamed: 0', axis = 1)
dm_qemistree_weighted.index = dm_qemistree_weighted.columns
dm_qemistree_weighted = dm_qemistree_weighted[samples]
dm_qemistree_weighted = dm_qemistree_weighted.loc[samples]

dm_qemistree_unweighted = pd.read_csv("../01_input_data/02_qemistree/qemistree_memo_cscs_comparison/unweighted_unifrac_dm.csv", sep=',')
dm_qemistree_unweighted = dm_qemistree_unweighted.drop('Unnamed: 0', axis = 1)
dm_qemistree_unweighted.index = dm_qemistree_unweighted.columns
dm_qemistree_unweighted = dm_qemistree_unweighted[samples]
dm_qemistree_unweighted = dm_qemistree_unweighted.loc[samples]

    #MEMO aligned
df_fingerprints = memo_qe.memo_matrix.loc[samples]
start = time.time()
dm_memo = sp.spatial.distance.pdist(df_fingerprints, 'braycurtis') 
end = time.time()
print(f"computing dm from aligned took {end-start} seconds")

    #MEMO unaligned
df_fingerprints2 = memo_merged.memo_matrix.loc[samples]
start = time.time()
dm_memo2 = sp.spatial.distance.pdist(df_fingerprints2, 'braycurtis') 
end = time.time()
print(f"computing dm from unaligned took {end-start} seconds")


''' Plotting '''


results = pd.DataFrame()
results_expvar = {}
metrics = ['Feature Table (Bray-Curtis)', 'MEMO unaligned (Bray-Curtis)', 'Qemistree weighted UniFrac', 'CSCS weighted', 'MEMO aligned (Bray-Curtis)', 'Qemistree unweighted UniFrac', 'CSCS unweighted']
pos = [[1,1], [1,2], [1,3], [1,4], [2,2], [2,3], [2,4]]
legend = [True, False, False, False, False, False, False]

dms = [dm_braycurtis,dm_memo2, dm_qemistree_weighted,dm_cscs_weighted, dm_memo, dm_qemistree_unweighted, dm_cscs_unweighted]

for metric, dm in zip(metrics, dms):
    pcoa_results = pcoa(dm)

    results[metric+'_PC1'] = pcoa_results.samples['PC1']
    results[metric+'_PC2'] = pcoa_results.samples['PC2']

    results_expvar[metric] = {}
    results_expvar[metric]['varPC1'] = round(100*pcoa_results.proportion_explained[0], 1)
    results_expvar[metric]['varPC2'] = round(100*pcoa_results.proportion_explained[1], 1)

results['contains'] = list(df_meta['contains'])
results['Experiment'] = list(df_meta['Experiment'])
results['Samplename'] = list(df_meta['Samplename'])
results['Samplename'] = results['Samplename'].str.replace('QE_C18_', '')
results['Samplename'] = results['Samplename'].str.replace('QE_C18-RTshift_', '')
results['Triplicate_number'] = list(df_meta['Triplicate_number'])
results['Proportion_Fecal_1'] = list(df_meta['Proportion_Fecal_1'].apply(str))
results['Proportion_Fecal_2'] = list(df_meta['Proportion_Fecal_2'].apply(str))
results['Proportion_Tomato'] = list(df_meta['Proportion_Tomato'].apply(str))
results['Proportion_NIST_1950_SRM'] = list(df_meta['Proportion_NIST_1950_SRM'].apply(str))
results['size'] = np.where(results['contains'].isin(['Fecal_2', 'Fecal_1', 'Tomato', 'Plasma']), 11, 7)
results['line'] = np.where(results['contains'].isin(['Fecal_2', 'Fecal_1', 'Tomato', 'Plasma']), 0.8, 0)

# Choos cat to plot
plot = 'Experiment' # Experiment contains Proportion_Fecal_1

if plot == 'contains':
    colorsIdx = {
        'Fecal_1': 'rgb(0, 95, 115)',
        'Fecal_2': 'rgb(148, 210, 189)',
        'Plasma': 'rgb(233, 216, 166)',
        'Tomato': 'rgb(155, 34, 38)',
        'Fecal_1 + Fecal_2': 'rgb(10, 147, 150)',
        'Fecal_1 + Tomato': 'rgb(202, 103, 2)',
        'Fecal_1 + Plasma': 'rgb(0, 150, 199)',
        'Fecal_1 + Fecal_2 + Tomato + Plasma': 'rgb(0, 18, 25)',
        'Fecal_2 + Tomato': 'rgb(187, 62, 3)',
        'Fecal_2 + Plasma': 'rgb(72, 202, 228)',
        'Tomato + Plasma': 'rgb(238, 155, 0)',
    }
    x_legend = 0.1
elif (plot == 'Proportion_Fecal_1') | (plot == 'Proportion_NIST_1950_SRM') | (plot == 'Proportion_Tomato')|(plot == 'Proportion_Fecal_2') :
    colorsIdx = {
        '0': 'rgb(233, 216, 166)',
        '20': 'rgb(144, 224, 239)',
        '25': 'rgb(72, 202, 228)',
        '40': 'rgb(0, 180, 216)',
        '50': 'rgb(0, 150, 199)',
        '75': 'rgb(0, 119, 182)',
        '100': 'rgb(2, 62, 138)'
    }
else:
    colorsIdx = {
        'C18': 'rgb(238, 155, 0)',
        'C18-RTshift': 'rgb(0, 95, 115)'
    }
    x_legend = 0.42

fig = make_subplots(rows=2, cols=4,
                    shared_xaxes=False,
                    vertical_spacing=0.13,
                    horizontal_spacing=0.1,
                    subplot_titles=('Feature Table (Bray-Curtis)', 'MEMO unaligned (Bray-Curtis)', 'Qemistree weighted UniFrac', 'CSCS weighted', '', 'MEMO aligned (Bray-Curtis)', 'Qemistree unweighted UniFrac', 'CSCS unweighted'))

cats = list(df_meta[plot].unique())

# Plot  lines

result_cat = results.copy()
for metric, position in zip(metrics, pos):
    for sample in list(result_cat['Samplename'].unique()):
        for triplicate in [1,2,3]:
            from_x = list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['Experiment'] == 'C18') & (result_cat['Triplicate_number'] == triplicate)][metric+'_PC1'])[0]
            from_y = list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['Experiment'] == 'C18') & (result_cat['Triplicate_number'] == triplicate)][metric+'_PC2'])[0]
            to_x =  list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['Experiment'] == 'C18-RTshift') & (result_cat['Triplicate_number'] == triplicate)][metric+'_PC1'])[0]
            to_y =  list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['Experiment'] == 'C18-RTshift') & (result_cat['Triplicate_number'] == triplicate)][metric+'_PC2'])[0]

            fig.add_trace(
                go.Scatter(
                    x= [from_x, to_x],
                    y = [from_y, to_y],
                    mode='lines', showlegend = False,
                    line=dict(color='black', width=0.07)
                    ),
                row=position[0], col=position[1])

for cat in cats:  
    result_cat = results[results[plot] == cat]
    for metric, position, legend_bool in zip(metrics, pos, legend):
        fig.add_trace(
            go.Scatter(
                marker=dict(size=list(result_cat['size']), line_width=list(result_cat['line']), line_color ='black', opacity=0.95), 
                text=result_cat['Samplename'],
                x= result_cat[metric+'_PC1'], y = result_cat[metric+'_PC2'],
                mode='markers', marker_color= colorsIdx[cat],
                name=cat, legendgroup=cat, showlegend=legend_bool
                ),
            row=position[0], col=position[1])


# Update xaxis properties
fig.update_xaxes(title_text=f"PC1 ({results_expvar['Feature Table (Bray-Curtis)']['varPC1']} %)", row=1, col=1)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['MEMO unaligned (Bray-Curtis)']['varPC1']} %)", row=1, col=2)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['Qemistree weighted UniFrac']['varPC1']} %)", row=1, col=3)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['CSCS weighted']['varPC1']} %)", row=1, col=4)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['MEMO aligned (Bray-Curtis)']['varPC1']} %)", row=2, col=2)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['Qemistree unweighted UniFrac']['varPC1']} %)", row=2, col=3)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['CSCS unweighted']['varPC1']} %)", row=2, col=4)

# Update yaxis properties
fig.update_yaxes(title_text=f"PC2 ({results_expvar['Feature Table (Bray-Curtis)']['varPC2']} %)", row=1, col=1)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['MEMO unaligned (Bray-Curtis)']['varPC2']} %)", row=1, col=2)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['Qemistree weighted UniFrac']['varPC2']} %)", row=1, col=3)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['CSCS weighted']['varPC2']} %)", row=1, col=4)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['MEMO aligned (Bray-Curtis)']['varPC2']} %)", row=2, col=2)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['Qemistree unweighted UniFrac']['varPC2']} %)", row=2, col=3)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['CSCS unweighted']['varPC2']} %)", row=2, col=4)

fig.update_layout(height=900, width=1800, template = 'simple_white')

fig.update_layout(legend=dict(
    orientation="h",
    font=dict(
        size=18,
        color="black"
        ),
))

fig.update_annotations(font_size=20)

fig.show()

fig.write_html(f"plot/qemistree_dataset_color_{plot}.html")
fig.write_image(f"plot/qemistree_dataset_color_{plot}.svg")
fig.write_image(f"plot/qemistree_dataset_color_{plot}.jpeg")


'''PERMANOVA'''


group = 'Experiment' # Experiment contains

sk_bio_dm_bray = DistanceMatrix(dm_braycurtis, samples)
sk_bio_dm_memo = DistanceMatrix(dm_memo, samples)
sk_bio_dm_memo2 = DistanceMatrix(dm_memo2, samples)
sk_bio_dm_qemistree_w = DistanceMatrix(dm_qemistree_weighted, samples)
sk_bio_dm_qemistree_u = DistanceMatrix(dm_qemistree_unweighted, samples)
sk_bio_dm_cscs_w = DistanceMatrix(dm_cscs_weighted, samples)
sk_bio_dm_cscs_u = DistanceMatrix(dm_cscs_unweighted, samples)

metrics = [
    'Feature Table (Bray-Curtis)', 'MEMO aligned (Bray-Curtis)', 'MEMO unaligned (Bray-Curtis)' , 'Qemistree weighted UniFrac',
    'Qemistree unweighted UniFrac', 'CSCS weighted', 'CSCS unweighted'
    ]

dms = [
    sk_bio_dm_bray, sk_bio_dm_memo, sk_bio_dm_memo2, sk_bio_dm_qemistree_w,
    sk_bio_dm_qemistree_u, sk_bio_dm_cscs_w, sk_bio_dm_cscs_u
    ]

permanova_results = {}
for metric, dm in zip(metrics, dms):
    permanova_results[metric] = permanova(dm, df_meta, column=group, permutations=999)

results = pd.DataFrame(permanova_results).transpose().rename(columns={'test statistic': 'pseudo-F'})[['pseudo-F', 'p-value']]
results

results.to_csv(f"permanova/qemistree_dataset_group_{group}.csv")


''' Co-analysis of q-ToF and Q-Exactive data '''


df_meta = pd.read_csv("../01_input_data/02_qemistree/1901_gradient_benchmarking_dataset_v4_sample_metadata.txt", sep='\t')
df_meta['Samplename'] = df_meta['Samplename'].str[:-6]
df_meta['Samplename'] = df_meta['Samplename'].str.replace('BLANK_', 'BLANK')
df_meta = df_meta[['Filename', 'Experiment', 'Samplename', 'Triplicate_number', 'Proportion_Fecal_1', 'Proportion_Fecal_2', 'Proportion_Tomato', 'Proportion_NIST_1950_SRM']]
df_meta['contains'] = df_meta.apply(conditions, axis=1)
df_meta['instrument'] = np.where(df_meta['Samplename'].str.contains('qTOF'), 'qTOF', 'QE')
df_meta['blank_qc'] = np.where(df_meta['Samplename'].str.contains('blank|qcmix', case = False), 'yes', 'no')
df_meta

    # QE
feat_table_qe = memo.FeatureTable(path="../01_input_data/02_qemistree/qe_qtof_coanalysis/qe_quant_nogapF.csv", software='mzmine')
spectra_qe = memo.SpectraDocuments(path="../01_input_data/02_qemistree/qe_qtof_coanalysis/qe_spectra_nogapF.mgf", min_relative_intensity = 0.01,
            max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)
memo_qe = memo.MemoMatrix()
memo_qe.memo_from_aligned_samples(feat_table_qe, spectra_qe)

    # Qtof
feat_table_qtof = memo.FeatureTable(path="../01_input_data/02_qemistree/qe_qtof_coanalysis/qtof_quant_nogapF.csv", software='mzmine')
spectra_qtof = memo.SpectraDocuments(path="../01_input_data/02_qemistree/qe_qtof_coanalysis/qtof_spectra_nogapF.mgf", min_relative_intensity = 0.01,
            max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)

memo_qtof = memo.MemoMatrix()
memo_qtof.memo_from_aligned_samples(feat_table_qtof, spectra_qtof)

memo_qe = memo_qe.filter(samples_pattern='blank|Qcmix')
memo_qtof = memo_qtof.filter(samples_pattern='blank|Qcmix')

memo_merged = memo_qe.merge_memo(memo_qtof, drop_not_in_common=False)

    # Filter and sort metadata + memo matrix
df_meta = df_meta[df_meta['Filename'].isin(list(memo_merged.memo_matrix.index))]
memo_merged.memo_matrix = memo_merged.memo_matrix.loc[list(df_meta['Filename'])]

    # Sample-wise normalization
memo_merged.memo_matrix = memo_merged.memo_matrix.div(memo_merged.memo_matrix.sum(axis=1), axis=0)

dm_memo = sp.spatial.distance.pdist(memo_merged.memo_matrix , 'braycurtis') 


''' Plot '''


results = pd.DataFrame()
results_expvar = {}
PCs = [['PC1', 'PC2'], ['PC2','PC3']]
pos = [[1,1],[1,2]]
legend = [True, False]

pcoa_results = pcoa(dm_memo)

results['PC1'] = pcoa_results.samples['PC1']
results['PC2'] = pcoa_results.samples['PC2']
results['PC3'] = pcoa_results.samples['PC3']

varPC1 = round(100*pcoa_results.proportion_explained[0], 1)
varPC2 = round(100*pcoa_results.proportion_explained[1], 1)
varPC3 = round(100*pcoa_results.proportion_explained[2], 1)

results['contains'] = list(df_meta['contains'])
results['instrument'] = list(df_meta['instrument'])
results['Samplename'] = list(df_meta['Samplename'])
results['Samplename'] = results['Samplename'].str.replace('QE_C18_', '')
results['Samplename'] = results['Samplename'].str.replace('qTOF_C18_', '')
results['Triplicate_number'] = list(df_meta['Triplicate_number'])
results['size'] = np.where(results['contains'].isin(['Fecal_2', 'Fecal_1', 'Tomato', 'Plasma']), 11, 7)
results['line'] = np.where(results['contains'].isin(['Fecal_2', 'Fecal_1', 'Tomato', 'Plasma']), 0.8, 0)

    # Choose cat to plot
plot = 'instrument' # instrument contains

if plot == 'contains':
    colorsIdx = {

        'Fecal_1 + Fecal_2': 'rgb(10, 147, 150)',
        'Fecal_1 + Tomato': 'rgb(202, 103, 2)',
        'Fecal_1 + Plasma': 'rgb(0, 150, 199)',
        'Fecal_1 + Fecal_2 + Tomato + Plasma': 'rgb(0, 18, 25)',
        'Fecal_2 + Tomato': 'rgb(187, 62, 3)',
        'Fecal_2 + Plasma': 'rgb(72, 202, 228)',
        'Tomato + Plasma': 'rgb(238, 155, 0)',
        'Fecal_1': 'rgb(0, 95, 115)',
        'Fecal_2': 'rgb(148, 210, 189)',
        'Plasma': 'rgb(233, 216, 166)',
        'Tomato': 'rgb(155, 34, 38)',
    }
else:
    colorsIdx = {
        'QE': 'rgb(238, 155, 0)',
        'qTOF': 'rgb(0, 95, 115)'
    }

fig = make_subplots(rows=1, cols=2,
                    shared_xaxes=False,
                    vertical_spacing=0.13,
                    horizontal_spacing=0.2,
                    subplot_titles=('PC1 and PC2', 'PC2 and PC3'))

cats = list(colorsIdx.keys())

# Plot  lines

result_cat = results.copy()
for pc, position in zip(PCs, pos):
    for sample in list(result_cat['Samplename'].unique()):
        for triplicate in [1,2,3]:
            from_x = list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['instrument'] == 'QE') & (result_cat['Triplicate_number'] == triplicate)][pc[0]])[0]
            from_y = list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['instrument'] == 'QE') & (result_cat['Triplicate_number'] == triplicate)][pc[1]])[0]
            to_x =  list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['instrument'] == 'qTOF') & (result_cat['Triplicate_number'] == triplicate)][pc[0]])[0]
            to_y =  list(result_cat[(result_cat['Samplename'] == sample) & (result_cat['instrument'] == 'qTOF') & (result_cat['Triplicate_number'] == triplicate)][pc[1]])[0]

            fig.add_trace(
                go.Scatter(
                    x= [from_x, to_x],
                    y = [from_y, to_y],
                    mode='lines', showlegend = False,
                    line=dict(color='black', width=0.07)
                    ),
                row=position[0], col=position[1])

for cat in cats:  
    result_cat = results[results[plot] == cat]
    for pc, position, legend_bool in zip(PCs, pos, legend):
        fig.add_trace(
            go.Scatter(
                marker=dict(size=list(result_cat['size']), line_width=list(result_cat['line']), line_color ='black', opacity=0.95), 
                x= result_cat[pc[0]], y = result_cat[pc[1]],
                mode='markers', marker_color= colorsIdx[cat], name=cat, legendgroup=cat, showlegend=legend_bool
                ),
            row=position[0], col=position[1])

# Update xaxis properties
fig.update_xaxes(title_text=f"PC1 ({varPC1} %)", row=1, col=1)
fig.update_xaxes(title_text=f"PC2 ({varPC2} %)", row=1, col=2)

# Update yaxis properties
fig.update_yaxes(title_text=f"PC2 ({varPC2} %)", row=1, col=1)
fig.update_yaxes(title_text=f"PC3 ({varPC3} %)", row=1, col=2)

fig.update_layout(height=550, width=1000, template = 'simple_white')

fig.update_layout(legend=dict(
    orientation="h",
    yanchor = 'bottom', 
    y = -0.5,
    font=dict(
        size=12,
        color="black"
        ),
))
fig.update_annotations(font_size=20)

fig.show()

fig.write_html(f"plot/qemistree_dataset_qe_vs_qtof_color_{plot}.html")
fig.write_image(f"plot/qemistree_dataset_qe_vs_qtof_color_{plot}.svg")
fig.write_image(f"plot/qemistree_dataset_qe_vs_qtof_color_{plot}.jpeg")


'''PERMANOVA'''


group = 'instrument' # instrument contains

samples = memo_merged.memo_matrix.index.to_list()

sk_bio_dm_memo = DistanceMatrix(dm_memo, samples)

df_meta_2 = df_meta.set_index('Filename')

permanova(sk_bio_dm_memo, df_meta_2, column='instrument', permutations=999)
permanova(sk_bio_dm_memo, df_meta_2, column='contains', permutations=999)

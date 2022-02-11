import memo_ms as memo
import pandas as pd
import datatable as dt
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import time
import scipy as sp
import umap
from skbio.stats.ordination import pcoa
import plotly.graph_objects as go
import pickle

########################################################################################################################

def conditions(df_meta):
    if (df_meta['sample_plate_id'] in ['VGF138', 'VGF139','VGF140','VGF141','VGF142', 'VGF143', 'VGF144', 'VGF145', 'VGF146', 'VGF147']):
        return 'Batch 1'
    elif (df_meta['sample_plate_id'] in ['VGF150', 'VGF151','VGF152','VGF153','VGF154', 'VGF155', 'VGF156', 'VGF157', 'VGF158', 'VGF159']):
        return 'Batch 2'
    else:
        return 'What is it?'

########################################################################################################################


""" PART 1: IMPORT AND PREPARE METADATA"""


df_meta = pd.read_csv("../01_input_data/03_plant_extract_dataset/metadata/plant_extract_dataset_metadata.tsv", sep='\t')
df_meta.sort_values(['ms_filename'], inplace=True)

df_meta['batch'] = df_meta.apply(conditions, axis=1)

df_meta['ms_injection_date_dateformat'] = pd.to_datetime(df_meta['ms_injection_date'], format = '%Y-%m-%d', infer_datetime_format=True)
df_meta['before_after'] = np.where(df_meta.ms_injection_date_dateformat < '2017-10-23', 'Before 2017-10-23', 'After 2017-10-26')         

df_meta_samples = df_meta[df_meta['sample_type'] == 'sample']
df_meta_qc = df_meta[df_meta['sample_type'] == 'qc']
df_meta_blanks = df_meta[df_meta['sample_type'] == 'blank']

df_meta_samples['goniothalamus'] = np.where(df_meta_samples['organism_species'] == 'Goniothalamus gabriacianus', 'Goniothalamus gabriacianus', 'Other')
df_meta_samples['genus_selected'] = np.where((df_meta_samples['organism_genus'] == 'Limonia')|(df_meta_samples['organism_genus'] == 'Mangifera')|(df_meta_samples['organism_genus'] == 'Zanthoxylum'), df_meta_samples['organism_genus'], 'Other')
df_meta_samples['species_selected'] = np.where(
    (df_meta_samples['organism_species'] == 'Limonia acidissima')|(df_meta_samples['organism_species'] == 'Baliospermum solanifolium')|
    (df_meta_samples['organism_species'] == 'Bombax anceps')|(df_meta_samples['organism_species'] == 'Cerbera manghas')|
    (df_meta_samples['organism_species'] == 'Goniothalamus gabriacianus'),
    df_meta_samples['organism_species']
    , 'Other')

df_meta_samples['species_organ'] = df_meta_samples['organism_species'] + ' ' + df_meta_samples['organism_organ']
df_meta_samples['species_organ_selected'] = np.where(
    (df_meta_samples['species_organ'] == 'Hymenocardia punctata leaves')|(df_meta_samples['species_organ'] == 'Goniothalamus gabriacianus roots')|
    (df_meta_samples['species_organ'] == 'Limonia acidissima leaves')|(df_meta_samples['species_organ'] == 'Bombax anceps roots')|(df_meta_samples['species_organ'] == 'Bombax anceps leaves')|
    (df_meta_samples['species_organ'] == 'Pithecellobium dulce leaves')|(df_meta_samples['species_organ'] == 'Goniothalamus gabriacianus leaves')|(df_meta_samples['species_organ'] == 'Flacourtia jangomas leaves'),
    df_meta_samples['species_organ']
    , 'Other')

df_meta_samples['species_organ'] = df_meta_samples['organism_species'] + ' ' + df_meta_samples['organism_organ']
df_meta_samples['species_organ_selected_size'] = np.where(
    df_meta_samples['species_organ_selected'] == 'Other',
    1, 10)

dic_cat = {}
categories = ['species_organ_selected', 'before_after', 'genus_selected', 'goniothalamus', 'species_selected', 'ms_injection_date', 'tcruzi_activity_class']

for cat in categories:
    if cat == 'before_after':
        categorical = True
        colorsIdx = {
            'Before 2017-10-23': 'rgb(238, 155, 0)', #'rgb(174, 32, 18)',
            'After 2017-10-26': 'rgb(0, 95, 115)'
        }
        title = 'Injection date'
    elif cat == 'goniothalamus':
        categorical = True
        colorsIdx = {
        'Goniothalamus gabriacianus': 'rgb(0, 95, 115)',
        'Other': 'rgb(233, 216, 166)'
        }
        title = 'Species'
    elif cat == 'genus_selected':
        categorical = True
        colorsIdx = {
        'Mangifera': 'rgb(10, 147, 150)',
        'Limonia': 'rgb(238, 155, 0)',
        'Zanthoxylum': 'rgb(187, 62, 3)',
        'Other': 'rgb(233, 216, 166)'
        }
        title = 'Genus'
    elif cat == 'species_selected':
        categorical = True
        colorsIdx = {
        'Limonia acidissima': 'rgb(0, 18, 25)',
        'Baliospermum solanifolium': 'rgb(10, 147, 150)',
        'Goniothalamus gabriacianus': 'rgb(238, 155, 0)',
        'Cerbera manghas': 'rgb(187, 62, 3)',
        'Bombax anceps': 'rgb(148, 210, 189)',
        'Other': 'rgb(233, 216, 166)'
        }
        title = 'Species'
    elif cat == 'ms_injection_date':
        categorical = True
        colorsIdx = {
        '2017-04-14': '#006f96', '2017-04-22': '#008080',
        '2017-04-29': '#00906c', '2017-04-30': '#0f9e66',
        '2017-05-05': '#36a978', '2017-05-07': '#59b488',
        '2017-08-26': '#7abd96', '2017-10-21': '#99c7a4',
        '2017-10-22': '#b7d0b2', '2017-10-27': '#eebd00', 
        '2017-10-28': '#dda300', '2017-11-03': '#cd8900',
        '2017-11-10': '#ba7101', '2017-11-12': '#9b5f03',
        '2017-11-18': '#7e4f06', '2017-11-19': '#613e09',
        '2017-12-01': '#462f0c', '2017-12-02': '#2c200e',
        '2018-03-02': '#111111'
        }
        title = 'Injection date'

    elif cat == 'species_organ_selected':
        categorical = True
        colorsIdx = {
        'Hymenocardia punctata leaves': 'rgb(155, 34, 38)',
        'Limonia acidissima leaves': 'rgb(174, 32, 18)',
        'Bombax anceps roots': 'rgb(187, 62, 3)',
        'Bombax anceps leaves': 'rgb(202, 103, 2)',
        'Pithecellobium dulce leaves': 'rgb(238, 155, 0)',
        'Goniothalamus gabriacianus leaves': 'rgb(148, 210, 189)',
        'Goniothalamus gabriacianus roots': 'rgb(10, 147, 150)',
        'Flacourtia jangomas leaves': 'rgb(0, 95, 115)',
        'Other': 'rgb(233, 216, 166)'
        }
        title = 'Species Organ Selected'

    elif cat == 'tcruzi_activity_class':
        categorical = True
        colorsIdx = {
        'Active': 'rgb(6, 214, 160)',
        'Inactive': 'rgb(233, 216, 166)'
        }
        
    dic_cat[cat] = {}
    dic_cat[cat]['categorical'] = categorical
    dic_cat[cat]['colorsIdx'] = colorsIdx
    dic_cat[cat]['title'] = title

samples = df_meta_samples['sample_id'].to_list()

""" PART 2: GENERATE MEMO MATRIX FROM UNALIGNED SAMPLES"""

    # Path to folder containing the mgf
path_to_file = "../01_input_data/03_plant_extract_dataset/individual_mgf_files"

    # Generating memo matrix
memo_unaligned = memo.MemoMatrix()

start = time.process_time()
memo_unaligned.memo_from_unaligned_samples(path_to_file, min_relative_intensity = 0.01,
            max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)
print(f'Computing MEMO matrix from unaligned samples took: {time.process_time() - start} seconds')

memo_unaligned.memo_matrix.index = memo_unaligned.memo_matrix.index.str.replace("_features_ms2_pos", "")

memo_unaligned_filtered = memo_unaligned.filter(samples_pattern='01')
memo_unaligned_filtered = memo_unaligned_filtered.filter(samples_pattern='12', max_occurence=0)


""" PART 3: GENERATE MEMO MATRIX FROM ALIGNED SAMPLES"""

    # Path to .mgf and .csv of spectra and feature table
path_mgf = "../01_input_data/03_plant_extract_dataset/aligned_feat_table_and_spectra/210302_feature_list_filtered_strong_aligned_cosine03.mgf"
path_quant = "../01_input_data/03_plant_extract_dataset/aligned_feat_table_and_spectra/210302_feature_list_filtered_strong_aligned_cosine03_quant.csv"

start = time.process_time()
feat_table = memo.FeatureTable(path=path_quant, software='mzmine')

spectra = memo.SpectraDocuments(path=path_mgf, min_relative_intensity = 0.01,
            max_relative_intensity = 1, min_peaks_required=10, losses_from = 10, losses_to = 200, n_decimals = 2)

memo_aligned = memo.MemoMatrix()
memo_aligned.memo_from_aligned_samples(featuretable= feat_table, spectradocuments= spectra)
print(f'Computing MEMO matrix from aligned samples took: {time.process_time() - start} seconds')

memo_aligned.memo_matrix.index = memo_aligned.memo_matrix.index.str.replace(".mzXML", "")
feat_table.feature_table.index = feat_table.feature_table.index.str.replace(".mzXML", "")

memo_aligned = memo_aligned.filter(samples_pattern='01')
memo_aligned = memo_aligned.filter(samples_pattern='12', max_occurence=0)


""" SAVE/LOAD/EXTRACT MATRICES"""


    # Extract the MEMO matrix from the MemoMatrix dataclass for plotting
memo_aligned_matrix = memo_aligned.memo_matrix
memo_unaligned_matrix = memo_unaligned_filtered.memo_matrix

    # Sort matrix to match metadata
memo_aligned_matrix = memo_aligned_matrix.loc[samples]
memo_unaligned_matrix = memo_unaligned_matrix.loc[samples]

    # Remove blanks and QC from feature table abds sort to match metadata
feat_samples = feat_table.feature_table[~feat_table.feature_table.index.str.contains("01|12")].sort_index()
feat_samples = feat_samples.loc[samples]

    # Save MEMO matrices for TMAP computation
#memo_aligned.memo_matrix.to_csv('memo_matrix/memo_aligned_filtered_memo_matrix.csv')
#memo_unaligned.memo_matrix.to_csv('memo_matrix/memo_unaligned_filtered_memo_matrix.csv')

    # Save feature table for TMAP computation
#feat_samples.to_csv('feature_table/plant_extract_dataset_feature_table_samples_only.csv')

    # Load pre-computed MEMO matrices
memo_aligned_matrix = dt.fread('memo_matrix/memo_aligned_filtered_memo_matrix.csv').to_pandas()
memo_aligned_matrix.set_index('filename', inplace=True)
memo_unaligned_matrix = dt.fread('memo_matrix/memo_unaligned_filtered_memo_matrix.csv').to_pandas()
memo_unaligned_matrix.set_index('C0', inplace=True)

    # Load pre-computed feature table
feat_samples = dt.fread('feature_table/plant_extract_dataset_feature_table_samples_only.csv').to_pandas()
feat_samples.columns = feat_samples.iloc[0]
feat_samples = feat_samples[1:]
feat_samples.set_index('filename', inplace=True)


""" PLOT PCOA"""


    # distance matrix generation
print('start dm calculation')
dm_feature_samples = sp.spatial.distance.pdist(feat_samples, 'braycurtis')
print('1/3 done')
dm_memo_aligned = sp.spatial.distance.pdist(memo_aligned_matrix, 'braycurtis')
print('2/3 done')
dm_memo_unaligned = sp.spatial.distance.pdist(memo_unaligned_matrix, 'braycurtis')
print('3/3 done')

    # Convert to DataFrame and save as .csv
# df_dm_feature_samples= pd.DataFrame(sp.spatial.distance.squareform(dm_feature_samples), index= feat_samples.index, columns = feat_samples.index)
# df_dm_feature_samples.to_csv('dm_feature.csv')

# df_dm_memo_aligned = pd.DataFrame(sp.spatial.distance.squareform(dm_memo_aligned), index= memo_aligned.index, columns = memo_aligned.index)
# df_dm_memo_aligned.to_csv('distance_matrix/dm_memo_aligned.csv')

# df_dm_memo_unaligned = pd.DataFrame(sp.spatial.distance.squareform(dm_memo_unaligned), index= memo_unaligned.filtered_memo_matrix.index, columns = memo_unaligned.filtered_memo_matrix.index)
# df_dm_memo_unaligned.to_csv('distance_matrix/dm_memo_unaligned.csv')

    # Load pre-computed distance matrices
dm_feature_samples = pd.read_csv('distance_matrix/dm_feature.csv').set_index('filename')
dm_memo_aligned = pd.read_csv('distance_matrix/dm_memo_aligned.csv').set_index('filename')
dm_memo_unaligned = pd.read_csv('distance_matrix/dm_memo_unaligned.csv').set_index('C0')

# Plotting

    # 1. Generate df used for plotting  
results_pcoa = df_meta_samples.copy()
results_expvar = {}

metrics = ['Feature Table (Bray-Curtis)', 'MEMO aligned (Bray-Curtis)', 'MEMO unaligned (Bray-Curtis)']
pos = [[1,1], [1,2], [1,3]]
legend = [True, False, False]

dms = [dm_feature_samples, dm_memo_aligned, dm_memo_unaligned]

for metric, dm in zip(metrics, dms):
    pcoa_results = pcoa(dm)

    results_pcoa[metric+'_PC1'] = list(pcoa_results.samples['PC1'])
    results_pcoa[metric+'_PC2'] = list(pcoa_results.samples['PC2'])

    results_expvar[metric] = {}
    results_expvar[metric]['varPC1'] = round(100*pcoa_results.proportion_explained[0], 1)
    results_expvar[metric]['varPC2'] = round(100*pcoa_results.proportion_explained[1], 1)

    # 2. Select category to plot
category = 'ms_injection_date' # tcruzi_activity_class 'bio_tryp_cruzi_10ugml_inhibition' #before_after #genus_selected #goniothalamus #species_selected

fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    subplot_titles=('PCoA Feature Table (Bray-Curtis)', 'PCoA MEMO aligned (Bray-Curtis)', 'PCoA MEMO unaligned (Bray-Curtis)')
                    )

cats = list(results_pcoa[category].unique())

if dic_cat[category]['categorical'] is True:
    for cat in cats:  
        result_cat = results_pcoa[results_pcoa[category] == cat]
        if (category == 'species_organ_selected')|(category == 'tcruzi_activity_class'):
            if (cat == 'Other') | (cat == 'Inactive') :
                size=4
            else:
                size=10
        else:
            size = 5
        for metric, position, legend_bool in zip(metrics, pos, legend):
            fig.add_trace(
                go.Scatter(
                    x= result_cat[metric+'_PC1'], y = result_cat[metric+'_PC2'],
                    mode='markers',
                    marker=dict(
                        size = size,
                        opacity=0.85,
                        line_width=0.5,
                        line_color ='grey'
                    ),
                    marker_color= dic_cat[category]['colorsIdx'][cat], name=cat, legendgroup=cat, showlegend=legend_bool,
                    hovertext = result_cat['sample_id'], line_width=1,
                    line_color ='grey'
                    ),
                row=position[0], col=position[1])

if dic_cat[category]['categorical'] is False:
    for metric, position, legend_bool in zip(metrics, pos, legend):
        fig.add_trace(
            go.Scatter(
                opacity=0.85,
                mode='markers',
                x= results_pcoa[metric+'_PC1'], y = results_pcoa[metric+'_PC2'],
                marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='lightgray',
                    color=results_pcoa[category], 
                    colorscale='YlOrRd', 
                    showscale=True
                ),
                name=category, showlegend=legend_bool,
                hovertext = results_pcoa['organism_species']
                ),
            row=position[0], col=position[1])

# Update xaxis properties
fig.update_xaxes(title_text=f"PC1 ({results_expvar['Feature Table (Bray-Curtis)']['varPC1']} %)", row=1, col=1)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['MEMO aligned (Bray-Curtis)']['varPC1']} %)", row=1, col=2)
fig.update_xaxes(title_text=f"PC1 ({results_expvar['MEMO unaligned (Bray-Curtis)']['varPC1']} %)", row=1, col=3)

# Update yaxis properties
fig.update_yaxes(title_text=f"PC2 ({results_expvar['Feature Table (Bray-Curtis)']['varPC2']} %)", row=1, col=1)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['MEMO aligned (Bray-Curtis)']['varPC2']} %)", row=1, col=2)
fig.update_yaxes(title_text=f"PC2 ({results_expvar['MEMO unaligned (Bray-Curtis)']['varPC2']} %)", row=1, col=3)


fig.update_layout(height=600, width=1500, template = 'simple_white')

fig.update_layout(
    legend=dict(
        title = dic_cat[category]['title'],
        orientation="h",
        y=-0.3,
        x = -0.1,
        bordercolor="Black",
        borderwidth=2
        ),
    font=dict(
        family="Arial",
        size=16,
        color="black"
        )    
    )

fig.update_annotations(font_size=20)

fig.show()

    # Save plots 
fig.write_html(f"plot/pcoa/pcoa_vgf_color_{category}.html")
fig.write_image(f"plot/pcoa/pcoa_vgf_color_{category}.jpeg",  scale=3)
fig.write_image(f"plot/pcoa/pcoa_vgf_color_{category}.svg",  scale=3)


""" PLOT: UMAP """


    # Sample-wise normalization
memo_aligned_matrix = memo_aligned_matrix.div(memo_aligned_matrix.sum(axis=1), axis=0)
memo_unaligned_matrix = memo_unaligned_matrix.div(memo_unaligned_matrix.sum(axis=1), axis=0)

# Plotting

    # 1. Generate df used for plotting 
results_umap = df_meta_samples.copy()

metrics = ['Feature Table (Bray-Curtis)', 'MEMO aligned (Bray-Curtis)', 'MEMO unaligned (Bray-Curtis)']
pos = [[1,1], [1,2], [1,3]]
legend = [True, False, False]

matrices = [feat_samples, memo_aligned_matrix, memo_unaligned_matrix]

for metric, matrix in zip(metrics, matrices):
    reducer = umap.UMAP(metric='braycurtis')
    embedding = reducer.fit_transform(matrix.values)

    results_umap[metric+'_x'] = embedding[:, 0]
    results_umap[metric+'_y'] = embedding[:, 1]

    # 2. Select category to plot
category = 'tcruzi_activity_class' #tcruzi_activity_class' #'bio_tryp_cruzi_10ugml_inhibition' #before_after #genus_selected #goniothalamus #species_selected

fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    subplot_titles=('UMAP Feature Table', 'UMAP MEMO aligned', 'UMAP MEMO unaligned'))

cats = list(results_umap[category].unique())

if dic_cat[category]['categorical'] is True:
    for cat in cats:  
        result_cat = results_umap[results_umap[category] == cat]
        if (category == 'species_organ_selected')|(category == 'tcruzi_activity_class'):
            if (cat == 'Other') | (cat == 'Inactive') :
                size=4
            else:
                size=10
        else:
            size = 5
        for metric, position, legend_bool in zip(metrics, pos, legend):
            fig.add_trace(
                go.Scatter(
                    x= result_cat[metric+'_x'], y = result_cat[metric+'_y'],
                    mode='markers',
                    marker=dict(
                        size = size,
                        opacity=0.85,
                        line_width=0.5,
                        line_color ='grey'
                    ),
                    marker_color= dic_cat[category]['colorsIdx'][cat], name=cat, legendgroup=cat, showlegend=legend_bool,
                    hovertext = result_cat['organism_species'], line_width=1,
                    line_color ='grey'
                    ),
                row=position[0], col=position[1])

if dic_cat[category]['categorical'] is False:
    for metric, position, legend_bool in zip(metrics, pos, legend):
        fig.add_trace(
            go.Scatter(
                opacity=0.85,
                mode='markers',
                x= results_umap[metric+'_x'], y = results_umap[metric+'_y'],
                marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='lightgray',
                    color=results_umap[category], 
                    colorscale='YlOrRd', 
                    showscale=True
                ),
                name=category, showlegend=legend_bool,
                hovertext = results_umap['organism_species']
                ),
            row=position[0], col=position[1])

fig.update_layout(height=600, width=1500, template = 'simple_white')

fig.update_layout(
    legend=dict(
        title = dic_cat[category]['title'],
        orientation="h",
        y=-0.3,
        x = -0.1,
        bordercolor="Black",
        borderwidth=2
        ),
    font=dict(
        family="Arial",
        size=16,
        color="black"
        )    
    )

fig.update_annotations(font_size=20)

fig.show()

    # Save plots 
fig.write_html(f"plot/umap/umap_vgf_color_{category}.html")
fig.write_image(f"plot/umap/umap_vgf_color_{category}.jpeg",  scale=3)
fig.write_image(f"plot/umap/umap_vgf_color_{category}.svg",  scale=3)


""" PLOT TMAP """


    # Plot TMAP with coordinates outputted from tmap_plotter.py
path_to_tmap_coordinates = [
    'tmap_coordinates/tmap_coordinates_feature_table',
    'tmap_coordinates/tmap_coordinates_memo_aligned',
    'tmap_coordinates/tmap_coordinates_memo_unaligned',
]

# Plotting

    # 1. Generate df used for plotting
results_tmap = df_meta_samples.copy()

metrics = ['Feature Table (Bray-Curtis)', 'MEMO aligned (Bray-Curtis)', 'MEMO unaligned (Bray-Curtis)']
pos = [[1,1], [1,2], [1,3]]
legend = [True, False, False]

edges = {}

for metric, path in zip(metrics, path_to_tmap_coordinates):
    x, y, s, t = pickle.load(open(path + ".dat", "rb"))  
    edges_list = []
    for source, target in zip(s,t):
        x_source = x[source]
        y_source = y[source]
        x_target = x[target]
        y_target = y[target]
        coord = ((x_source, y_source), (x_target, y_target))
        edges_list.append(coord)
        edges[metric] = edges_list

    results_tmap[metric+'_x'] = x
    results_tmap[metric+'_y'] = y

    # 2. Select category to plot
category = 'ms_injection_date' #'tcruzi_activity_class' #before_after #genus_selected #goniothalamus #species_selected tcruzi_activity_class

fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    subplot_titles=('TMAP Feature Table', 'TMAP MEMO aligned', 'TMAP MEMO unaligned'))

for metric, position, legend_bool in zip(metrics, pos, legend):
    for edge in edges[metric]:
        from_x = edge[0][0]
        from_y = edge[0][1]
        to_x = edge[1][0]
        to_y = edge[1][1]
        fig.add_trace(
            go.Scatter(
                x= [from_x, to_x],
                y = [from_y, to_y],
                mode='lines', showlegend = False,
                line=dict(color='black', width=0.07)
            ),
            row=position[0], col=position[1]
        )

cats = list(results_tmap[category].unique())

if dic_cat[category]['categorical'] is True:
    for cat in cats:  
        result_cat = results_tmap[results_tmap[category] == cat]
        if (category == 'species_organ_selected')|(category == 'tcruzi_activity_class'):
            if (cat == 'Other') | (cat == 'Inactive') :
                size=4
            else:
                size=10
        else:
            size = 5
        for metric, position, legend_bool in zip(metrics, pos, legend):
            fig.add_trace(
                go.Scatter(
                    x= result_cat[metric+'_x'], y = result_cat[metric+'_y'],
                    mode='markers',
                    marker=dict(
                        size = size,
                        opacity=0.85,
                        line_width=0.5,
                        line_color ='grey'
                    ),
                    marker_color= dic_cat[category]['colorsIdx'][cat], name=cat, legendgroup=cat, showlegend=legend_bool,
                    hovertext = result_cat['organism_species'], line_width=1,
                    line_color ='grey'
                    ),
                row=position[0], col=position[1])

elif dic_cat[category]['categorical'] is False:
    for metric, position, legend_bool in zip(metrics, pos, legend):
        fig.add_trace(
            go.Scatter(
                opacity=0.85,
                mode='markers',
                x= results_tmap[metric+'_x'], y = results_tmap[metric+'_y'],
                marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='lightgray',
                    color=results_tmap[category], 
                    colorscale='YlOrRd', 
                    showscale=True
                ),
                name=category, showlegend=legend_bool,
                hovertext = results_tmap['organism_species']
                ),
            row=position[0], col=position[1])

fig.update_layout(height=600, width=1500, template = 'simple_white')

fig.update_layout(
    legend=dict(
        title = dic_cat[category]['title'],
        orientation="h",
        y=-0.3,
        x = -0.1,
        bordercolor="Black",
        borderwidth=2
        ),
    font=dict(
        family="Arial",
        size=16,
        color="black"
        )    
    )

fig.update_annotations(font_size=20)

fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False)

fig.show()

    # Save plots 
fig.write_html(f"plot/tmap/tmap_vgf_color_{category}.html")
fig.write_image(f"plot/tmap/tmap_vgf_color_{category}.jpeg",  scale=3)
fig.write_image(f"plot/tmap/tmap_vgf_color_{category}.svg",  scale=3)


""" PART 5: GENERATE MEMO FOR WALTHERIA INDICA SAMPLES """


    # Adapt metadata

metadata = pd.read_csv("../01_input_data/03_plant_extract_dataset/metadata/plant_extract_dataset_metadata.tsv", sep='\t')
metadata_sub = metadata[['sample_id', 'sample_type', 'organism_species', 'organism_family','organism_genus', 'organism_organ']]
metadata_sub = metadata_sub[metadata_sub['sample_type'] == 'sample']
metadata_sub['species_organ'] = metadata_sub['organism_species'] + ' ' +  metadata_sub['organism_organ']

metadata_sub['species_organ_selected'] = np.where(
    metadata_sub['species_organ'] == 'Melochia umbellata green stems',
    metadata_sub['species_organ'], 'Other')

metadata_sub['species_selected'] = np.where(
    metadata_sub['organism_species'] == 'Melochia umbellata',
    metadata_sub['organism_species'], 'Other')

metadata_sub = metadata_sub.drop(['organism_organ'], axis=1) 
df_meta_waltheria = pd.read_csv("../01_input_data/03_plant_extract_dataset/metadata/metadata_waltheria.csv")
df_meta_with_waltheria = metadata_sub.append(df_meta_waltheria)

dic_cat_waltheria = {}
categories = ['species_selected', 'species_organ_selected']

for cat in categories:
    if cat == 'species_selected':
        categorical = True
        colorsIdx = {
        'Melochia umbellata': 'rgb(155, 34, 38)',
        'Waltheria indica': 'rgb(174, 32, 18)',
        'Other': 'rgb(233, 216, 166)'
        }
        title = 'Species Organ Selected'

    elif cat == 'species_organ_selected':
        categorical = True
        colorsIdx = {
        'Melochia umbellata green stems': 'rgb(69, 123, 157)',
        'Other': 'rgb(233, 216, 166)',
        'Waltheria indica aerial parts 2014': 'rgb(155, 34, 38)',
        'Waltheria indica aerial parts 2018': 'rgb(187, 62, 3)',
        'Waltheria indica roots 2018': 'rgb(238, 155, 0)',
        }
        title = 'Sample Identity'
        
    dic_cat_waltheria[cat] = {}
    dic_cat_waltheria[cat]['categorical'] = categorical
    dic_cat_waltheria[cat]['colorsIdx'] = colorsIdx
    dic_cat_waltheria[cat]['title'] = title

    # Generate MEMO matrix
    # Path to folder containning the mgf
path_to_file = "../01_input_data/03_plant_extract_dataset/waltheria_indica_mgf_files/"

    # Generating memo matrix
memo_waltheria = memo.MemoMatrix()
memo_waltheria.memo_from_unaligned_samples(path_to_file)

    # Merge with MEMO matrix of the plant extract dataset
memo_merged = memo_unaligned.merge_memo(memo_waltheria, drop_not_in_common=False)

memo_merged_filtered = memo_merged.filter(samples_pattern='01')
memo_merged_filtered = memo_merged_filtered.filter(samples_pattern='12', max_occurence=0)

memo_with_waltheria = memo_merged_filtered.memo_matrix
memo_with_waltheria = memo_with_waltheria.loc[df_meta_with_waltheria.sample_id]

    # Save as .csv
#memo_with_waltheria.to_csv('memo_matrix/memo_unaligned_filtered_memo_matrix_with_waltheria.csv')

    # Load preccomputed merged MEMO matrix
memo_with_waltheria = dt.fread('memo_matrix/memo_unaligned_filtered_memo_matrix_with_waltheria.csv').to_pandas()
memo_with_waltheria.set_index('C0', inplace=True)


""" UMAP WITH WALTHERIA"""


# Plotting

    # 1. Generate df used for plotting  
results_umap_waltheria = df_meta_with_waltheria.copy()

memo_with_waltheria = memo_with_waltheria.div(memo_with_waltheria.sum(axis=1), axis=0)

reducer = umap.UMAP(metric='braycurtis')
embedding = reducer.fit_transform(memo_with_waltheria.values)
results_umap_waltheria['umap_x'] = embedding[:, 0]
results_umap_waltheria['umap_y'] = embedding[:, 1]

    # 2. Select category to plot
category = 'species_organ_selected'

cats = list(results_umap_waltheria[category].unique())

fig = go.Figure()
for cat in cats:  
    result_cat = results_umap_waltheria[results_umap_waltheria[category] == cat]
    if (category == 'species_organ_selected')|(category == 'species_selected'):
        if (cat == 'Other') | (cat == 'Inactive') :
            size=4
        else:
            size=10
    else:
        size = 5
    fig.add_trace(
        go.Scatter(
            x= result_cat['umap_x'], y = result_cat['umap_y'],
            mode='markers',
            marker=dict(
                size = size,
                opacity=0.85,
                line_width=0.5,
                line_color ='grey'
            ),
            marker_color= dic_cat_waltheria[category]['colorsIdx'][cat], name=cat, legendgroup=cat,
            showlegend=True,
            hovertext = result_cat['species_organ'],
            line_width=1,
            line_color ='grey'
            )
        )
          
fig.update_layout(
    title = 'UMAP',
    height=600, width=600,
    template = 'simple_white',
    legend=dict(
        title = category,
        orientation="h",
        y=-0.3,#1.4
        font=dict(
            size=12,
            color="black"
            ),
    x = -0.1,
    bordercolor="Black",
    borderwidth=2
))
fig.show()

    # Save plots
fig.write_html(f"plot/waltheria/umap_vgf_with_waltheria_color_{category}.html")
fig.write_image(f"plot/waltheria/umap_vgf_with_waltheria_color_{category}.jpeg",  scale=3)
fig.write_image(f"plot/waltheria/umap_vgf_with_waltheria_color_{category}.svg",  scale=3)


""" TMAP WITH WALTHERIA """


path_to_tmap_coordinates = 'tmap_coordinates/tmap_coordinates_memo_with_waltheria'

# Plotting

    # 1. Generate df used for plotting
results_tmap_waltheria = df_meta_with_waltheria.copy()
edges = {}

x, y, s, t = pickle.load(open(path_to_tmap_coordinates + ".dat", "rb"))  
edges_list = []
for source, target in zip(s,t):
    x_source = x[source]
    y_source = y[source]
    x_target = x[target]
    y_target = y[target]
    coord = ((x_source, y_source), (x_target, y_target))
    edges_list.append(coord)

    results_tmap_waltheria['tmap_x'] = x
    results_tmap_waltheria['tmap_y'] = y

    # 2. Select category to plot
category = 'species_organ_selected'

fig = go.Figure()
for edge in edges_list:
    from_x = edge[0][0]
    from_y = edge[0][1]
    to_x = edge[1][0]
    to_y = edge[1][1]
    fig.add_trace(
        go.Scatter(
            x= [from_x, to_x],
            y = [from_y, to_y],
            mode='lines', showlegend = False,
            line=dict(color='black', width=0.07)
        )
    )

cats = list(results_tmap_waltheria[category].unique())

for cat in cats:  
    result_cat = results_tmap_waltheria[results_tmap_waltheria[category] == cat]
    if (category == 'species_organ_selected')|(category == 'species_selected'):
        if (cat == 'Other') | (cat == 'Inactive') :
            size=4
        else:
            size=10
    else:
        size = 5
    fig.add_trace(
        go.Scatter(
            x= result_cat['tmap_x'], y = result_cat['tmap_y'],
            mode='markers',
            marker=dict(
                size = size,
                opacity=0.85,
                line_width=0.5,
                line_color ='grey'
            ),
            marker_color= dic_cat_waltheria[category]['colorsIdx'][cat], name=cat, legendgroup=cat,
            showlegend=True,
            hovertext = result_cat['species_organ'],
            line_width=1,
            line_color ='grey'
            )
        )
 
fig.update_layout(
    title = 'TMAP',
    height=700, width=600,
    template = 'simple_white',
    legend=dict(
        title = category,
        orientation="h",
        y=-0.3,#1.4
        font=dict(
            size=12,
            color="black"
            ),
    x = -0.1,
    bordercolor="Black",
    borderwidth=2
))

fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False)

fig.show()

    # Save plots
fig.write_html(f"plot/waltheria/tmap_vgf_with_waltheria_color_{category}.html")
fig.write_image(f"plot/waltheria/tmap_vgf_with_waltheria_color_{category}.jpeg",  scale=3)
fig.write_image(f"plot/waltheria/tmap_vgf_with_waltheria_color_{category}.svg",  scale=3)

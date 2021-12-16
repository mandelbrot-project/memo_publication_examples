import qiime2.plugins.diversity.actions as q2_diversity
import qiime2.plugins.feature_table.actions as q2_feature_table
import qiime2 as q2

from skbio import DistanceMatrix, OrdinationResults
import pandas as pd
import itertools
from skbio import OrdinationResults
from skbio.tree import TreeNode

def load_mf(fn, index='#SampleID', sep='\t'):
    _df = pd.read_csv(fn, sep=sep, dtype='str', na_values=[], keep_default_na=False)
    _df.set_index(index, inplace=True)
    return _df

def summarize_distances(dm, metric, metadata, grouping='ATTRIBUTE_Proportions', 
                        catdiff='ATTRIBUTE_Experiment', comparison='across'):
    
    if comparison not in {'across', 'within'}:
        raise ValueError('Unsupported comparison style: %s' % comparison)
    
    groups = metadata[grouping].unique()
    pair_dist = {}
    
    for val in groups:
        # exclude blanks
        if val == '0_0_0_0':
            continue
            
        # subset to only the samples in this grouping category
        sampls = set(dm.ids).intersection(metadata[metadata[grouping] == val].index)

        if comparison == 'across':
            subset_md = metadata.loc[sampls] # only samples with grouping = val
            categories = subset_md[catdiff].unique() # catdiffs per grouping = val
            cat_samp_list = [] # list of samples per catdiff, len of this list = num_categories

            for cat in categories:
                cat_samp_list.append(list(subset_md[subset_md[catdiff] ==  cat].index))
            pair_cats = list(itertools.combinations(range(len(categories)), 2)) # pairs of categories

            for pairs in pair_cats:
                res = list(itertools.product(*[cat_samp_list[pairs[0]],cat_samp_list[pairs[1]]])) # all pairs across categories
                for samp_pair in res:
                    pair_dist[samp_pair] = [dm[samp_pair], val, metric, 'across']
        elif comparison == 'within':
            
            for _, _df in metadata.loc[sampls].groupby(catdiff):
            
                res = list(itertools.combinations(_df.index, 2)) # all pairs for a given proportion
                
                for samp_pair in res:
                    pair_dist[samp_pair] = [dm[samp_pair], val, metric, 'within']
            

    return pd.DataFrame.from_dict(data=pair_dist, orient='index',
                                  columns=['distance', 'grouping', 'metric', 'comparison'])
    
mapping = load_mf('data/qe_qert_qiime2_metadata.tsv', 'sample_name')

[i for i in mapping.index if 'SPE-blanc' in i]

_mapping = mapping.copy()
results = {}

data_locations = [
    ('ATTRIBUTE_Experiment',
     'data/qe_qert_qiime2_metadata.tsv',
     'data/ftable_qe_qert.qza',
     #'data/210721_output.qza',
     'data/tree_qe_qert.qza'),
]

for benchmark, mapping_path, table_path, tree_path in  data_locations: 
    
    table = q2.Artifact.load(table_path)
    tree = q2.Artifact.load(tree_path)
    _mapping = load_mf(mapping_path, 'sample_name')
    
    df = table.view(pd.DataFrame)
    
    results[benchmark] = {}

    # List the number of blanks in this dataset:
    # there's some samples that aren't identified as blanks but that have to be removed
    blanks = set(_mapping.query('ATTRIBUTE_Sample_type != "Sample"').index)
    blanks = set(df.index.intersection(blanks))
    
    samples = set(df.index) - set(blanks)

    print('[%s] blanks: %d, samples: %d' % (benchmark, len(blanks), len(samples)))

    # Add a new column that describes if the sample in this dataset is a blank.

    _mapping['is_blank'] = "False"
    _mapping.loc[blanks, 'is_blank'] = "True"

    _mapping.replace(to_replace='', value='NA', inplace=True)

    # first remove the blank samples, note this depends in a column that was created a few cells above
    _table, = q2_feature_table.filter_samples(table=table, metadata=q2.Metadata(_mapping), where='is_blank = "False"')
    
    results[benchmark]['mapping'] = _mapping.loc[samples].copy()
    results[benchmark]['table'] = _table.view(pd.DataFrame)
    results[benchmark]['tree'] = tree.view(TreeNode)

    for metric in ['Bray-Curtis', 'Weighted UniFrac', 'Unweighted UniFrac']:

        if metric == 'Bray-Curtis':
            dm, = q2_diversity.beta(table=_table, metric='braycurtis')
        elif metric == 'Weighted UniFrac':
            dm, = q2_diversity.beta_phylogenetic(table=_table, phylogeny=tree,
                                                 metric='weighted_normalized_unifrac') #'weighted_normalized_unifrac'
        elif metric == 'Unweighted UniFrac':
            dm, = q2_diversity.beta_phylogenetic(table=_table, phylogeny=tree,
                                                 metric='unweighted_unifrac')                          
        else:
            raise ValueError("Not expecting this metric %s" % metric)

        results[benchmark][metric] = {}
        results[benchmark][metric]['dm'] = dm.view(DistanceMatrix)
        results[benchmark][metric]['dm_view'] = dm.view(DistanceMatrix).to_data_frame()
        results[benchmark][metric]['pcoa'] = q2_diversity.pcoa(dm)[0].view(OrdinationResults)
    
results[benchmark]['Bray-Curtis']['dm_view'].to_csv('/mnt/c/users/gaudrya.FARMA/Desktop/qemistree_benchmark/bray_curtis_dm.csv')
results[benchmark]['Weighted UniFrac']['dm_view'].to_csv('/mnt/c/users/gaudrya.FARMA/Desktop/qemistree_benchmark/weighted_unifrac_dm.csv')
results[benchmark]['Unweighted UniFrac']['dm_view'].to_csv('/mnt/c/users/gaudrya.FARMA/Desktop/unweighted_unifrac_dm.csv')
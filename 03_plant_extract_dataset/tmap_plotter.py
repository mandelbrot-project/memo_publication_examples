from os import path
import pandas as pd
import tmap as tm
import pickle
import datatable as dt


##############################################################################################

def compute_tmap(table):
    
    dims = 1024
    enc = tm.Minhash(len(table.columns), 42, dims)

    i = 0
    fps = []
    for _, row in table.iterrows():
        if i != 0 and i % 10 == 0:
                print(100 * i / len(table))
        fps.append(tm.VectorFloat(list(row))) #VectorFloat or VectorUint
        i+=1

    lf = tm.LSHForest(dims * 2, 128, store=True, weighted=True)
    lf.batch_add(enc.batch_from_weight_array(fps, method="I2CWS")) #int
    lf.index()

    CFG_TMAP = tm.LayoutConfiguration()
    CFG_TMAP.k = 10
    CFG_TMAP.kc = 10

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, CFG_TMAP)
    return x, y, s, t

def store_coordinates(x, y, s, t, path):
        # To store coordinates
    x = list(x)
    y = list(y)
    s = list(s)
    t = list(t)
    pickle.dump(
        (x, y, s, t), open("tmap_coordinates/" + path, "wb+"), protocol=pickle.HIGHEST_PROTOCOL
    )
    return None

##############################################################################################

    # Load MEMO for plant extract dataset
memo_aligned_filtered_memo_matrix = dt.fread('memo_matrix/memo_aligned_filtered_memo_matrix.csv').to_pandas()
memo_aligned_filtered_memo_matrix.set_index('filename', inplace=True)

memo_unaligned_filtered_memo_matrix = dt.fread('memo_matrix/memo_unaligned_filtered_memo_matrix.csv').to_pandas()
memo_unaligned_filtered_memo_matrix.set_index('C0', inplace=True)

    # Load Feature Table for plant extract dataset
feature_table = dt.fread('feature_table/plant_extract_dataset_feature_table_samples_only.csv').to_pandas()
feature_table.columns = feature_table.iloc[0]
feature_table = feature_table[1:]
feature_table.set_index('filename', inplace=True)

    # Load MEMO for plant extract dataset with Waltheria
memo_vgf_waltheria = dt.fread('memo_matrix/memo_unaligned_filtered_memo_matrix_with_waltheria.csv').to_pandas()
memo_vgf_waltheria.set_index('C0', inplace=True)
memo_vgf_waltheria = memo_vgf_waltheria.div(memo_vgf_waltheria.sum(axis=1), axis=0)

tables = [memo_aligned_filtered_memo_matrix, memo_unaligned_filtered_memo_matrix, feature_table, memo_vgf_waltheria]
output_paths = ['tmap_coordinates_memo_aligned.dat', 'tmap_coordinates_memo_unaligned.dat', 'tmap_coordinates_feature_table.dat', 'tmap_coordinates_memo_with_waltheria.dat']

for table, path in zip(tables, output_paths):
    x, y, s, t = compute_tmap(table)
    store_coordinates(x, y, s, t, path)

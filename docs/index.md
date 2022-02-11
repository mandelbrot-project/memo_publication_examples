---
layout: default
---
![memo_logo](memo_logo.jpg)


This website gathers interactive visualizations of the plots presented in MEMO article's figures.

MEMO analysis is made available to the community through a python package. For any suggestion/remark/bug, you can raise an issue on the package github page [here](https://github.com/mandelbrot-project/memo). 

# Interactive plots

## Evaluation dataset: Q-Exactive RT vs Q-Exactive RT-shift

Classical feature-table (Bray-Curtis), MEMO, Qemistree and CSCS comparison using PCoA visualization on the [Qemistree Evaluation Dataset](https://www.nature.com/articles/s41589-020-00677-3).

Comparative PCoA with samples colored according to their:

*   [**Content**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_color_contains.html)
*   [**LC experiment**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_color_Experiment.html)
*   [**Proportion of Fecal-1**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_color_Proportion_Fecal_1.html)
*   [**Proportion of Fecal-2**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_color_Proportion_Fecal_2.html)
*   [**Proportion of Tomato**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_color_Proportion_Fecal_1.html)
*   [**Proportion of Plasma**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_color_Proportion_NIST_1950_SRM.html)

## Evaluation dataset: Q-Exactive vs Q-ToF

Co-analysis of the same samples as above, analyzed with the same LC-method on **2 different Mass Spectrometers (MS)**.

PCoA (PC1 vs PC2 _and_ PC2 vs PC3) with samples colored according to their:

*   [**Content**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_qe_vs_qtof_color_contains.html)
*   [**MS instrument**](https://mandelbrot-project.github.io/memo_publication_examples/benchmark/qemistree_dataset_qe_vs_qtof_color_instrument.html)
  
## Plant Extracts Dataset (1600 plant extracts)
Comparison of 3 different visualization (PCoA, UMAP and TMAP) performed on 3 different matrices: 
*    the feature-table (Bray-Curtis)
*    MEMO from aligned samples matrix
*    and MEMO from unaligned samples matrix
on the plant extract dataset, consisting of 1600 plant extracts-

### PCoA

Comparative PCoA with samples colored according to their:

*   [**Injection date (2 groups)**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/pcoa_vgf_color_before_after.html)
*   [**Injection date**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/pcoa_vgf_color_ms_injection_date.html)
*   [**_Trypanosoma cruzi_ activity**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/pcoa_vgf_color_tcruzi_activity_class.html)

### UMAP

Comparative UMAP with samples colored according to their:

*   [**Injection date (2 groups)**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/umap_vgf_color_before_after.html)
*   [**Injection date**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/umap_vgf_color_ms_injection_date.html)
*   [**_Trypanosoma cruzi_ activity**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/umap_vgf_color_tcruzi_activity_class.html)

### TMAP

Comparative TMAP with samples colored according to their:

*   [**Injection date (2 groups)**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/tmap_vgf_color_before_after.html)
*   [**Injection date**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/tmap_vgf_color_ms_injection_date.html)
*   [**_Trypanosoma cruzi_ activity**](https://mandelbrot-project.github.io/memo_publication_examples/plant_extract_dataset/tmap_vgf_color_tcruzi_activity_class.html)

## Plant Extract Dataset and three _Waltheria indica_(Malvaceae) samples (1603 plant extracts)
UMAP and TMAP visulalization of MEMO from unaligned samples matrix of the Plant extract dataset after addition of 3 _Waltheria indica_ samples: 2 extracts of aerial parts and roots analyzed with same method as the plant extract dataset, and 1 extract of aerial parts analyzed with a different LC-method and Orbitrap (Plus vs Q-Exactive).

*   [**UMAP**](https://mandelbrot-project.github.io/memo_publication_examples/waltheria_indica/umap_vgf_with_waltheria_color_species_organe_selected.html)

*   [**TMAP**](https://mandelbrot-project.github.io/memo_publication_examples/waltheria_indica/tmap_vgf_with_waltheria_color_species_organe_selected.html)

## Credits
Plots were generated using [Plotly](https://plot.ly.) (Plotly Technologies Inc. Collaborative data science. Montréal, QC, 2015)

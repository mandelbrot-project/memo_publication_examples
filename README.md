# memo_publication_examples

![memo_logo](docs/memo_logo.jpg)

# Code and Data to reproduce analyses of MEMO paper

## Interactive HTML plots
Interactive plots are available [here](https://mandelbrot-project.github.io/memo_publication_examples/).

# Steps to Reproduce

### 1.  clone this repo

`git clone https://github.com/mandelbrot-project/memo_publication_examples.git`
or via SSH

### 2.  set the correct environmenet

The requested package can be installed by creating the conda environment and activating it.

Use the following line for environment creation 

`conda env create -f environment.yml`

And this one to activate it 

`conda activate memo_examples_env`

If you need to update the environment run 

`conda env update --file environment.yml`

### 3. To use tmap_plotter.py

TMAP package is only avaialble on MacOS and Linux. If you are using Windows, you can use WSL.
Install TMAP package using:

`conda install -c tmap tmap`

### 4. Download the Plant extract dataset data on MASSive repository

Large files are available on [MASSive repository](https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=b753bf1e39cb4875bdf3b786e747bc15#%7B%22table_sort_history%22%3A%22main.collection_dsc%22%2C%22main.collection_input%22%3A%22other%7C%7CEXACT%22%7D):

You need to download:
- the aligned feature table and spectra: see details [here](https://github.com/mandelbrot-project/memo_publication_examples/tree/main/01_input_data/03_plant_extract_dataset/aligned_feat_table_and_spectra)
- the individual .mgf spectra files: see detailes [here](https://github.com/mandelbrot-project/memo_publication_examples/tree/main/01_input_data/03_plant_extract_dataset/individual_mgf_files)

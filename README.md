# CORNETO_GEM_RNA

## Description

This repo provides a simple tutorial with custom functions to run MultiSample iMAT through [CORNETO](https://corneto.org/dev/index.html) with a custom genome-scale metabolic model and RNA-Seq dataset.

The input would be :
 - RNA-Seq normalized matrix of genes X samples (`.csv`)
 - Sample metadata spreadsheet (`.csv`) having 2 columns
   - Sample_ID (Should correspond to column names of RNA-Seq data)
   - Condition (Grouping column)
 - GEM model (`.xml/.sbml`)

The output is a table of reactions as rows and sample fluxes as columns grouped by condition based on metadata.

## Installation

1. Clone the repository
   ```
   git clone https://github.com/Bisho2122/CORNETO_GEM_RNA.git
   cd CORNETO_GEM_RNA
   ```
2. Create and activate conda environment
   ```
   conda env create -f environment.yml
   conda activate corneto_gemsembler
   ```

## Usage

To interactively inspect the results, you can run the [jupyter notebook](https://github.com/Bisho2122/CORNETO_GEM_RNA/blob/main/tutorial.ipynb) tutorial which is an extension of CORNETO's
context-specific-networks [tutorial](https://corneto.org/dev/tutorials/context-specific-metabolic-omics.html).

As another option, run the following commands in your terminal after activating the conda environment to use the example data
```
cd Scripts
```
```
python Run_corneto_imat.py \
-rna ../test_data/ecoli_example_data.csv \
-meta ../test_data/ecoli_example_metadata.csv \
-gem ../test_data/ecoli_core.xml \
-disc quantile \
-o ../test_data
```
To get more information about CLI options :
```
python Run_corneto_imat.py --help
```

## Contributing
1. Fork this repo
2. Create a new branch (`git checkout -b feature/my-feature`)
3. Make your changes
4. Commit and push (`git commit -m "Add new feature"` → `git push origin feature/my-feature`)
5. Open a pull request

## License
This project is licensed under the MIT License. See the [LICENSE](https://github.com/Bisho2122/CORNETO_GEM_RNA/blob/main/LICENSE) file for details.


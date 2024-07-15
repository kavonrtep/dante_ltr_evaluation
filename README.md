# Comparison of DANTE_LTR with Inpactor2 and EDTA

This repository include scripts used for evaluation of LTR-RT annotation method. The scripts were used for comparison of performance of DANTE_LTR, Inpactor2 and EDTA against the standard annotations.


## Requirements

Required packages are listed in `requirements.yaml` file. They can be installed using conda or mamba. To create conda environment with required packages, run the following command:
 
```bash
 conda env create -f requirements.yaml
```

## Data

Data used in evaluation are available in zenodo repository https://zenodo.org/doi/10.5281/zenodo.10891048. The data used in analysis are stored in `reference_genomes` directory which can be downloaded from zenodo repository as zip archive. 


## Included scripts

There are five scripts included in this repository:
- Scripts `compare_annotation_maize.R` and `compare_annotation_rice.R` are used to compare annotation obtained from DANTE_LTR, Inpactor2 and EDTA with reference annotations for maize and rice genomes. 
   - These script are used to convert annotation to standard format suitable for comparison. Each base of the genome is assigned to one of the following four categories: LTR, LTR|Ty1/copia, LTR|Ty3/gypsy, no_annotation. This annoptation is compared with reference annotation for each method and stored in `rice_v7_plots` and `maize_B73_plots` directories as `annot_pairs_*_.csv`.
   - Comparison with reference annotation is then used to calculate FP, FN, TP, TN, sensitivity, specificity, precision, F1 score each method. The results are stored in `stat_LTR_RT.csv` in corresponding directories. These statistices are then used to create radar plots (files `annot_stat_comparison_full_elements_radarchart_ltr.pdf`).
   - Individual elements obtained from each methods are compared with each other to get number of unique and shared elements. The results are experted as Venn diagrams in `annotation_overlaps.pdf` files.   

- `compare_annotation_maize_rm.R`
- `compare_annotation_rice_rm.R`
- `compare_annot_utils.R`
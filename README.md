# Dataset Title
"Comparing DNA Metabarcoding & morphological seed analysis of Neotropical bird feces"

## Author(s)
Carina Motta (carina.i.motta@unesp.br; carinaisabellamotta@gmail.com)

## Description
This dataset contains information on bird fecal samples collected in Corumbataí River Basin, São Paulo, Brazil between August 2022 and August 2023. Seeds found in feces samples (n = 64) were morphologically classified as morphotypes and material from the feces was sequenced using metabarcoding.  Assembled reads are provided here (data > processed_data > metabarcoding) but raw sequence data are available in the NCBI BioProject and in the Sequence Read Archive repository under Accession number PRJNA1203590

## Folder Structure & Files Included
analyses/
├── data/
	- droppings_metabarcoding.csv
	- droppings_seeds.csv
	- droppings_seeds_morphotypes.csv 
├── figures_and_results/
└── scripts/
	- gen_morph_script.R
data/
├── raw_data/ 
	- MottaCarina_bird_captures.xlsx: Capture data and corresponding dropping collection
	- MottaCarina_droppings_seeds.xlsx: Seed IDs based on morphology
	- MottaCarina_droppings_metabarcoding.xlsx: Plant material identified using metabarcoding
├── external_data/
	- newfor_plots_metadata.csv
└── processed_data/
	└── metabarcoding/
		├── OTUs/
			- Assembled OTUs of metabarcoding results; file names correspond to full dataset (data > raw_data > MottaCarina_droppings_metabarcoding.xlsx). Control samples are denoted as "_control".
		└── reads_per_OTU/
			- Number of reads per assembled OTU; file names correspond to full dataset (data > raw_data > MottaCarina_droppings_metabarcoding.xlsx). Control samples are denoted as "_control".


## Methods

Full DNA extraction protocol can be found at dx.doi.org/10.17504/protocols.io.8epv527x5v1b/v1.
Briefly describe:
- Sample collection procedures
- Morphological identification criteria
- DNA extraction, sequencing, and bioinformatics pipeline
- Any filtering or data cleaning steps


## Licensing
This dataset is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)

## Contact
For questions, please contact Carina Motta: carina.i.motta@unesp.br; carinaisabellamotta@gmail.com

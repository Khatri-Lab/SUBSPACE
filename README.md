# SUBSPACE
The SUBSPACE consortium was founded to identify the consensus biologic endotypes in sepsis. Consortium members have shared transcriptomic, proteomic, and clinical data from their respective sepsis biobanks.

The key goals of the SUBSPACE consortium in using this data is to:
1) Validate and test available sepsis endotyping schema
2) Integrate these sepsis endotyping schema to identify commonalities and develop overarching consensus endotypes
3) Evaluate the clinical implications of these endotypes
4) Identify the underlying biology of these endotypes by integrating multi-omic data
5) Develop a simplified consensus endotyping scores to better evaluate the host response to infection

# Datasets:
•	public_genes.csv: the co-normalized genes used for generation of signatures scores
•	public_score_table.csv: the scores created from public_genes. These are what was used for subsequent analyses (clustering and evaluation of outcomes)
•	subspace_genes.csv: the co-normalized genes used for generation of signature scores
•	subspace_score_table.csv: the scores created from subspace_genes. These are what was used for subsequent analyses (clustering and evaluation of outcomes)
•	messi_score_table.csv: the scores generated from MESSI gene expression data
•	glue_score_table.csv: the scores created from co-normalized glue grant data
•	vanish_score_table.csv: the scores created from VANISH public data
•	single_cell_scores_seurat.RDS: integrated Seurat object used for single cell analyses

# Code:
•	subspace_load_data.R: Code sourced to load the necessary data and libraries.
•	tpm_normalization.R: Template used for tpm normalization of count data (this is upstream of datasets provided and is provided for reference only)
•	genescore_calc_template.R: The template used to calculate signature scores. Note for public genesets, did not have genes needed for davenport and cano-gamez SRS scores (this is upstream of score_tables and utilizes genes tables)

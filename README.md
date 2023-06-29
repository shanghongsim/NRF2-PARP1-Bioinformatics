# HU-ATRi-Code-Repository

Heatmap code is in Rcode/Heatmap.Rmd

analysis.py contains all the functions used to analyse data

gene_correlation_screen.ipynb contains script that calculated r value of each siganture/gene against another signature/gene of interest for every cancer in TCGA

single_gene_screen.ipynb contains code to find the correlation of expresssion of every gene in a gene set eg 125 genes in oxidative stress signature to a target gene of interest\

oxstress_factor_boxplot.ipynb plots r values of each signature against RRM2B as a boxplot

calculate_overlap.ipynb compares 125 oxstress gene set aganst 42 antioxidant gene set and calculates the overlap between these two sets.

preprocess_for_gsea.ipynb splits data into low vs high RRM2B cohorts for GSEA analysis

venn_diagram_plotter.ipynb plots Venn diagram with 2-6 sections

gsea_results_plotter.ipynb plots GSEA/STRING results as bar or scatterplot

signature_screen.ipynb finds the r value of every gene signature against RRM2B


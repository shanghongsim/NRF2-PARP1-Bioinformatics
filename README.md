# p53-dependent crosstalk between DNA replication integrity and redox metabolism mediated through a NRF2-PARP1 axis 

> 📣 This paper has been accepted at Nucleic Acids Research!

This is the official repository containing the bioinformatics scripts used in p53-dependent crosstalk between DNA replication integrity and redox metabolism mediated through a NRF2-PARP1 axis (published at Nucleic Acids Research). 🧬

## Figures

| ![Correlation plot](assets/p53%20vs%20rrm2b.png)                     | ![Heatmap](assets/RRM2B%20-%20125%20GO.png)                           | ![Enrichment plot](assets/STRING%20nrf2%20wikipathways.png)            |
|----------------------------------------------------------------------|----------------------------------------------------------------------|----------------------------------------------------------------------|
| Plot showing correlation between normalized log2 RRM2B expression and p53 gene signature activity | Hierarchical heatmap showing the correlation of RRM2B with other genes. | Plot showing the STRING enrichment of the largest cluster derived from hierarchical heatmap |



## Repository structure

```text
Project File Structure:
├── Rcode/                            # Folder containing all bioinformtics scripts in R
│   ├── Heatmap.Rmd                   # Script to plot heatmap figures in the paper
├── assets/                          # Folder containing key figures used for the paper
├── analysis.py                       # Main script for data analysis: fetch, process, analyze and plot
├── calculate_overlap.ipynb           # Script to compute overlap between oxidative stress gene set and antioxidant gene set
├── gene_correlation_screen.ipynb     # Script to compute the pairwise correlation (r-value) of every signature/gene and cancer in TCGA
├── gsea_results_plotter.ipynb        # Script to plot GSEA/STRING results as bar or scatterplot
├── oxstress_factor_boxplot.ipynb     # Script to plot r-values of each signature against RRM2B as a boxplot
├── preprocess_for_gsea.ipynb         # Script to split data into low vs high RRM2B expression cohorts for GSEA analysis
├── signature_screen.ipynb            # Script to find the r-value of every gene signature against RRM2B
├── single_gene_screen.ipynb          # Script to correlate expression of every gene in the oxidative stress signature to a target gene of interest
├── venn_diagram_plotter.ipynb        # Script to plot Venn diagram with 2-6 sets

```


## Citation

Please cite out paper if you found our work useful for your research!

```bibtext
@article{elfar2024p53redoxmetabolism,
      title={p53-dependent crosstalk between DNA replication integrity and redox metabolism mediated through a NRF2-PARP1 axis}, 
      author={Gamal Ahmed Elfar and Obed Aning and Tsz Wai Ngai and Pearlyn Yeo and Joel Wai Kit Chan and Shang Hong Sim and Leonard Goh and Ju Yuan and Cheryl Zi Jin Phua and Joanna Zhen Zhen Yeo},
      journal={Nucleic Acids Research}, 
      volume={52},
      number={20},
      pages={12351--12377},
      year={2024},
      doi={10.1093/nar/gkae811},
      url={https://doi.org/10.1093/nar/gkae811},
}

```

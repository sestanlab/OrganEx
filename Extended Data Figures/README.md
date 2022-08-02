
## Visualize sample size and data quality
plot.quality.R \
#### Heatmap showing cell type heterogeneity
plot.heatm.R

## Hippocampus
#### UMAP visualization & heatmap showing dataset alignment
hip_plot.match.human.R \

#### Augur analysis
- Calculation \
hip_calc.augur.v2.R \

- Visualization \
hip_plot.augur.v2.R \

#### Monocle pseudotime analysis
The following script is also applicable to other organs: \
pseudotimeLiver_monocleSubtype.R

#### Volcano plots showing DEG (all organs)
volcano.plot.R
volcano.plot.ctrl.R - ischemic conditions

#### Gene Ontology analysis (all organs)
plot.GO.ctrl.R

#### Cell death analysis (all organs)
cell.death.enrich.R
celldeath.dotplot.R

#### Cell-Cell interaction analysis
cellchat.hipo.R


## Heart
#### UMAP visualization & heatmap showing dataset alignment
heart_plot.match.human.R \

#### Augur analysis
- Calculation \
heart_calc.augur.R \

- Visualization \
heart_plot.augur.R \

#### Cell-Cell interaction analysis
cellchat.heart.R


## Liver
#### UMAP visualization & heatmap showing dataset alignment
liver_plot.match.public.R \

#### Augur analysis
- Calculation \
liver_calc.augur.R \

- Visualization \
liver_plot.augur.R \

#### Cell-Cell interaction analysis
cellchat.liver.R
The following script is also applicable to other organs: \
cellchat.plot.R 



## Kidney
#### UMAP visualization & heatmap showing dataset alignment
kidney_plot.match.human.R \

#### Augur analysis
- Calculation \
kidney_calc.augur.R \

- Visualization \
kidney_plot.augur.R \

#### Cell-Cell interaction analysis
cellchat.kidney.R
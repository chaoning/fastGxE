---
layout: full
homepage: true
disable_anchors: true
description: A scalable and effective method for genome-wide GxE association analysis 
---
## fastGxE Overview
![iDEA\_pipeline](./images/overview.v1.1.png)
fastGxE is a scalable and effective method designed for genome-wide GxE association analysis. fastGxE handles multiple environmental factors and examines one SNP at a time, decomposing the phenotype into SNP main effect, environmental main effects, GxE interaction effects, while controlling for polygenic effects, polygenic interaction effects, and noise heterogeneity (a). fastGxE evaluates various GxE effect size configurations and combines the resulting *p*-values into a single *p*-value to test whether the SNP interacts with at least one environmental factor (b). By explicitly modeling polygenic background and heteroscedastic noise, fastGxE generates calibrated *p*-values for identifying candidate GxE loci (c). Additionally, it utilizes mmSuSiE, an extension of the SuSiE algorithm, to identify the environmental factors driving the detected GxE interactions and employs the stratified Wald test to visualize and support these interactions (d). fastGxE is implemented as an open-source C++ package, freely available at [Software - Xiang Zhou Lab Website](https://xiangzhou.github.io/software/). 

## User's Guide with fastGxE: [Tutorial](https://chaoning.github.io/fastGxE/documentation/03_Tutorial.html).

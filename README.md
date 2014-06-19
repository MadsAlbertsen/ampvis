ampvis
========

Tools for visualising amplicon sequencing data

## Changelog

### 1.4.0
#### Enhancements
 - amp_ordinate: You can now order the colors using "plot.color.order".
 - amp_ordinate: You can now scale the abundance counts by a sample variable using "scale".
 
#### Bugfixes
 - amp_ordinate: The displayed constrained variance was not displaying correct. 

### 1.3.1
#### Bugfixes
 - amp_rabund: Fixed small error when OTUs was missing Phylum level assignment.
 - amp_heatmap: Fixed small error when OTUs was missing Phylum level assignment.

### 1.3.0
#### Enhancements
 - amp_rabund: You can now order groups using the variable group.order.

#### Bugfixes
 - amp_rabund: Fixed error when trying to display Phylum only information.

### 1.2.1
#### Bugfixes
  - amp_heatmap: didn't handle datasets without genus level classification proberly.
  - amp_rabund: didn't handle datasets without genus level classification proberly.

### 1.2.0
#### Enhancements
 - amp_heatmap: Now supports transformation ("log" or "sqrt") of the background color using the variable plot.colorscale.

### 1.1.0
#### Enhancements
 - General: 10.000 sequences is now the default vaule of scale.seq.
 - amp_heatmap: Supports additional taxonomic labels using the tax.add variable.
 - amp_heatmap: Default color scaling is now square root instead of log10 to better display 0 counts.
 - amp_heatmap: If no taxonomic level is present at the tax.aggregate level the function now reports the best classification along with the OTU name (used to be just the OTU name).
 - amp_rabund: If no taxonomic level is present at the tax.aggregate level the function now reports the best classification along with the OTU name (used to be just the OTU name).
 -amp_ordinate: Supports labeling the plotted OTUs (plot.nspecies) with any taxonomic level through the plot.nspecies.tax variable.

### 1.0.0
First release of the ampvis package
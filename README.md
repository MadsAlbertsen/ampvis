ampvis
========

Tools for visualising amplicon sequencing data.

[Start here](http://madsalbertsen.github.io/ampvis)

## Changelog

### 1.9.2
#### Bugfixes
  - amp_rabund: ggplot2 removes 0 values when using log10 scales - as you cant take log10(0). This means that the boxplots are not showing the correct median of the complete dataset. The parameter `adjust.zero` can be used to add a small constant to all OTU abundances that are 0 in order to display the correct median values when using log scales.

### 1.9.1
#### Enhancements
  - amp_ordinate: `plot.shape` can now be used to assign different shapes to points based on a sample variable.
  - amp_ordinate: Now includes a `plot.theme` option that can be used to make publication friendly images fast.
  - amp_ordinate: The option `envfit.show` can be used to disable showing envfit vectors on the plot.
  - amp_test_speices: Now includes a `plot.theme` option that can be used to make publication friendly images fast.
  - amp_heatmap: Now includes a `plot.theme` option that can be used to make publication friendly images fast.
  - amp_test_cluster: Now includes a `plot.theme` option that can be used to make publication friendly images fast.

### 1.9.0
#### Enhancements
  - data: The newest version of the MiDAS data have been included as `data(MiDAS_1.20)`. It includes small updates to the taxonomy and a single sample have been removed (2591) due to wierdness. The OTU names are now locked! Hence, OTUs in this release will have the same name in the next release.
  - data: Note that the MiDAS_1.20 data is not rarefied by default anymore.
  - data: The DNA extraction data is now included as `data(DNAext_1.0)`.
  - amp_test_species: Updated to work with the same input as the rest of the functions.
  - amp_test_cluster: Updated to work with the same input as the rest of the functions.
  
#### Bugfixes
  - amp_load: the rarefy parameter only worked with 10000 seqences.
  - amp_heatmap: order.x and order.y was broken. fixed.

### 1.8.3
#### Bugfixes
  - amp_ordinate: The envfit.factor function was broken. Fixed.

### 1.8.2
#### Bugfixes
  - amp_rename: Now can convert more than 1 phylum to class level.

### 1.8.1
#### Bugfixes
  - amp_rabund: Now sorts the boxplots by median as default.

#### Enhancements
  - amp_rabund: THe new option "sort.by" can be used to sort the boxplots by Median, Mean or Total. 

### 1.8.0
#### Enhancements 
 - Major speed improvement to generation of data.frames from phyloseq objects. amp_convert is not needed anymore. Updated in all functions.
 - amp_core: Now calculates absolute frequencies instead of relative.

#### New functions
 - amp_export: export reference sequences contained in a phyloseq object to a file.

#### Bugfixes
 - amp_core: A bug was present in the "frequency" plot that made rare species more abundant than they were supposed to be! If you've used the function, then please recalculate using the updated function.

### 1.7.1
#### Bugfixes
 - Plyr and dplyr loading order messed up the functions. Fixed.

### 1.7.0
#### Enhancements
 - data.table and dplyr have been implemented to replace plyr and phyloseq functions. This results in a massive speed increase.
 - All functions now accept a list of data.frames instead of a phyloseq object. This significantly improves speed.

#### New functions
 - amp_convert: Convert a phyloseq object to a list of dataframes.
 - amp_load: Load data and convert it to a phyloseq object.

### 1.5.0
#### Enhancements
 - amp_rabund: You can now flip the axis using plot.flip = T.

### 1.4.1
#### Bugfixes
 - amp_core: Fixed small bug in the reporting of the mean abundance in output data.

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
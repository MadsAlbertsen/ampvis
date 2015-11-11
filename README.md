ampvis
========

Tools for visualising amplicon sequencing data.

[Start here](http://madsalbertsen.github.io/ampvis)

To cite ampvis please use:

Albertsen M, Karst SM, Ziegler AS, Kirkegaard RH, Nielsen PH (2015) Back to Basics – The Influence of DNA Extraction and Primer Choice on Phylogenetic Analysis of Activated Sludge Communities. PLoS ONE 10(7): e0132783. doi:10.1371/journal.pone.0132783

## Changelog
### 1.21.0
#### Enhancements
 - New function: `amp_venn` generate venn diagrams of core species.
 
#### Bugfixes
 - amp_load: Fixed problems with taxonomies that had leading and trailing spaces (@Kirk3gaard).

### 1.20.0
#### Enhancements
 - amp_function: Updated MiDAS functional data to the current version. It's now possible to use all the functional data in MiDAS.
 - amp_cleanMiF: Small function for internal use to format raw MiDAS functional data to ampvis format.

### 1.19.0
#### Enhancements
 - amp_rarecurve: Now supports a basic legend (@Kirk3gaard).
 - amp_heatmap: Now supports the use of other color scales directly using the `color.vector` option (@Kirk3gaard).
 
#### Bugfixes
 - amp_rarecurve: Colors were messed up. Fixed now.

### 1.18.0
#### Enhancements
 - amp_rarecurve: Now supports coloring by groups (@Kirk3gaard).
 - amp_rabund: Now supports manual sorting of the y-axis by `order.y` (@Kirk3gaard).
 - amp_heatmap: The new option `sort.by` enables to show the x most abundant taxa from a specific group. Default option is to show the x most abundant taxa on average.

### 1.17.1
#### Bugfixes
 - amp_heatmap: Fixed an error related to the new `min.abundance` option. The default value is now 0.1 which fixes the problems with log10 and 0 values.
 - amp_rabund: `scale.seq` is now 100 by default. Like all other functions.


### 1.17.0
#### Enhancements
 - New function: `amp_stats` can be used to generate a table with metadata, number of sequences and alpha-diversity indices.
 - New function: `amp_rarecurve` can be used to generate rarefaction curves.

### 1.16.1
#### Bugfixes
 - Added the `scales` library as a dependency.
 - amp_heatmap: Various bugfixes in the clustering option. It should now work with all other settings.

### 1.16.0
#### Bugfixes
 - amp_rabund: `tax.show` now correctly shows that 50 are displayed by default (@Kirk3gaard).
#### Enhancements
 - amp_heatmap: Added the ability to calculate and display median values in the `calc` parameter (@Kirk3gaard).

### 1.15.0
#### Enhancements
 - amp_heatmap: `order.y` now supports the input "cluster". Which will do basic hclust on the presented data and arrange it according to the clustering.
 - amp_heatmap: The new options `min.abundance` and `max.abundance` can be used to limit the values displayed in the plot. This also impacts the new clustering option in `order.y`.
 - amp_heatmap: `scale.seq` is now 100 by default. Hence, assuming that the data have been converted to a percentage.

### 1.14.0
#### Enhancements
 - amp_ordinate: `plot.group.manual` can now be used to group samples independently on the `plot.color` option.

### 1.13.2
#### Bugfixes
 - Fixed integer vs. numeric error when not transforming data.

### 1.13.1
#### Bugfixes
 - Fixed wrong link in the citation information.

### 1.13.0
#### Enhancements
 - amp_export_table: New function to export an otutable from a phyloseq object.
 - amp_rabund: Added plot.theme functionallity.
 - added citation information.

### 1.12.0
#### Enhancements
  - amp_stability: `group` now supports a vector of parameters.
  - amp_stability: The individual parametes in group can now be used in `facet_wrap`.
  - amp_test_species: Enabled the parallel option in DESeq2 (thanks to: @mdehollander).
  - amp_heatmap: The individual parametes in group can now be used in `facet_wrap`.

#### Bug fixes
  - amp_stability: Now works as intended with replicates.

### 1.11.0
#### Enhancements
  - amp_stability: Calculate and plot overall similarity between samples in timeseries.

### 1.10.0
#### Enhancements
  - amp_function: Make a function table sorted by an input heatmap.
  - data: The current version of the online MiDAS function data is now included (2015-04-24). It's used by `amp_function`, but can also be loaded manually using `data(MiF)`.

### 1.9.2
#### Enhancements
  - amp_heatmap: Now allows calculation of max values using the `calc` option. Default is still mean.
  - amp_export: Now supports export of the taxonomic string in the header using the `tax` option. 

#### Bugfixes
  - amp_rabund: ggplot2 removes 0 values when using log10 scales - as you cant take log10(0).. This means that the boxplots are not showing the correct median of the complete dataset. The parameter `adjust.zero` can be used to add a small constant to all OTU abundances that are 0 in order to display the correct median values when using log scales.

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
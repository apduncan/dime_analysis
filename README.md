# *D*ietary B*I*oactives and *M*icrobiome Div*E*rsity (DIME) Analysis

Dietary bioactives have been associated with positive health effects such as 
cardiovascular health, or having anti-diabetic properties. Many bioactive 
compounds survive digestion arriving in the colon, where the gut microbiome can 
further metabolise them to forms more available to the host. Previous work has 
shown diets such as the Mediterranean diet, enriched in plant bioactives, 
affecting the composition and diversity of the gut microbiome. However, it is 
unclear which factors in such whole diet interventions underlie these changes, 
such as the bioactives themselves, the percentage of dietary fibres or 
macronutrient composition. In the *D*ietary B*I*oactives and *M*icrobiome 
*D*ivErsity (DIME) study we investigated the impact of a diet rich in a wide 
range of plant bioactives compared to a low bioactive diet with matched total 
fibre intake.

This repository contains plotting and analysis scripts supporting the
manuscript analysing the DIME study.

## Data
All neccesary data is included in `data.tar.gz`, available from the associated
figshare repository due to size: [10.6084/m9.figshare.29860841](https://doi.org/10.6084/m9.figshare.29860841).
Before running any scripts, download and extract this archive: `tar -xvf data.tar.gz`

## Running
### Restore environment
Dependencies are managed using 
[renv](https://rstudio.github.io/renv/articles/renv.html).
To install the required packages, in an `R` session, run `renv::restore()`.
`R` version 4.3.3 was used for the analysis.

### targets
The analysis is written as a [targets](https://books.ropensci.org/targets/)
pipeline. To produce any of the results, you can use the `targets` functions
```
> tar_make(fig_three)
✔ skipping targets (1 so far)...
▶ dispatched target tbl_metab_additional_annotations
● completed target tbl_metab_additional_annotations [0.148 seconds, 3.44 kilobytes]
✔ skipping targets (9 so far)...
▶ dispatched target tbl_metab_peaks
● completed target tbl_metab_peaks [20.852 seconds, 211.887 megabytes]
✔ skipping targets (12 so far)...
▶ dispatched target plt_metab_volcano_revised_ident
Joining with `by = join_by(feature_id)`
● completed target plt_metab_volcano_revised_ident [0.41 seconds, 428.539 megabytes]
✔ skipping targets (16 so far)...
▶ dispatched target fig_three
● completed target fig_three [0.095 seconds, 430.788 megabytes]
▶ ended pipeline [3.488 minutes]
> tar_load(fig_three)
> fig_three
```
Output will typically be produced in `/outputs/{figures|tables|models}`.

#### Network targets
Code is provided to generate networks in `R/network_construct.R`.
However, the networks we generated and used for our results are distributed
in the figshare repository, and are used for all figure generation. 

To recreate networks, run
```
tar_make(vct_net_high_paths)
tar_make(vct_net_low_paths)
```

This will make new networks in `output/network_rerun`, producing three `Rds`
files for each network, and igraph format network, and association matrix,
and the SPIEC-EASI output.
These could be moved to `data/network` to use in plots, though be aware
the netowrk plotting functions are quite specific to the bundled structures,
so if for some reason you learn a quite different structure it may not
be useful.

Networks are built from a subset of PFAM and metabolites.
Metabolites are those identified as significant using t-tests, and similarly
we identify PFAMs siginificantly associated to diet using non-parametric
tests, implemented in `network_prevalence_filter` in `R/network_construct.R`.
To change the filtering, you can change the implementation of this function.

To change parameters passed to SPIEC-EASI, edit the function `se_params` in
`R/network_construct.R`.

#### Run all
To run the full pipeline, which will produce output figures and tables, run
```
library(targets)
tar_make()
```
Individual results objects can then be accessed with `tar_load()`, and 
tables and figures found in `output`.
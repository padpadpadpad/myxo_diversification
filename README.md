
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Looking at adaptive diversification of Myxobacteria in natural environments

## Summary

This repository is a shared resource for us to play about with data and
share methods in relation to the project looking at adaptive
diversification of Myxobacteria in natural environments.

## Where we are at

  - [x] Sampled different habitats and sites across Cornwall **(Autumn
    2020)**
  - [x] DNA extracted from these samples
  - [x] Designed primers for amplifying Myxobacteria **(January
    2021-June 2021)**
      - [ ] Scripts uploaded
      - [ ] Bought primers uploaded
  - [x] Samples 16S sequenced **(May 2021)**
      - [x] Processed data uploaded
      - [ ] Analysis scripts uploaded
      - [ ] Clustering analysis to look at habitat separation and
        potential definitions
      - [ ] Data archived and backed up
  - [ ] Samples amplicon sequenced to amplify rpoB gene **(October
    2021)**
      - [ ] Processed data uploaded
      - [ ] Exploratory analysis scripts

## To use

The GitHub repo should be super easy to use on everyone’s machine using
RStudio projects. RStudio projects automatically assigns the root
directory to the directory in which the project resides. Consequently
all of the analyses should be runnable without altering paths. Simply
double-click on the `.Rproj` folder that is in the base folder of this
repository and it should open the project in RStudio.

I have avoided uploading overly large files to this repository for
obvious reasons. The sequencing data is stored as
[phyloseq](https://joey711.github.io/phyloseq/) objects and should be
easy to query.

## Folder structure

### data

**data** includes the processed data from each part of the project.
`sequencing_16s` contains the data for the 16s sequencing etc.

  - **metadata.csv** - completed metadata for each sample
  - **sequencing 16s** - folder containing data from the 16s sequencing

<!-- end list -->

1.  `ps_16s_complete.rds` - the original phyloseq object from the
    processing pipeline
2.  `ps_16s_prev_filtered.rds` - low abundance and low prevalence ASVs
    removed
3.  `ps_16s_low_depth_removed.rds` - low depth samples removed **no
    samples were too low here**
4.  `ps_16s_rarefied.rds` - rarefied samples

### scripts

**scripts** includes the scripts from each part of the project so far.
`sequencing_16s` contains the scripts for the 16s sequencing etc. At the
top of each script there is a summary of what has been done in each
script. Order of scripts to be run is the same as numbered here.

  - **sequencing\_16s**

<!-- end list -->

1.  `prevalence_filtering.R`
2.  `rarefy_data.R`
3.  `first_clustering.R`
4.  `second_clustering.R`

<!-- end list -->

  - **sequencing\_rpoB**

### plots

**plots** includes the plots from each part of the project so far.
`sequencing_16s` contains the plots for the 16s sequencing etc.

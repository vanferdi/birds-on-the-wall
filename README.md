# Overview

This repository contains all of the code and figures found in the paper:

[put final citation here]
Ferdinand, Pattenden, Brightsmith, &amp; Hobson (2022)


# Code by section in paper:

## Summarizing the extent to which species mix

`solo_presence.R` calculates column 5 of Figure 3.

`[Liz please upload]` contains everything else in Figure 3.

## Copresence patterns

`copresence_counts.R` creates the heatmap of copresence counts (Figure 4) and the community detection analysis on these counts.

`copresence_network.R` contains the glmer analysis of significant ties and creates the network in the top two panels of Figure 5a & 5b. The statistical results of the glmers are saved in `copresence_glmers_output.txt`. It also runs the community detection analysis and centrality analyses on the ties in 5a.

`body_size.R` creates Figure 6.

## Dynamic joining decisions

`joining_network.R` creates Figure 5c & 5d and runs the community detection analysis on the ties in 5c. The statistical results of the glmers are saved in `joining_glmers_output.txt`.

## Inferring categorization schemas from joining decisions

`coarsegraining_analysis.R` finds the best-fit categorization schema per species using model selection and AIC, reported in Figure 7. The statistical results of the glmers are saved in `coarsegraining_glmers_output.txt`

`complexity_and_centrality.R` creates Figure 8.


## Supplemental: Overview of clay lick use

Figure S1 is created in ` `.

Figure S2 is created in ` `.

Figure S3 is created in ` `.

The data in Table S1 can be found in `coarsegraining_glmers_output.txt`.

All of the information in Figure S4 can be inferred directly from Figure 7.

# Overview

This is a copy of the main repo where I'm saving the results of the validations.

I erased everything in the directories Plots and Rdata (except `circle_layout.Rdata`). <br>
Then I ran everything below from scratch using the validation data set called odddays_dataset.csv.

Below is a list of the dependencies you need to run to do each validation.

## 1) Body size validation

1) run `copresence_network.R` to get `copresence_network_ties.Rdata`

2) run `body_size.R`


## 2) Complexity vs centrality validation

1) run `copresence_network.R` to get `copresence_network_ties.Rdata`

2) run `adjacentizer.R` to get `adjacentized_df.RData`

3) run `joining_network.R` to get `joining_columns.Rdata`

4) run `coarsegraining_analysis.R` to get `ncats.RData`

5) run `complexity_and_centrality.R` to get `complexity_centrality.pdf`
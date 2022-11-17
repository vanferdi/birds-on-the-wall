##########################################################################################
# overview
##########################################################################################

# count the number of pairwise co-occurrences between species in the dataset
# and run a community detection algorithm on them

# input:
# evendays_dataset.csv

# output:
# copresence_counts_heatmap.pdf


##########################################################################################
# setup
##########################################################################################

# load common functions
source("./Code/function_library.R")

options(scipen=999)  # run this once to disable scientific notation in print to screen

##########################################################################################
# load data
##########################################################################################

db <- read.csv("./Data/evendays_dataset.csv", header = TRUE, sep = ",")
nrow(db)  # 7362 rows

# remove data from unstandardized spots
db <- filter_spots(db,spots)
nrow(db)  # 7357, correct

# change columns to the shorter species codes
colnames(db)[4:15] <- ss  # change cols to shorter species codes

# convert data to binary presence/absence of species
dbb <- binarize_df(db)

# should look like this:
head(dbb,3)
"       Rdate Xhr.dectime Spot BH BY CF DH ME OC RB RG SC WB WE YC
1 2002-10-02    5.500000   3B  0  1  1  0  0  0  1  0  0  0  0  1
2 2002-10-02    5.583333   3A  1  0  1  0  0  0  1  0  0  0  0  1
3 2002-10-02    5.583333   3B  0  1  1  0  0  0  0  0  1  0  0  1"


##########################################################################################
# get the raw number of times each pair of species co-occurred on the same spot
##########################################################################################

# get the number of times each pair of species was present together on the smae spot
get_copresence_counts <- function(df) {  # df must be binary
	M <- matrix(0,12,12)
	colnames(M) <- ss
	rownames(M) <- ss
	
	# got through every row in the dataframe
	for (r in 1:nrow(df)) {
		print(r)
		slice <- df[r,][4:15]
		species_present <- which(slice==TRUE)
		
		# if more than one species was present
		if (length(species_present) > 1) {
			# get all combinations of species pairs there
			pair_combos <- combn(species_present,2)

			# go through each combination present and save one tally for it to the matrix
			for (i in 1:ncol(pair_combos)) {
				one <- pair_combos[,i][1]
				two <- pair_combos[,i][2]
				M[one,two] <- M[one,two]+1  # this only fills in the top triangle of M
			}
		}
	}
	return(M)
}

M_counts <- get_copresence_counts(dbb)

M_counts
"  BH  BY   CF  DH   ME   OC   RB  RG   SC  WB  WE  YC
BH  0 947  953  58 1258 1071  879 143  589 111 387 334
BY  0   0 1015  41 1152  753  943 386 1101  77 325 368
CF  0   0    0 108 1267  705 1343 118  620 164 423 492
DH  0   0    0   0   70  279   80   7   26 379 341  11
ME  0   0    0   0    0 1019 1063 189  767  92 410 303
OC  0   0    0   0    0    0  635 122  464 294 468 235
RB  0   0    0   0    0    0    0 118  563 129 410 413
RG  0   0    0   0    0    0    0   0  355  11  53  44
SC  0   0    0   0    0    0    0   0    0  31 207 200
WB  0   0    0   0    0    0    0   0    0   0 607  46
WE  0   0    0   0    0    0    0   0    0   0   0 109
YC  0   0    0   0    0    0    0   0    0   0   0   0"

# make it into a symmetrical matrix for plotting to Figure
M_full <- t(M_counts)+M_counts

# set diag to NA so it will plot in grey (as opposed to count = 0)
diag(M_full) <- NA

###################################################################################
# create heatmap
###################################################################################

require(ggplot2)

row <- rep(ss,12)
col <- rep(ss, each=12)
count <- as.numeric(M_full)
data <- data.frame(row,col,count)

# change the order of species on axes to their body size order
data$row <- factor(data$row, levels=c("RG","BY","SC","ME","YC","CF","RB","BH","OC","WE","WB","DH"))
data$col <- factor(data$col, levels=c("RG","BY","SC","ME","YC","CF","RB","BH","OC","WE","WB","DH"))

co <- "#333333"  # set border and axis text color to dark grey

pdf(file = "./Plots/copresence_counts_heatmap.pdf", width = 7, height = 6)
par(mfrow=c(1,1), mar=c(1,1,1,1))
ggplot(data, aes(row, col, fill=count)) + 
	geom_tile() + 
	xlab("") + ylab("") +
	geom_text(aes(label = count), col="black") +
	scale_fill_gradient(low = "white", high = myblue) +
	scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + # remove padding
	theme(axis.text = element_text(size = 15, color=co), 
		  panel.border = element_rect(color = co, fill=NA, size=1)) # border around whole heatmap
dev.off()


###################################################################################
# community detection on copresence counts
###################################################################################

require(igraph)

# get all ties and their weights (from the full version of the matrix)
all_ties <- edgelist_from_matrix(matrix=M_full)
tie_weights <- weightlist_from_matrix(M_full,all_ties)

all_ties
tie_weights

# make graph object
g <- make_empty_graph(n = 12, directed=F) + edges(all_ties, weight=tie_weights)

# run the community detection algorithm
result <- leading.eigenvector.community(g, weights=E(g)$weight)
result
"IGRAPH clustering leading eigenvector, groups: 2, mod: 0.06
+ groups:
$`1`
[1]  1  2  3  5  7  8  9 12

$`2`
[1]  4  6 10 11"

# the two groups are:
ss[result[1]$`1`]  # "BH" "BY" "CF" "ME" "RB" "RG" "SC" "YC"
ss[result[2]$`2`]  # "DH" "OC" "WB" "WE"

# That's what I expected, but I wasn't sure where OC would go because it looks like a broker

# leading.eigenvector.community() implements this:
# MEJ Newman: Finding community structure using the eigenvectors of matrices, Phys Rev E 74:036104 (2006).
# https://www.rdocumentation.org/packages/igraph/versions/0.4.3/topics/leading.eigenvector.community

# leading.eigenvector.community() converts everything to undirected ties 
# works with weighted and non-weighted networks (must be non-negative weights)


###################################################################################
# END
###################################################################################


###################################################################################
# Not in paper: plot heatmap showing expected copresence by chance

# I want to make sure there isn't anything in the baseline presence of each species
# that is creating any community structure by itself.

# compare the above heatmap to one showing chance copresence
# based on the baseline probability of observing each species in any given scan
# (a scan = one row of the dataset)

baseline <- get_baseline(dbb)
# species_n: shows total number of scans (i.e. row of dataframe) in which the species appeared
# species_p: gives the probability of observing the species in a scan

# get the probability of each species occurring together by chance
# if we merely randomized all the entries within each species column
Mb <- baseline$species_p %o% baseline$species_p  # get the outer product

# calculate expected number of co-occurrences from the probabilities
Me <- matrix(0,12,12)
N <- nrow(dbb)                       # total number of scans
for (i in 1:12) {  
	for (j in 1:12) {
		pA <- baseline$species_p[i] # probability of species i appearing in a scan
		pB <- baseline$species_p[j] # probability of species j appearing in a scan
		theta <- pA*pB              # independent probability of the two appearing together by chance
		Me[i,j] <- round(theta*N)   # shows theta as expected number of events between the pair by chance
	}
}
Me

# create dataframe from matrix
diag(Me) <- NA
row <- rep(ss,12)
col <- rep(ss, each=12)
count <- as.numeric(Me)
data <- data.frame(row,col,count)
# change the order of species on axes to their body size order
data$row <- factor(data$row, levels=c("RG","BY","SC","ME","YC","CF","RB","BH","OC","WE","WB","DH"))
data$col <- factor(data$col, levels=c("RG","BY","SC","ME","YC","CF","RB","BH","OC","WE","WB","DH"))

#pdf(file = "./Plots/copresence_counts_heatmap_bychance.pdf", width = 7, height = 6)
#par(mfrow=c(1,1), mar=c(1,1,1,1))
ggplot(data, aes(row, col, fill=count)) + 
	geom_tile() + 
	xlab("") + ylab("") +
	geom_text(aes(label = count), col="black") +
	scale_fill_gradient(low = "white", high = myblue) +
	scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + # remove padding
	theme(axis.text = element_text(size = 15, color=co), 
		  panel.border = element_rect(color = co, fill=NA, size=1)) # border around whole heatmap
#dev.off()

# RG, YC, DH have such light rows because they're the 3 rarest visitors on the clay lick.
baseline[order(baseline$species_p),]
"     s species_n  species_p
8  RGMA       491 0.06673916
12 YCPA       738 0.10031263
4  DHPA       861 0.11703140
9  SCMA      1464 0.19899416
6  OCPA      1740 0.23650945
1  BHPA      1776 0.24140275
11 WEPA      1903 0.25866522
10 WBPA      2030 0.27592769
7  RBMA      2111 0.28693761
2  BYMA      2185 0.29699606
3  CFMA      2333 0.31711295
5  MEPA      2376 0.32295773"

# Do community detection on the chance encounters 
# sanity check: it should be one community

# get all ties and their weights from Me
all_ties <- edgelist_from_matrix(Me)
tie_weights <- weightlist_from_matrix(Me,all_ties)

# make graph object
g <- make_empty_graph(n = 12, directed=F) + edges(all_ties, weight=tie_weights)

# run the community detection algorithm
result <- leading.eigenvector.community(g, weights=E(g)$weight)
result
"groups:
$`1`
[1]  1  2  3  4  5  6  7  8  9 10 11 12"

# yep it is one community
# this means there's nothing in the baseline presence of species that would create community structure




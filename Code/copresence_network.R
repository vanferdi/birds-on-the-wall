##########################################################################################
# overview
##########################################################################################

# create a network graph of significant copresence ties between species

# evendays_dataset.csv data is processed into two matrices: copr_est, copr_pvals
# which are saved as Rdata: copresence_glmer_matrices.RData
# the network is created from copr_est, copr_pvals

# significant ties are computed from glmers run pairwise between all species
# copr_est contains the estimate of the effect
# copr_pvals contains the pvalue of the comparison
# pvals are corrected for multiple comparisons and only the ones that remain significant are plotted
# results from the glmer are logged in file copresence_glmers_output.txt

# you can skip to section "load results from Rdata" to load these matrices and go from there

# input:
# evendays_dataset.csv
# circle_layout.RData

# output:
# copresence_glmers_output.txt
# copresence_glmer_matrices.RData
# copresence_network_ties.RData
# copresence_network_positive.pdf
# copresence_network_negative.pdf


##########################################################################################
# setup
##########################################################################################

# load common functions
source("./Code/function_library.R")

require(lme4)  # package for running mixed effects regression analyses

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


###################################################################################
# run glmers
###################################################################################

copr_pvals <- matrix(0,12,12)  # save pvals here
colnames(copr_pvals) <- ss
rownames(copr_pvals) <- ss

copr_est <- matrix(0,12,12)    # save estimates here
colnames(copr_est) <- ss
rownames(copr_est) <- ss

# print results to file - remember to run sink() to end this process and print to console again
sink(file = "./Code/Output/copresence_glmers_output.txt")

for (i in 1:12) {
	focal <- ss[i]  # get current focal species
	print(noquote("==========================="))
	print(noquote(paste("CURRENT FOCAL SPECIES IS",focal)))
	
	# create the 12 full formulas, one for each focal species as the dependent variable
	# use the 11 non-focal species as the predictor species, plus random effects for spot, time and date
	predictor_species <- paste0(ss[-i], sep="+", collapse="")  # remove the focal species from predictors
	f <- formula(paste(focal,"~", predictor_species,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	full <- glmer(f, family="binomial", nAGQ=0, data=dbb)
	fs <- summary(full)
	print(fs)
	
	# given the current focal species, create the 11 reduced models
	for (j in 1:11) {
		x <- ss[-i]  # remove the focal
		one_to_remove <- x[j]  # and remove the predictor species being reduced
		x <- x[-j]
		predictor_species2 <- paste0(x, sep="+", collapse="")
		r <- formula(paste(focal,"~", predictor_species2,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
		
		print(noquote("==========================="))
		print(noquote(paste("=== removing species",one_to_remove,"===")))
		
		# conduct anova on the full and reduced model
		reduced <- glmer(r, family="binomial", nAGQ=0, data=dbb)
		a <- anova(full,reduced)
		print(a)  # will print to sink file instead of screen
		
		# get pvalue and save to copr_pvals
		pval <- a$`Pr(>Chisq)`[2]
		copr_pvals[focal,one_to_remove] <- round(pval,10)  # focal is on row, predictor species is on col
		
		# get estimate of effect and save to copr_est
		value <- fs$coefficients[one_to_remove,][1]
		copr_est[focal,one_to_remove] <- value
		print(noquote(paste("Estimate:",value)))  # print the estimate for the current species off the full model
	}
	if (i == 12 & j == 11) {  # when everything has finished running
		sink()  # close the sink connection and create the file
		print("DONE")
	}
}

###################################################################################
# save results as Rdata
###################################################################################

# save the two resulting matrices as Rdata
save(copr_est, copr_pvals, file = "./Rdata/copresence_glmer_matrices.RData")


###################################################################################
# load results from Rdata
###################################################################################

load("./Rdata/copresence_glmer_matrices.RData")

round(copr_est,3)
round(copr_pvals,4)


###################################################################################
# extract significant ties
###################################################################################

# gets the significant ties with FDR correction for multiple comparisons

# save ties as a matrix where:
# 1 = significant affiliative tie
# -1 = significant avoidant tie
# 0 = tie not significant

ties <- get_sig_ties(copr_pvals, copr_est)
"  BH BY CF DH ME OC RB RG SC WB WE YC
BH  0  1  1 -1  1  1  1 -1  0 -1  1  1
BY  1  0  1 -1  1  1  0  1  1 -1 -1  1
CF  1  1  0  0  1  0  1 -1  0 -1  0  1
DH -1 -1  0  0 -1  1 -1 -1 -1 -1  0 -1
ME  1  1  1 -1  0  1  1  0  1 -1  0 -1
OC  1  1  0  1  1  0  0  0  0  0  1  0
RB  1  0  1  0  1  0  0 -1 -1 -1  0  0
RG  0  1 -1 -1  0  0 -1  0  1 -1 -1  0
SC  0  1  0 -1  1  0  0  1  0 -1 -1  0
WB  0 -1 -1 -1 -1  0 -1 -1 -1  0 -1 -1
WE  1 -1  0  0  0  1  0 -1 -1 -1  0  0
YC  1  1  1 -1 -1  0  1  0  0 -1 -1  0"

save(ties, file = "./Rdata/copresence_network_ties.RData")

# explanation of resulting matrix:
# the glmers checked if the focal species presence was predicted by each of the other species' presence
# each row of the matrix is one focal species and each column is a predictor species
# all non-zero entries in "ties" denote significant relationships between species
# and the sign (plus or minus) denotes the direction of the relationship
# these relationships are directed, meaning BH ~ RG doesn't necessarily equal RG ~ BH

# for example,
# the presence of BH is negatively predicted by RG: ties["BH","RG"] = -1
# check "copresence_glmers_output.txt" and you will see the estimate for BH ~ RG is -0.43690 and the pval is 0.00564
# also, the presence of RG is not significantly prediced by BH: ties["RG","BH"] = 0
# when RG is the focal species, the estimate for RG ~ BH is copr_est["RG","BH"] = -0.25461 and the pval is 0.1391


###################################################################################
# show ties as a binary heatmap and check reciprocity
###################################################################################

require(ggplot2)

ties_heatmap <- function(M_ties) {
	row <- rep(ss,12)
	col <- rep(ss, each=12)
	val <- as.numeric(M_ties)
	dir <- sign(val)  # convert tie value to its direction only (if not already)
	
	color <- rep("#FFFFFF",length(dir))  # white
	for (i in 1:length(dir)) {
		if (dir[i] > 0) { color[i] <- "#3ec9f7" }  # blue
		if (dir[i] < 0) { color[i] <- "#cf1c08" }  # red
	}
	
	data <- data.frame(row,col,val,dir,color)
	
	ggplot(data, aes(col, row, col=color, fill= color)) + 
		geom_tile() + 
		xlab("predicted by") + ylab("presence") +
		scale_color_identity() +
		scale_fill_identity()
}

ties_heatmap(ties)

# check reciprocity of copresence ties
144-sum(ties==t(ties))  # there are 12 unreciprocated ties

# they are:
# BH->RG  negative
# RG->BH  n.s.

# BH->WB  negative
# WB->BH  n.s.

# DH->RB  negative
# RB->DH  n.s.

# RB->SC  negative
# SC->RB  n.s.

# RB->YC  n.s.
# YC->RB  postive

# WE->YC  n.s.
# YC->WE  negative

sum(ties==1)  # 45 positive ties
sum(ties==-1) # 47 negative ties

sum(ties+t(ties)==2)   # 44 reciprocated positive ties
sum(ties+t(ties)==-2)  # 42 reciprocated negative ties

sum(ties+t(ties)==-1)  # 10/2 = 5 of the unreciprocated ties are n.s. and negative
sum(ties+t(ties)==1)   # 2/2 = 1 of the unreciprocated ties are n.s. and positive


###################################################################################
# plot network
###################################################################################

require(igraph)

# load the main layout
load("./Rdata/circle_layout.RData")  # loads circle_layout, 
circle_layout  # ith row corresponds to ith species in ss (row 1 is BH)

# get your edgelist from a matrix - takes matrix with pos and neg ties
edgelist_from_matrix <- function(matrix,tie_type="all") {
	if (tie_type == "all") { t <- which(matrix != 0,arr.ind=TRUE) }
	if (tie_type == "positive") { t <- which(matrix > 0,arr.ind=TRUE) }
	if (tie_type == "negative") { t <- which(matrix < 0,arr.ind=TRUE) }
	return(as.numeric(t(t)))
}

plotit <- function(filename,M_ties,tie_type="both") {
	pos_ties <- edgelist_from_matrix(matrix=M_ties,tie_type="positive")
	neg_ties <- edgelist_from_matrix(matrix=M_ties,tie_type="negative")
	
	if (tie_type == "both") {
		g <- make_empty_graph(n = 12, directed=T) +
			edges(pos_ties, color="#3ec9f7", weight=3) +
			edges(neg_ties, color="#cf1c08", weight=1)
	}
	if (tie_type == "positive") {
		g <- make_empty_graph(n = 12, directed=T) +
			edges(pos_ties, color="#3ec9f7", weight=2)
	}
	if (tie_type == "negative") {
		g <- make_empty_graph(n = 12, directed=T) + 
			edges(neg_ties, color="#cf1c08", weight=2)
	}
	
	pdf(file = filename, width = 8, height = 6)
	par(mfrow=c(1,1), mar=c(1,1,1,1))
	plot(g, layout=circle_layout, rescale=T, asp=0, edge.arrow.size=1, 
		 vertex.label.cex=1.2, vertex.label.family="Helvetica", vertex.label.font=2,
		 vertex.label=ss, vertex.shape="circle", vertex.size=12, 
		 vertex.color="grey", vertex.label.color="black", edge.width=E(g)$weight)
	dev.off()
}

plotit("./Plots/copresence_network_positive.pdf",ties,"positive")
plotit("./Plots/copresence_network_negative.pdf",ties,"negative")

# positive tie YC -> RB means "YC is more likely to be present if RB is present"
# negative tie BH -> RG means "BH is less likely to be present if RG is present"


###################################################################################
# community detection
###################################################################################

# leading.eigenvector.community() only works for undirected ties 
# and non-negative weights, or no weights too.

M <- ties
"   BH BY CF DH ME OC RB RG SC WB WE YC
BH  0  1  1 -1  1  1  1 -1  0 -1  1  1
BY  1  0  1 -1  1  1  0  1  1 -1 -1  1
CF  1  1  0  0  1  0  1 -1  0 -1  0  1
DH -1 -1  0  0 -1  1 -1 -1 -1 -1  0 -1
ME  1  1  1 -1  0  1  1  0  1 -1  0 -1
OC  1  1  0  1  1  0  0  0  0  0  1  0
RB  1  0  1  0  1  0  0 -1 -1 -1  0  0
RG  0  1 -1 -1  0  0 -1  0  1 -1 -1  0
SC  0  1  0 -1  1  0  0  1  0 -1 -1  0
WB  0 -1 -1 -1 -1  0 -1 -1 -1  0 -1 -1
WE  1 -1  0  0  0  1  0 -1 -1 -1  0  0
YC  1  1  1 -1 -1  0  1  0  0 -1 -1  0"

# run the community detection algorithm for the positive ties only
pos_ties <- edgelist_from_matrix(M,tie_type="positive")
g <- make_empty_graph(n = 12, directed=F) + edges(pos_ties)

# community detection on the copresence network
result <- leading.eigenvector.community(g, weights=E(g)$weight)
ss[result[1]$`1`]  # "BH" "OC" "WE"
ss[result[2]$`2`]  # "WB"
ss[result[3]$`3`]  # "BY" "ME" "RG" "SC"
ss[result[4]$`4`]  # "CF" "RB" "YC"
ss[result[5]$`5`]  # "DH"


###################################################################################
# END
###################################################################################

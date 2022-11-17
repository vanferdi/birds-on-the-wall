##########################################################################################
# overview
##########################################################################################

# create a species-level network of who joins whom
# each tie is directed and shows if species A joins species B more than expected by chance
# ties are ruled in using glmers and correcting for multiple comparisons

# inputs:  
# adjacentized_df.RData
# circle_layout.RData

# outputs: 
# joining_columns.Rdata
# joining_glmers_output.txt
# joining_glmer_matrices.RData
# joining_network_ties.RData
# joining_network_positive.pdf
# joining_network_negative.pdf


##########################################################################################
# setup
##########################################################################################

source("./Code/function_library.R") # load common functions
require(lme4)  # package for running mixed effects regression analyses
options(scipen=999)  # run this once to disable scientific notation in print to screen


##########################################################################################
# load data
##########################################################################################

load("./Rdata/adjacentized_df.RData")
# loads one dataframe called "dfa"

# the column dfa$t2 shows whether the row was in a t2 of an event chunk
# see definition of a "joining event" in the paper

nrow(dfa)             # 6586 scans that were taken 5-minutes apart
sum(dfa$t2 == FALSE)  # there are 1289 separate runs of 5-min adjacent scans

# convert data to binary presence/absence of species
dfb <- binarize_df(dfa)

# put the t2 column back on
t2 <- dfa$t2
dfb <- cbind(dfb,t2)

# dfb should look like this:
head(dfb,3)
"      Rdate Xhr.dectime Spot BH BY CF DH ME OC RB RG SC WB WE YC    t2
1 2002-10-02    5.500000   3B  0  1  1  0  0  0  1  0  0  0  0  1 FALSE
3 2002-10-02    5.583333   3B  0  1  1  0  0  0  0  0  1  0  0  1  TRUE
5 2002-10-02    5.666667   3B  0  0  0  0  0  0  1  0  0  0  0  0  TRUE"

###################################################################################
# identify all events where the focal species joins the wall
###################################################################################

# use function get_join_events() in function_library.R

# for each species, 
# create a column that shows whether or not the species just joined the spot
# 1 = it just joined the spot, 0 = it didn't

# takes about a minute to run
BHj <- get_join_events(dfb,focal_ID="BHPA")
BYj <- get_join_events(dfb,focal_ID="BYMA")
CFj <- get_join_events(dfb,focal_ID="CFMA")
DHj <- get_join_events(dfb,focal_ID="DHPA")
MEj <- get_join_events(dfb,focal_ID="MEPA")
OCj <- get_join_events(dfb,focal_ID="OCPA")
RBj <- get_join_events(dfb,focal_ID="RBMA")
RGj <- get_join_events(dfb,focal_ID="RGMA")
SCj <- get_join_events(dfb,focal_ID="SCMA")
WBj <- get_join_events(dfb,focal_ID="WBPA")
WEj <- get_join_events(dfb,focal_ID="WEPA")
YCj <- get_join_events(dfb,focal_ID="YCPA")

# sanity check - there shouldn't be any BHj=1 where t2=FALSE
table(BHj,dfb$t2)
"BHj FALSE TRUE
   0  1289 4844
   1     0  453" 
# BH joined the wall 453 times
# NOTE: each join event can correspond to zero or more chunks, 
# depending on how many species where on the spot where BH landed

# save the resulting columns as Rdata
save(BHj,BYj,CFj,DHj,MEj,OCj,RBj,RGj,SCj,WBj,WEj,YCj, file = "./Rdata/joining_columns.RData")


###################################################################################
# look at the nubmer of join events per species
###################################################################################

# number of join events per species, in order of alphabetical species code
N_joins <- c(sum(BHj),sum(BYj),sum(CFj),sum(DHj),sum(MEj),sum(OCj),sum(RBj),sum(RGj),sum(SCj),sum(WBj),sum(WEj),sum(YCj))

Njoin <- data.frame(ss,N_joins)
Njoin <- Njoin[order(N_joins),]  # sort
Njoin
"  ss   joins
8  RG     165
4  DH     241
12 YC     245
10 WB     252
7  RB     307
11 WE     328
5  ME     388
1  BH     453
3  CF     458
9  SC     463
6  OC     501
2  BY     515"


###################################################################################
# create predictor columns
###################################################################################

# use function create_t1t2_col() in function_library.R

# the "focal species" is species A in the joining event
# here we make a new column for each species B, where 1 means species B was present at t1 and t2

# example usage:
RB_isB <- create_t1t2_col(dfb,"RB")
# 1 means that RB did qualify as a species B in chunk #14 on this row of the dataframe
# 0 means that RB didn't qualify as species B in chunk #14 on this row
# remember, these values occur on the row that corresponds to t2 of the chunk

# look at RB and RB_isB below - in row 5, RB qualifies as species B in a #14 chunk
head(cbind(dfb,RB_isB))
"      Rdate Xhr.dectime Spot BH BY CF DH ME OC RB RG SC WB WE YC    t2 RB_isB
1 2002-10-02    5.500000   3B  0  1  1  0  0  0  1  0  0  0  0  1 FALSE      0
3 2002-10-02    5.583333   3B  0  1  1  0  0  0  0  0  1  0  0  1  TRUE      0
5 2002-10-02    5.666667   3B  0  0  0  0  0  0  1  0  0  0  0  0  TRUE      0
2 2002-10-02    5.583333   3A  1  0  1  0  0  0  1  0  0  0  0  1 FALSE      0   <- not here, because it spans a chunk
4 2002-10-02    5.666667   3A  1  1  1  0  1  1  1  0  0  0  0  1  TRUE      1   <- but here yes, because it's within a chunk
6 2002-10-02    5.750000   3A  1  1  1  0  1  1  0  0  0  0  0  0  TRUE      0"


###################################################################################
# create new data frame, dfc, for running glmers on

# BHj is the variable to predict (whether BH just joined the spot)
# BH, BY, etc are the predictor variables (whether BH, BY, etc were on the spot at t1 and t2 when BH joined)

# create each new column 
BH <- create_t1t2_col(dfb,"BHPA")
BY <- create_t1t2_col(dfb,"BYMA")
CF <- create_t1t2_col(dfb,"CFMA")
DH <- create_t1t2_col(dfb,"DHPA")
ME <- create_t1t2_col(dfb,"MEPA")
OC <- create_t1t2_col(dfb,"OCPA")
RB <- create_t1t2_col(dfb,"RBMA")
RG <- create_t1t2_col(dfb,"RGMA")
SC <- create_t1t2_col(dfb,"SCMA")
WB <- create_t1t2_col(dfb,"WBPA")
WE <- create_t1t2_col(dfb,"WEPA")
YC <- create_t1t2_col(dfb,"YCPA")

dfc <- cbind(dfb[1:3],BHj,BYj,CFj,DHj,MEj,OCj,RBj,RGj,SCj,WBj,WEj,YCj,BH,BY,CF,DH,ME,OC,RB,RG,SC,WB,WE,YC)
head(dfc)


###################################################################################
# run glmers
###################################################################################

# First, let's look at an example testing one species pair 

# are BH join events predicted by the presence of BY on the spot?
full <- glmer(BHj ~ BY + CF + DH + ME + OC + RB + RG + SC + WB + WE + YC + (1|Spot) + (1|Xhr.dectime) + (1|Rdate), family="binomial", nAGQ=0, data=dfc)
reduced <- glmer(BHj ~ CF + DH + ME + OC + RB + RG + SC + WB + WE + YC + (1|Spot) + (1|Xhr.dectime) + (1|Rdate), family="binomial", nAGQ=0, data=dfc)
# remove BY and determine whether it was a significant predictor of BHj
a <- anova(full,reduced)
"       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
reduced 14 3082.8 3177.9 -1527.4   3054.8                           
full    15 3080.7 3182.6 -1525.4   3050.7 4.0176      1    0.04503 "  

# get the direction of the effect
summary(full)$coefficients["BY","Estimate"]  # 0.2586285

# get the pvalue 
a$`Pr(>Chisq)`[2]  # 0.04502889

# after we run all of these models, we correct for multiple comparisons


###################################################################################
# run all glmers

# loop through all lmer comparisons, print the results to screen, save effects and pvals to matrix

join_pvals <- matrix(0,12,12) # save the pvalue here
join_est <- matrix(0,12,12)   # save the direction of the effect here
colnames(join_pvals) <- ss
rownames(join_pvals) <- ss
colnames(join_est) <- ss
rownames(join_est) <- ss

# print results to file - remember to run sink() to end this process and print to console again
sink(file = "./Code/Output/joining_glmers_output.txt")

for (i in 1:12) {
	focal <- ss[i]  # get current focal species
	print(noquote("==========================="))
	print(noquote(paste("CURRENT FOCAL SPECIES IS",focal)))
	
	# create the 12 full formulas, one for each focal species as the dependent variable
	# use the 11 non-focal species as the predictor species, plus random effects for spot, time and date
	predictor_species <- paste0(ss[-i], sep="+", collapse="")  # remove the focal species from predictors
	f <- formula(paste(paste0(focal,"j"),"~", predictor_species,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	full <- glmer(f, family="binomial", nAGQ=0, data=dfc)
	fs <- summary(full)
	print(fs)
	
	# given the current focal species, create the 11 reduced models
	for (j in 1:11) {
		x <- ss[-i]  # remove the focal
		one_to_remove <- x[j]  # and remove the predictor species being reduced
		x <- x[-j]
		predictor_species2 <- paste0(x, sep="+", collapse="")
		r <- formula(paste(paste0(focal,"j")," ~", predictor_species2,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
		
		print(noquote("==========================="))
		print(noquote(paste("=== removing species",one_to_remove,"===")))
		
		# conduct anova on the full and reduced model
		reduced <- glmer(r, family="binomial", nAGQ=0, data=dfc)
		a <- anova(full,reduced)
		print(a)  # will print to sink file instead of screen
		
		# get pvalue and save to join_pvals
		pval <- a$`Pr(>Chisq)`[2]
		join_pvals[focal,one_to_remove] <- round(pval,10)
		
		# get estimate of effect and save to join_est
		value <- fs$coefficients[one_to_remove,][1]
		join_est[focal,one_to_remove] <- value
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

# save the resulting matrices as Rdata
save(join_est, join_pvals, file = "./Rdata/joining_glmer_matrices.RData")


###################################################################################
# load results from Rdata
###################################################################################

load("./Rdata/joining_glmer_matrices.RData")

join_pvals	
join_est


###################################################################################
# extract significant ties
###################################################################################

# gets the significant ties with FDR correction for multiple comparisons

# save ties as a matrix where:
# 1 = significant affiliative tie
# -1 = significant avoidant tie
# 0 = tie not significant

# each entry in the matrix is for the row species as the focal species
# because in the glmers, presence of row species is predicted by presence of column species

join_ties <- get_sig_ties(join_pvals, join_est)
"  BH BY CF DH ME OC RB RG SC WB WE YC
BH  0  0  1  0  0  0  1  0  0  0  0  1
BY  0  0  0  0  0  0  1  0  0 -1  0  0
CF  0  0  0  0  0  0  1  0  0  0  0  0
DH  0  0  0  0  0  0  0  0  0  1  1  0
ME  0  0  0  0  0  0  1 -1  0  0  0  1
OC  0  0  1  1  0  0  0 -1  0  1  1  0
RB  0  0  1  0  0  0  0  0  0 -1  0  0
RG  0  1  0  0  0  0  0  0  1  0  0  0
SC  0  1  0  0  0  0  0  0  0 -1  0  0
WB  0  0  0  0  0  0  0  0  0  0  1  0
WE  0  0  0  0  0  0  0  0  0  1  0  0
YC  0  0  1  0 -1  0  1  0  0  0  0  0"

# examples of how to read the join_ties matrix:
# ME is significantly more likely to join if YC is present
# YC is significantly less likely to join if ME is present

save(join_ties, file="./Rdata/joining_network_ties.RData")


###################################################################################
# plot network
###################################################################################

require(igraph)

# load the main layout
load("./Rdata/circle_layout.RData")  # loads circle_layout
circle_layout

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

plotit("./Plots/joining_network_positive.pdf",join_ties,"positive")
plotit("./Plots/joining_network_negative.pdf",join_ties,"negative")

# positive tie BH -> RB means "BH is more likely to join the spot if RB is present on the spot"
# negative tie RB -> WB means "RB is less likely to join the spot if WB is present on the spot"


###################################################################################
# community detection
###################################################################################

M <- join_ties

# run the community detection algorithm for the positive ties only
pos_ties <- edgelist_from_matrix(M,tie_type="positive")
g <- make_empty_graph(n = 12, directed=F) + edges(pos_ties)

# community detection on the network of joining behavior
result <- leading.eigenvector.community(g, weights=E(g)$weight)
ss[result[1]$`1`]  # "BH" "CF" "ME" "RB" "YC"
ss[result[2]$`2`]  # "DH" "OC" "WB" "WE"
ss[result[3]$`3`]  # "BY" "RG" "SC"


###################################################################################
# END
###################################################################################

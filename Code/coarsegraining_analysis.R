##########################################################################################
# overview
##########################################################################################

# See if any categorization systems describe joining behavior better than species categories do

# input: 
# adjacentized_df.RData
# joining_columns.RData

# output:
# coarsegraining_glmers_output.txt

source("./Code/function_library.R")  # load common functions
require(lme4)  # package for running mixed effects regression analyses
options(scipen=999)  # run this once to disable scientific notation in print to screen


##########################################################################################
# load data
##########################################################################################

load("./Rdata/adjacentized_df.RData")  # loads one dataframe called "dfa"
dfb <- binarize_df(dfa)  # convert data to binary presence/absence of species
t2 <- dfa$t2  # put the t2 column back on dfb
dfb <- cbind(dfb,t2)
nrow(dfb)  # 6586, correct

# load 12 columns that were created in joining_network.R
# loads variables BHj,BYj,CFj,DHj,MEj,OCj,RBj,RGj,SCj,WBj,WEj,YCj
load("./Rdata/joining_columns.RData")


##########################################################################################
# create predictor columns for species
##########################################################################################

# create a column per species showing whether it was there at t1 and t2
# 1 = true, 0 = false, and the number is labelled on the t2 row
# aka 1 = "I'm here now and I was also here in the previous row (and both rows are in a valid chunk)"

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

t1t2_cols <- data.frame(BH,BY,CF,DH,ME,OC,RB,RG,SC,WB,WE,YC)

# create a new dataframe
# BHj is a binary column saying whether BH just joined the wall or not
# BH is a binary column saying whether BH was there at t1 and t2 (where current row is t2)
dfc <- cbind(dfb[1:3],BHj,BYj,CFj,DHj,MEj,OCj,RBj,RGj,SCj,WBj,WEj,YCj,BH,BY,CF,DH,ME,OC,RB,RG,SC,WB,WE,YC,t2)
head(dfc)

##########################################################################################

# look at the total number of join events per species
N_joins <- colSums(dfc[4:15])
Njoin <- data.frame(ss,N_joins)
Njoin <- Njoin[order(N_joins),]  # sort
Njoin
"  ss N_joins
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

# look at how many times each species was there at a t1 and t2
N_theres <- colSums(dfc[16:27])
Nthere <- data.frame(ss,N_theres)
Nthere <- Nthere[order(N_theres),]  # sort
Nthere
"   
ss N_theres
RG      235
YC      292
DH      362
SC      780
OC     1051
WB     1094
WE     1126
BH     1127
BY     1307
RB     1323
CF     1380
ME     1713"


##########################################################################################
# create predictor columns for various coarse-grainings
##########################################################################################

# takes a set of species and tallies them up into a larger category
# requires a binary dataframe with columns named by the short code (will use dfc)

category_columns <- function(df) {
	
	##############################
	macaw_L <- as.numeric(df$BY==1 | 
						  	df$RG==1 | 
						  	df$SC==1)
	
	macaw_N <- as.numeric(df$BH==1 | 
						  	df$CF==1 | 
						  	df$DH==1 | 
						  	df$ME==1 | 
						  	df$OC==1 | 
						  	df$RB==1 | 
						  	df$WB==1 | 
						  	df$WE==1 | 
						  	df$YC==1)
	
	##############################
	clade_SP <- as.numeric(df$BH==1 | 
						   	df$OC==1 | 
						   	df$WB==1)
	
	clade_LP <- as.numeric(df$ME==1 | 
						   	df$YC==1)
	
	clade_K <- as.numeric(df$WE==1 | 
						  	df$DH==1)
	
	clade_SM <- as.numeric(df$RB==1 | 
						   	df$CF==1)
	
	clade_LM <- as.numeric(df$BY==1 | 
						   	df$SC==1 | 
						   	df$RG==1)
	
	##############################
	size_S <- as.numeric(df$BH==1 | 
						 	df$DH==1 |
						 	df$OC==1 | 
						 	df$WB==1 |
						 	df$WE==1)
	
	size_M <- as.numeric(df$CF==1 | 
						 	df$ME==1 |
						 	df$RB==1 | 
						 	df$YC==1)
	
	size_L <- as.numeric(df$BY==1 | 
						 	df$RG==1 |
						 	df$SC==1)
	
	##############################
	back_G <- as.numeric(df$BH==1 | 
						 	df$CF==1 |
						 	df$DH==1 |
						 	df$ME==1 |
						 	df$OC==1 |
						 	df$RB==1 |
						 	df$WB==1 |
						 	df$WE==1 |
						 	df$YC==1)
	
	back_B <- as.numeric(df$BY==1)
	
	back_R <- as.numeric(df$RG==1 | 
						 	df$SC==1)
	
	##############################
	head_B <- as.numeric(df$BH==1 | 
						 	df$BY==1)
	
	head_G <- as.numeric(df$CF==1 | 
						 	df$DH==1 |
						 	df$ME==1 |
						 	df$RB==1 |
						 	df$WE==1 |
						 	df$YC==1)
	
	head_R <- as.numeric(df$RG==1 | 
						 	df$SC==1)
	
	head_O <- as.numeric(df$WB==1)
	
	head_K <- as.numeric(df$OC==1)
	
	##############################
	face_B <- as.numeric(df$BH==1)
	
	face_O <- as.numeric(df$WB==1)
	
	face_Y <- as.numeric(df$OC==1 | 
						 	df$YC==1)
	
	face_W <- as.numeric(df$BY==1 | 
						 	df$CF==1 |
						 	df$RB==1 |
						 	df$RG==1 |
						 	df$SC==1)
	
	face_G <- as.numeric(df$DH==1 | 
						 	df$ME==1 |
						 	df$WE==1)
	
	##############################
	# save all columns in a list structure  - creates 23 new columns
	
	x <- list(macaw_L,macaw_N,
			  clade_SP,clade_LP,clade_K,clade_SM,clade_LM,
			  size_S,size_M,size_L,
			  back_G,back_B,back_R,
			  head_B,head_G,head_R,head_O,head_K,
			  face_B,face_O,face_Y,face_W,face_G)
	return(x)
}

# create the category columns for the size categories that are relative to the size of the focal species
category_columns_byfocal <- function(df,focal) {
	
	# get the set of species that go in each category
	larger_set <- larger_smaller_sets(focal)[[1]]   # get all of the species that are larger the focal
	smaller_set <- larger_smaller_sets(focal)[[2]]  # get all of the species that are smaller the focal
	adja_set <- same_diff_sets(focal)[[1]]
	diff_set <- same_diff_sets(focal)[[2]]
	
	# get the corresponding t1t2 columns for the species in this set
	L <- t1t2_cols[larger_set]
	S <- t1t2_cols[smaller_set]
	A <- t1t2_cols[adja_set]
	D <- t1t2_cols[diff_set]
	
	# create the category columns for the larger and smaller categories
	larger <- as.numeric(rowSums(L)>0)  # if any of these columns have a one in them, larger=1, ortherwise larger=0
	smaller <- as.numeric(rowSums(S)>0)
	sim_size <- as.numeric(rowSums(A)>0)
	dif_size <- as.numeric(rowSums(D)>0)

	x <- list(larger,smaller,sim_size,dif_size)
	return(x)
}

# now to actually make these columns,
# you need to remove the focal species data from those columns
# for example, if we want to now whether BY is predicted by the "large macaw" category
# we should only tally up the presence of SC and RG as "large macaw" 
# and exclude BY's own presence from that tally

# so we make one dataframe per focal species
# and tally up the categories, customized to that focal species


###################################################################################
# create 12 dataframes, one per focal species
###################################################################################

# focal must be the short species code, because it refers to the species columns in dfb

make_focal_df <- function(df,focal) {  # ex: focal="BH" (needs the short code)
	
	# create the new predictor column (focal species joined vs not) 
	focal_joins <- get_join_events(df,focal)  
	
	# you've got to remove the focal species data before you create the categories
	dtemp <- df
	dtemp[,focal] <- 0  # zero out focal species data
	
	# create all categories
	macaw_L <- category_columns(dtemp)[[1]]
	macaw_N <- category_columns(dtemp)[[2]]
	clade_SP <- category_columns(dtemp)[[3]]
	clade_LP <- category_columns(dtemp)[[4]]
	clade_K <- category_columns(dtemp)[[5]]
	clade_SM <- category_columns(dtemp)[[6]]
	clade_LM <- category_columns(dtemp)[[7]]
	size_S <- category_columns(dtemp)[[8]]
	size_M <- category_columns(dtemp)[[9]]
	size_L <- category_columns(dtemp)[[10]]
	back_G <- category_columns(dtemp)[[11]]
	back_B <- category_columns(dtemp)[[12]]
	back_R <- category_columns(dtemp)[[13]]
	head_B <- category_columns(dtemp)[[14]]
	head_G <- category_columns(dtemp)[[15]]
	head_R <- category_columns(dtemp)[[16]]
	head_O <- category_columns(dtemp)[[17]]
	head_K <- category_columns(dtemp)[[18]]
	face_B <- category_columns(dtemp)[[19]]
	face_O <- category_columns(dtemp)[[20]]
	face_Y <- category_columns(dtemp)[[21]]
	face_W <- category_columns(dtemp)[[22]]
	face_G <- category_columns(dtemp)[[23]]
	
	larger <- category_columns_byfocal(dtemp,focal)[[1]]
	smaller <- category_columns_byfocal(dtemp,focal)[[2]]
	sim_size <- category_columns_byfocal(dtemp,focal)[[3]]
	dif_size <- category_columns_byfocal(dtemp,focal)[[4]]
	
	dtemp <- cbind(dtemp, focal_joins,
				   macaw_L,macaw_N,
				   clade_SP,clade_LP,clade_K,clade_SM,clade_LM,
				   size_S,size_M,size_L,
				   back_G,back_B,back_R,
				   head_B,head_G,head_R,head_O,head_K,
				   face_B,face_O,face_Y,face_W,face_G,
				   larger,smaller,sim_size,dif_size)
	
	dtemp[,focal] <- df[,focal]  # put all focal species data back
	
	return(dtemp)
}

# make all of the data frames
# use the new dataframe, dfc, where species cols (ex:BH) show t1t2 presence
df_BH <- make_focal_df(dfc,"BH")
df_BY <- make_focal_df(dfc,"BY")
df_CF <- make_focal_df(dfc,"CF")
df_DH <- make_focal_df(dfc,"DH")
df_ME <- make_focal_df(dfc,"ME")
df_OC <- make_focal_df(dfc,"OC")
df_RB <- make_focal_df(dfc,"RB")
df_RG <- make_focal_df(dfc,"RG")
df_SC <- make_focal_df(dfc,"SC")
df_WB <- make_focal_df(dfc,"WB")
df_WE <- make_focal_df(dfc,"WE")
df_YC <- make_focal_df(dfc,"YC")

# sanity check: these two columns should be identical
unique(df_BH$BHj == df_BH$focal_joins)  # TRUE, good

# save all dataframes in this list, indexable by same species code order
dfs <- list(df_BH,df_BY,df_CF,df_DH,df_ME,df_OC,df_RB,df_RG,df_SC,df_WB,df_WE,df_YC)


###################################################################################
# run the GLMERs
###################################################################################

evaluate_models <- function(dfs,focal) {  # focal = BHPA or BH
    # convert focal species to number 1-12 (1 = BHPA ... 12 = YCPA)
	species <- ssi(focal)
	
	# print the focal species
	print(noquote(s[species]))
	
	macaw_string <- "macaw_L + macaw_N +"
	clade_string <- "clade_SP + clade_LP + clade_K + clade_SM + clade_LM +"
	size_string <- "size_S + size_M + size_L +"
	back_string <- "back_G + back_B + back_R +"
	head_string <- "head_B + head_G + head_R + head_O + head_K +"
	face_string <- "face_B + face_O + face_Y + face_W + face_G +"
	relsize1_string <- "larger + smaller +"
	relsize2_string <- "sim_size + dif_size +"
	
	# get the focal species dataframe and the focal species column
	df <- dfs[[species]]   
	foc <- ss[species]  # short code of species
	
	# for the species-level model, remove the focal species from the list of predictor species
	predictor_species <- paste0(ss[-ssi(foc)], sep="+", collapse="")  
	
	species_formula <- formula(paste("focal_joins~", predictor_species,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	macaw_formula <- formula(paste("focal_joins~",macaw_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	clade_formula <- formula(paste("focal_joins~",clade_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	size_formula <- formula(paste("focal_joins~",size_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	back_formula <- formula(paste("focal_joins~",back_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))  
	head_formula <- formula(paste("focal_joins~",head_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))  
	face_formula <- formula(paste("focal_joins~",face_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	relsize1_formula <- formula(paste("focal_joins~",relsize1_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	relsize2_formula <- formula(paste("focal_joins~",relsize2_string,"(1|Spot) + (1|Xhr.dectime) + (1|Rdate)"))
	
	species_model <- glmer(species_formula, family="binomial", nAGQ=0, data=df)
	macaw_model <- glmer(macaw_formula, family="binomial", nAGQ=0, data=df)
	clade_model <- glmer(clade_formula, family="binomial", nAGQ=0, data=df)
	size_model <- glmer(size_formula, family="binomial", nAGQ=0, data=df)
	back_model <- glmer(back_formula, family="binomial", nAGQ=0, data=df)
	head_model <- glmer(head_formula, family="binomial", nAGQ=0, data=df)
	face_model <- glmer(face_formula, family="binomial", nAGQ=0, data=df)
	largersmaller_model <- glmer(relsize1_formula, family="binomial", nAGQ=0, data=df)
	simdifsize_model <- glmer(relsize2_formula, family="binomial", nAGQ=0, data=df)
	
	a <- anova(species_model,macaw_model,clade_model,size_model,back_model,head_model,face_model,largersmaller_model,simdifsize_model)
	
	# create your own output table, sorted by Akaike weights
	model <- rownames(a)
	param <- a$Df
	loglike <- a$logLik
	AIC <- a$AIC
	weight <- aicw(a$AIC)
	delta <- AIC-min(AIC)
	
	dx <- data.frame(model,param,loglike,AIC,delta,weight)
	dx <- dx[order(AIC),]
	
	# print stuff to screen so you can sink it to a log file
	print(a)  # anova
	print(dx) # your AIC table
	
	if (dx[1,]$model == "macaw_model") { best <- macaw_model }
	if (dx[1,]$model == "size_model") { best <- size_model }
	if (dx[1,]$model == "back_model") { best <- back_model }
	if (dx[1,]$model == "clade_model") { best <- clade_model }
	if (dx[1,]$model == "head_model") { best <- head_model }
	if (dx[1,]$model == "face_model") { best <- face_model }
	if (dx[1,]$model == "species_model") { best <- species_model }
	if (dx[1,]$model == "largersmaller_model") { best <- largersmaller_model }
	if (dx[1,]$model == "simdifsize_model") { best <- simdifsize_model }
	
	print(noquote(""))
	print(noquote(paste0(s[species]," >>> BEST FIT MODEL >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")))
	print(summary(best))
	
	print(noquote(""))
	print(noquote(paste0(s[species]," >>> SPECIES MODEL >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")))
	print(summary(species_model))
	
	# return the number of significant categories in the best-fit model
	x <- summary(best)$coefficients[,4]
	ncats <- sum(x[2:length(x)]<0.05) 
	return(ncats)
}

# print results to file - remember to run sink() to end this process and print to console again
sink(file = "./Code/Output/coarsegraining_glmers_output.txt")

ncats <- c()  # save number of significant categories in the best-fit category system here
for (i in 1:12) {
	print(noquote(""))
	print(noquote("==========================="))
	print(noquote(paste("CURRENT FOCAL SPECIES IS",ss[i])))
	print(noquote("==========================="))
	
	ncat <- evaluate_models(dfs,ss[i]) # run model evaluation and also save the ncats
	ncats <- c(ncats,ncat)
	
	if (i == 12) { 
		sink()  # close the sink connection
	}
}

# the rank-deficient warnings are fine, they refer to cases where the focal species was the only
# member of a category, and since we removed the focal species from the predictors,
# that means there's nothing to predict using that category, so it gets dropped as a predictor

# save variable ncats to use in the centrality and complexity analysis
save(ncats, file = "./Rdata/ncats.RData")


###################################################################################
# Results: see coarsegraining_glmers_output.txt
###################################################################################

# report all top models where delta AIC from the best model is < 2


###################################################################################
# END
###################################################################################

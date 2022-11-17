
# the 12 species, in alphabetical order
s <- c("BHPA","BYMA","CFMA","DHPA","MEPA","OCPA","RBMA","RGMA","SCMA","WBPA","WEPA","YCPA")  # long codes
ss <- c("BH","BY","CF","DH","ME","OC","RB","RG","SC","WB","WE","YC")                         # short codes

# body mass of each species, in order of s
# this data comes from [cited in paper]
grams <- c(247,1125,430,108,610,178,370,1250,1015,155,157,510)

# the list of standardized spots
spots <- c('1A','1B','1C','2A','2B','2C','3A','3B','3B1','3B2','3C') 

# remove rows in df that are from non-standardized spot  IDs
filter_spots <- function(df,spots) {
	df <- subset(df,Spot %in% spots)
	df$Spot <- factor(df$Spot)  # drop all factors that have been removed
	return(df)
}

# example usage:
# db <- filter_spots(db,spots)

# these are two colors I'll use across different plots and figures
myblue <- "#3ec9f7"
myred <- "#cf1c08"

# convert dataframe to binary presence/absence of species per row
binarize_df <- function(df) {
	df2 <- df[4:15]
	df2[df2 > 0] <- 1
	df2 <- cbind(df[1:3],df2)
	return(df2)
}

# example usage:
# db <- binarize_df(db)

# takes a matrix of pvalues and corrects for multiple comparisons with false discovery rate
# output is a binary matrix showing which cells are significant after correction
FDR_correction <- function(matrix_of_pvals,verbose=FALSE) {
	row <- rep(ss,12)
	col <- rep(ss, each=12)
	diag(matrix_of_pvals) <- NA 
	pval <- as.numeric(matrix_of_pvals)
	dx <- data.frame(row,col,pval)
	dx <- dx[order(pval),]
	dx <- na.omit(dx)
	rank <- seq(1,nrow(dx)) 
	dx <- cbind(dx,rank)
	cv <- (rank/nrow(dx))*0.05
	dx <- cbind(dx,cv)
	discoveries <- dx$pval <= dx$cv
	dx <- cbind(dx,discoveries)
	inR <- p.adjust(dx$pval,method="BH")  # inR < 0.05 is a discovery
	dx <- cbind(dx,inR)
	if(verbose == TRUE) { print(dx) }  # print the table
	join_sigs <- matrix(FALSE,12,12)
	colnames(join_sigs) <- ss
	rownames(join_sigs) <- ss
	for (i in 1:nrow(dx)) { join_sigs[as.character(dx$row[i]),as.character(dx$col[i])] <- dx$discoveries[i] }
	return(join_sigs)
}

# example usage:
# result <- FDR_correction(M_pvals)

# use an FDR correction to get significant ties from a matrix of pvalues and a matrix of estimates
# output is a matrix where 0 = not significant, 1 = significantly affiliative, -1 = significantly avoidant
get_sig_ties <- function(M_pvals,M_est) {
	M_sig <- FDR_correction(M_pvals)
	M_ties <- M_est
	M_ties[!M_sig] <- 0
	M_ties <- sign(M_ties)
	return(M_ties)
}

# convert a species code to it's numeric index in s
ssi <- function(species_code) {
	if (nchar(species_code)==4) { s <- c("BHPA","BYMA","CFMA","DHPA","MEPA","OCPA","RBMA","RGMA","SCMA","WBPA","WEPA","YCPA") }
	if (nchar(species_code)==2) { s <- c("BH","BY","CF","DH","ME","OC","RB","RG","SC","WB","WE","YC") }
	return(which(s==species_code))
}

# convert a species code to its numeric index in the dataframe (columns 4:15)
# accepts either the long (BYMA) or short (BY) form of the species code
species_to_column <- function(species_code) {
	if (nchar(species_code)==4) { s <- c("BHPA","BYMA","CFMA","DHPA","MEPA","OCPA","RBMA","RGMA","SCMA","WBPA","WEPA","YCPA") }
	if (nchar(species_code)==2) { s <- c("BH","BY","CF","DH","ME","OC","RB","RG","SC","WB","WE","YC") }
	return(which(s==species_code)+3)
}

# create array of social network ties formatted for use in graphing social networks with igraph
edgelist_from_matrix <- function(matrix,tie_type="all") {
	if (tie_type == "all") { t <- which(matrix != 0,arr.ind=TRUE) }
	if (tie_type == "positive") { t <- which(matrix > 0,arr.ind=TRUE) }
	if (tie_type == "negative") { t <- which(matrix < 0,arr.ind=TRUE) }
	return(as.numeric(t(t)))
}

# grabs all tie weights from a matrix, but only for the ties supplied in edgelist
weightlist_from_matrix <- function(M,edgelist) {
	tie_weights <- c()
	for (i in 1:(length(edgelist)/2)) {  # length(edgelist) will always be even
		current_pair <- edgelist[1:2]  # get the current edge (it's a pair of 2 nodes)
		edgelist <- edgelist[3:length(edgelist)]  # remove that pair
		a <- M[current_pair[1],current_pair[2]]
		tie_weights <- c(tie_weights,a)
	}
	return(tie_weights)
}

# example usage:
# weightlist <- weightlist_from_matrix(M,edgelist_from_matrix(M))

# get baseline probability that each species appears in an observation slot
# based on total number of appearances for all observations in the dataframe regardless of spot
get_baseline <- function(df) {
	N <- nrow(df)
	species_n <- c()
	for (i in s) {
		col <- species_to_column(i)
		n <- sum(df[,col] > 0)  # number of occurances
		species_n <- c(species_n,n)
	}
	species_p <- species_n/N
	species_n
	species_p
	lookup <- data.frame(s,species_n,species_p)
	return(lookup)
}

# for one species in the binary adjacentized dataframe
# creates a new column showing whether or not the species just joined the spot on row i+1
# 1 = it just joined, 0 = not
get_join_events <- function(df,focal_ID) {
	foc <- species_to_column(focal_ID)
	joined <- c(0)  # save joins here, per each row in df - start it with zero because the first row is t2=FALSE
	
	for (i in 1:(nrow(df)-1)) {                         
		rowt1 <- df[i,]    # extract entire row for time t1
		rowt2 <- df[i+1,]  # extract entire row for time t2
		
		# make sure join happened within a valid chunk
		if ( rowt2$t2 == TRUE ) {  # we're currently in a valid chunk
			# check if focal joined or not
			if (rowt1[foc] == 0 && rowt2[foc] == 1) { joined <- c(joined,1) } else { joined <- c(joined,0) }
		}
		if ( rowt2$t2 == FALSE ) {  # we're not in a valid chunk
			joined <- c(joined,0)   # mark it as not a join, regardless
		}
	}
	return(joined)
}

# for one species in the binary adjacentized dataframe
# creates a new column for that species, where 1 (on t2 row) means the species was present at t1 and t2
# explained more in joining_network.R
create_t1t2_col <- function(df,species_ID) {
	present <- rep(0,nrow(df))
	print(length(present))
	
	# get presence data for current species
	data <- df[,species_to_column(species_ID)]  # format is a vector, index like data[1], not data[1,]
	
	for (i in 1:(nrow(df)-1)) {
		if (data[i] == 1 && data[i+1] == 1) {  # the current species was there at t1 and t2
			# then set present[i+1] to 1 if we're in a valid chunk
			if (df$t2[i+1]  == TRUE) {  present[i+1] <- 1 } # otherwise leave the default 0
		}
	}
	return(present)
}

# compute the Akaike weights
aicw <- function(AICs) {
	deltas <- AICs-min(AICs)
	x <- exp(-0.5*deltas)
	return(x/sum(x))
}

# returns a list of 
# 1) larger: all species that are larger than the focal species
# 2) smaller: all species that are smaller than the focal species
larger_smaller_sets <- function(focal) {  # focal can be lond or short code
	mass <- data.frame(ss,grams)
	focal_mass <- mass[ssi(focal),]$grams
	larger <- ss[mass$grams > focal_mass]
	smaller <- ss[mass$grams < focal_mass]
	return(list(larger,smaller))
}

# returns a list of
# 1) same: the set of species that are adjacent sizes to the focal (will be 1 or 2 species)
# 2) diff: the set of species that are not adjacent in size to the focal (these are all other species)
same_diff_sets <- function(focal) {
	mass <- data.frame(ss,grams)
	mass <- mass[order(grams),]
	focal_index <- which(mass$ss == focal)  # get index of the focal species
	if (focal_index == 1) {  # edge case 1
		same <- mass$ss[2]
		diff <- mass$ss[3:12]
	}
	if (focal_index == 12) {  # edge case 2
		same <- mass$ss[11]
		diff <- mass$ss[1:10]
	}
	else {  # general case
		sames <- c(focal_index-1,focal_index+1)  # get indeces of the two adjacent species
		same <- mass$ss[sames]
		diff <- mass$ss[-c(sames,focal_index)]
	}
	
	return(list(same,diff))
}

###################################################################################
# END
###################################################################################

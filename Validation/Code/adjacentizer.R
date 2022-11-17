##########################################################################################
# overview
##########################################################################################

# This code creates a new dataframe containing only the scans that were taken in 5-minute runs
# If a scan had another scan that was either five minutes before it or after it, that scan is in this new dataset.

# output:
# dfa_binary - dataframe of 5-min adjacent scans, with species presence in binary
# dfa_counts - dataframe of 5-min adjacent scans, with species presence in terms of N birds
# both are Rdata object, stored in adjacentized_dfs.Rdata

# This code modifies the original dataframe (evendays_dataset.csv) like this:
# 1) sorts dataframe so all spots per date are grouped
# 2) adds a column saying whether the row is part of an adjacent chunk
#    adjacent chunk means row[i] and row[i+1] are 5 minutes apart
# see below for in-code explanations plus examples

# example:
# df$adjacent[3] == FALSE, df$adjacent[4] == TRUE, df$adjacent[5] == TRUE 
# here, 3 isn't adjacent to 2 or 3, but 4 and 5 are adjacent to one another

# load common functions
source("./Code/function_library.R")

##################################################################################
# load data
##################################################################################

db <- read.csv("./Data/odddays_dataset.csv", header = TRUE, sep = ",")
nrow(db)  # 7682 rows

# remove data from unstandardized spots
db <- filter_spots(db,spots)
nrow(db)  # 7673, correct

# change columns to the shorter species codes
colnames(db)[4:15] <- ss  # change cols to shorter species codes


##################################################################################
# some support functions
##################################################################################

is.increasing <- function(x) { all(diff(x) > 0) }

# example usage:
is.increasing(c(5,3,2,4,1))  # FALSE
is.increasing(c(1,2,3,4,5))  # TRUE

# convert a number of minutes to its dectime representations
minute_to_dectime <- function(minute) {
	return((minute/60))
}

# example usage:
minute_to_dectime(5)  #  5 minutes is 0.08333333 in dectime

# check whether two scans are 5 minutes apart or not
is.5min <- function(t1,t2) { # t1 and t2 in dectime format
	diff <- t2-t1
	# due to rounding errors, allow small window for 5-min hits
	if (diff > minute_to_dectime(4.9) & diff < minute_to_dectime(5.1)) {
		result <- TRUE
	} else { result <- FALSE }
	return(result)
}

# example usage:
is.5min(5.366667,5.450000)  # TRUE
is.5min(5.366667,5.850000)  # FALSE

##################################################################################
# check if all times with an observation date and spot are increasing
# if so, this is a regularity in the data that we can exploit
# ANSWER: yes all times are increasing

check_increasing <- function(df) {
	result <- c()
	unique_dates <- unique(df$Rdate)
	for (d in unique_dates) {
		date_subset <- subset(df,Rdate==d)
		unique_spots <- unique(date_subset$Spot) # grab zones for that date
		for (s in unique_spots) {
			spot_subset <- subset(date_subset,Spot==s)
			result <- c(result,is.increasing(spot_subset$Xhr.dectime))
		}
	}
	return(result)
}

unique(check_increasing(db))  # TRUE


##################################################################################
# MAIN CODE
##################################################################################

adjacentize_dataframe <- function(df) {
	
	# create new empty dataframe with same names as original dataframe
	new_df <- data.frame(matrix(ncol=ncol(df), nrow=0))
	colnames(new_df) <- names(df)
	
	# new column for boolean variables
	adjacent <- c()  # TRUE if row is involved in a chunk (t1 or t2 of 5min-adjacent time slots)
	t2 <- c()  # TRUE if this row is a t2 in a chunk (some are only t1) (use to count # of chunks)
	
	unique_dates <- unique(df$Rdate)[1:10]  # for development
	unique_dates <- unique(df$Rdate)
	
	# for each date
	for (d in unique_dates) {
		print(d)
		date_subset <- subset(df,Rdate==d)
		unique_spots <- unique(date_subset$Spot) # grab zones for that date
		
		# for each spot in that date
		for (s in unique_spots) {
			spot_subset <- subset(date_subset,Spot==s)
			
			# save this data and spot chunk to the new df
			new_df <- rbind(new_df,spot_subset)
			
			#################################################################
			# add the variable showing where the 5min-adjacent entries are
			
			# grab all the times (we already know they are all increasing)
			t <- spot_subset$Xhr.dectime
			
			# get time differences, x[1] is diff between t[2]-t[1]
			x <- diff(t) > minute_to_dectime(4.9) & diff(t) < minute_to_dectime(5.1)
			x <- c(FALSE,x) # add FALSE at front so x[i] will be properly positioned
			adjacent <- c(adjacent,x)
			t2 <- c(t2,x)
			
			#################################################################
		}
	}
	# all the rows corresponding to t2 in the chunk are labeled 1 if they're adjacent
	# now go back through and label change all the labels of t1 adjacents to TRUE also
	for (i in 2:length(adjacent)) {
		if (adjacent[i]==TRUE) {  # if i is TRUE
			adjacent[i-1] <- TRUE  # then set i-1 to TRUE also
		}
	}
	
	new_df <- cbind(new_df,adjacent,t2)
	return(new_df)
}

# create the adjacentized dataframe
dfa <- adjacentize_dataframe(db)

nrow(dfa)  # 7673 - same as original, good

# subset so only rows involved in an adjacent chunk are kept in the dataframe
dfa <- subset(dfa,adjacent==TRUE)
nrow(dfa)  # 6919 scans that are 5-min adjacent
head(dfa)
"       Rdate Xhr.dectime Spot BH BY CF DH ME OC RB RG SC WB WE YC adjacent    t2
1  2002-10-01    6.283333   2C 52 13  7  0 13  0 10  0  0  0  0  0     TRUE FALSE
3  2002-10-01    6.366667   2C 80 19  7  0 29  0  8  0  0  0  0  2     TRUE  TRUE
6  2002-10-01    6.450000   2C 47 16  3  0 55  0  1  1  5  0  0  4     TRUE  TRUE
9  2002-10-01    6.533333   2C  1  0  0  0  0  0  0  0  1  0  0  0     TRUE  TRUE
13 2002-10-01    6.616667   2C 14  1  0  0 39  0  0  0  0  0  0  0     TRUE  TRUE
16 2002-10-01    6.700000   2C 13  3  0  0 16  2  0  2  0  0  0  0     TRUE  TRUE"

# how many event chunks are there? A chunk is a pair of t1 and t2 rows.
sum(dfa$t2==TRUE)  # 5589 chunks

# how many scans are involved in a chunk?
sum(dfa$adjacent==TRUE)  # 6919

# save the new adjacentized dataframe as Rdata
save(dfa,file="./Rdata/adjacentized_df.Rdata")  


##################################################################################
# END
##################################################################################

##########################################################################################
# overview
##########################################################################################

# get the total number of times each species was present in a scan
# and compute the percentage of those scans in which it was the only species present
# will show how anti-social each species is in its clay lick usage

# input:
# evendays_dataset.csv


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


##########################################################################################
# presence alone
##########################################################################################

# How often is each species present on a spot alone?

# create a new column for the data frame showing when only one species is present
present <- rowSums(dbb[4:15])
dbb <- cbind(dbb,present)
head(dbb)

# get all rows where only one species is present
sub <- subset(dbb,present==1)
nrow(sub)  # 2422

alone <- rep(0,12)
for (i in 1:12) { alone[i] <- sum(sub[ss[i]]) }

# divide by the total number of times the species is present
total <- rep(0,12)
for (i in 1:12) { total[i] <- sum(dbb[ss[i]]) }

# put into table and order by percentage
percent <- alone/total
d <- data.frame(ss,alone,total,percent)
d[rev(order(percent)),]

"  ss alone total    percent
10 WB   941  2030 0.46354680
11 WE   433  1903 0.22753547
4  DH   171   861 0.19860627
12 YC    73   738 0.09891599
7  RB   201  2111 0.09521554
5  ME   158  2376 0.06649832
2  BY   144  2185 0.06590389
3  CF   131  2333 0.05615088
9  SC    68  1464 0.04644809
6  OC    59  1740 0.03390805
8  RG    14   491 0.02851324
1  BH    29  1776 0.01632883"


###################################################################################
# END
###################################################################################

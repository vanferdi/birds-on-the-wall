# this code shows how to load two main data sets for this project
# and helps you check whether it loaded correctly

# load common functions
source("./Code/function_library.R")

##########################################################################################
# load main dataset
##########################################################################################

df <- read.csv("./Data/evendays_dataset.csv", header = TRUE, sep = ",")
nrow(df)  # 7362 rows

# remove unstandardized spots (removes 5 entries)
df <- filter_spots(df,spots)
nrow(df)  # 7357

# check the spots
table(df$Spot)
" 1A   1B   1C   2A   2B   2C   3A   3B  3B1  3B2   3C 
2339  582   15  788  222 1962  135   50  108  138 1018"

# check the dates
sort(df$Rdate)[1]       # first date: 2002-10-02
tail(sort(df$Rdate),1)  # last date:  2012-12-20

# this is what the first 3 rows of data should look like:
head(df,3)
"       Rdate Xhr.dectime Spot BHPA BYMA CFMA DHPA MEPA OCPA RBMA RGMA SCMA WBPA WEPA YCPA
1 2002-10-02    5.500000   3B    0    6   28    0    0    0    8    0    0    0    0    8
2 2002-10-02    5.583333   3A    1    0    1    0    0    0    9    0    0    0    0    4
3 2002-10-02    5.583333   3B    0   22    1    0    0    0    0    0   10    0    0    8"


##########################################################################################
# load validation dataset
##########################################################################################

# only use this data for testing well-defined hypotheses that you developed 
# in response to observing the results of your main analyses on the main data set

dfv <- read.csv("./Data/odddays_dataset.csv", header = TRUE, sep = ",")
nrow(dfv)  # 7682

# remove unstandardized spots (removes 9 entries)
dfv <- filter_spots(dfv,spots)
nrow(dfv)  # 7673

# check the dates
sort(dfv$Rdate)[1]       # first date: 2002-10-01
tail(sort(dfv$Rdate),1)  # last date:  2012-12-19

# this is what the first 3 rows of data should look like:
head(dfv,3)
"       Rdate Xhr.dectime Spot BHPA BYMA CFMA DHPA MEPA OCPA RBMA RGMA SCMA WBPA WEPA YCPA
1 2002-10-01    6.283333   2C   52   13    7    0   13    0   10    0    0    0    0    0
2 2002-10-01    6.283333   3A   11    0    0    0    1   10    0    0    0    0    0    3
3 2002-10-01    6.366667   2C   80   19    7    0   29    0    8    0    0    0    0    2"


##########################################################################################
# END
##########################################################################################
##########################################################################################
# overview
##########################################################################################

# Test if tie type is predicted by the difference in body size between two species
# creates a barplot for visualizing results
# and runs a Mann Whitney U test

# input: 
# copresence_network_ties.RData

# output:
# massdiff_barplot.pdf

# load common functions
source("./Code/function_library.R")


###################################################################################
# load copresence ties matrix from Rdata
###################################################################################

load("./Rdata/copresence_network_ties.RData")

ties
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

# these are the copresence ties that were significant under correction for multiple comparisons
# using the false discovery rate BH method
# see "copresence_network.R" for the calculation

sum(ties==1)  # 45 positive ties
sum(ties==-1) # 47 negative ties


##########################################################################################
# construct matrix showing the body mass difference of every species pair

# pull in the species IDs and mass data
mass <- data.frame(ss,grams)
mass

# matrix of 144 body size differences, for each pair of 12 species
M <- matrix(grams,12,12)
M <- abs(M-t(M))
rownames(M) <- ss
colnames(M) <- ss
M

# all size differences
vals <- M[lower.tri(M)]

# check that all body size differences are unique
length(vals) == length(unique(vals))  # TRUE, good that makes things easier

# rank them
mass_ranked <- rev(sort(vals))

# plot them
barplot(mass_ranked,las=1)


###################################################################################
# get mass differences between species for each affiliative and avoidant tie
###################################################################################

pos_ties <- edgelist_from_matrix(matrix=ties,tie_type="positive")
neg_ties <- edgelist_from_matrix(matrix=ties,tie_type="negative")

pos_ties
" 
[1]  2  1  3  1  5  1  6  1  7  1 11  1 12  1  1  2  3  2  5  2  6  2  8  2  9  2 12  2  1
[30]  3  2  3  5  3  7  3 12  3  6  4  1  5  2  5  3  5  6  5  7  5  9  5  1  6  2  6  4  6
[59]  5  6 11  6  1  7  3  7  5  7 12  7  2  8  9  8  2  9  5  9  8  9  1 11  6 11  1 12  2
[88] 12  3 12"

neg_ties
"
[1]  4  1  4  2 10  2 11  2  8  3 10  3  1  4  2  4  5  4  8  4  9  4 10  4 12  4  4  5 10
[30]  5 12  5  4  7  8  7 10  7  1  8  3  8  4  8  7  8 10  8 11  8  4  9  7  9 10  9 11  9
[59]  1 10  2 10  3 10  4 10  5 10  7 10  8 10  9 10 11 10 12 10  2 11  8 11  9 11 10 11 12
[88] 11  4 12  5 12 10 12
"

length(pos_ties)/2  # 45 positive ties
length(neg_ties)/2  # 47 negative ties

# get the value from M of each positive tie
pos_tie_massdiff <- c()
odds <- seq(1,length(pos_ties),2) 
for (i in odds) {
	x <- M[pos_ties[i],pos_ties[i+1]]
	pos_tie_massdiff <- c(pos_tie_massdiff,x)
}
pos_tie_massdiff
length(pos_tie_massdiff) # 45
mean(pos_tie_massdiff)  # 300.8444

# get the value from M of each negative tie
neg_tie_massdiff <- c()
odds <- seq(1,length(neg_ties),2)
for (i in odds) {
	x <- M[neg_ties[i],neg_ties[i+1]]
	neg_tie_massdiff <- c(neg_tie_massdiff,x)
}
neg_tie_massdiff
length(neg_tie_massdiff) # 47
mean(neg_tie_massdiff)  # 607.6383


##########################################################################################
# now color in the bars that had affiliative ties blue and avoidant ties red
##########################################################################################

# if a bar had at least one positive or negative tie, color it by that tie
# (where bars had two tie type, there was never one positive and one negative)

pos_indexes <- which(mass_ranked %in% pos_tie_massdiff)
neg_indexes <- which(mass_ranked %in% neg_tie_massdiff)

length(pos_indexes)  # 23 pairs with at least one positive tie
length(neg_indexes)  # 26 pairs with at least one negative tie

# check that the sets of pos and neg ties do not overlap
# these were all directed ties, so if they don't overlap that means no pair had one positive and one negative tie
unique(neg_indexes %in% pos_indexes)  # FALSE, they don't overlap
unique(pos_indexes %in% neg_indexes)  # FALSE, they don't overlap

coloring <- rep("grey",length(mass_ranked))
coloring[pos_indexes] <- myblue
coloring[neg_indexes] <- myred
coloring

# plot bar with coloring
barplot(mass_ranked,las=1,col=coloring, ylab="body mass difference",xlab="all species pairs")
legend(50,1000,legend=c("affiliative ties","avoidant ties"),fill=c(myblue,myred))

# save barplot to file
pdf(file = "./Plots/massdiff_barplot.pdf", width = 5.2, height = 2.5)
par(mfrow=c(1,1), mar=c(2,5,1,0))  # (bottom,left,top,right)
barplot(mass_ranked,las=1,col=coloring, ylab="body mass difference \n(grams)",xlab="")
mtext("all species pairs", side = 1, line = .5)
#legend(50,1000,legend=c("affiliative ties","avoidant ties"),fill=c(myblue,myred),bty = "n")
dev.off()


##########################################################################################
# perform a Mann Whitney U test
##########################################################################################

# Does the difference in body mass between two species predict whether their tie type is affiliative or avoidant?
# Hypothesis: mass differences on the negative ties are greater than those on the positive ties (one-tailed test)

# null hypothesis: there is NO difference between the ranks of each condition
# alternative hypothesis: there IS a difference the ranks of each condition

length(pos_tie_massdiff)  # 45 directed positive ties
length(neg_tie_massdiff)  # 47 directed negative ties

Utest_dataframe <- function(pos_set,neg_set) {
	n1 <- length(pos_set)
	n2 <- length(neg_set)
	type <- c(rep("positive",n1),rep("negative",n2))
	value <- c(pos_set,neg_set)
	rank <- seq(1,n1+n2)
	d <- data.frame(type,value)
	d <- d[order(value),]
	d <- cbind(d,rank)
	return(d)
}

d <- Utest_dataframe(pos_tie_massdiff,neg_tie_massdiff)
nrow(d)  # 92 directed ties in total

wilcox.test(value ~ type, data=d, alternative=c("greater"))
# W = 1538, p-value = 0.00008849

# so we reject the null hypotheses 

# NOTE:
# the directed ties above (in d) contain a lot of pairs of identical values, which is ok for the wilcox.test
# it means wilcoxon.test() can't compute an exact p-value and has to do a normal approximation
# alternatively, we can run it on data converted to undirected ties as visualized in the colored barplot above
# this is a worse case scenario, so if it still comes out significant, then I'll be convinced of the result
pos_massdiff_uniques <- unique(pos_tie_massdiff)
neg_massdiff_uniques <- unique(neg_tie_massdiff)
length(pos_massdiff_uniques)  # = the 23 blue bars in the plot
length(neg_massdiff_uniques)  # = the 26 red bars in the plot

d <- Utest_dataframe(pos_massdiff_uniques,neg_massdiff_uniques)
nrow(d)  # 49 undirected ties in total

wilcox.test(value ~ type, data=d, alternative=c("greater"))
# W = 433, p-value = 0.003342

# we still reject the null hypotheses


##########################################################################################
# END
##########################################################################################

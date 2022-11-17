##########################################################################################
# overview
##########################################################################################

# input:
# ncats.RData
# copresence_network_ties.RData

# output:
# complexity_centrality.pdf

source("./Code/function_library.R")  # load common functions
options(scipen=999)

##########################################################################################
# load data
##########################################################################################

# load the number of significant categories that each species uses from the best-fit categorization system
# they are in order of ss (alphabetical order of species name)
load("./Rdata/ncats.RData")

# load the copresence network
load("./Rdata/copresence_network_ties.RData")


###################################################################################
# compute each species' centrality in the affiliative copresence network
###################################################################################

require(igraph)

# create igraph object
pos_ties <- edgelist_from_matrix(ties,tie_type="positive")
g <- make_empty_graph(n = 12, directed=T) + edges(pos_ties)

# degree centrality
# total in and out degree of each node in the network
degree <- round(degree(g, v = V(g), mode = c("all"), loops = TRUE, normalized = FALSE),5)

# eigen vector centrality (should highly correlate with degree)
# high scores mean the species is also connected to others with high scores
# it's related to the stationary distribution of a random walk on the network
eigen <- round(evcent(g)$vector,2)

# hub score (very similar to eigen vector centrality)
hub <- round(hub.score(g)$vector,2)

# betweenness 
# the number of geodesics (shortest paths) going through a vertex, for detecting brokers
between <- round(betweenness(g),2)


###################################################################################
# relationship between complexity and centrality
###################################################################################

d <- data.frame(ss,ncats,degree,eigen,hub,between)
"  ss ncats degree eigen  hub between
1  BH     4     14  1.00 1.00   17.45
2  BY     3     14  0.96 0.94   25.47
3  CF     2     10  0.83 0.83    2.12
4  DH     1      2  0.14 0.14    0.00
5  ME     4     12  0.91 0.93   11.82
6  OC     3     10  0.68 0.67   22.00
7  RB     2      7  0.62 0.55    0.25
8  RG     2      4  0.28 0.28    0.00
9  SC     2      6  0.43 0.43    1.45
10 WB     1      0  0.00 0.00    0.00
11 WE     1      4  0.34 0.34    0.00
12 YC     2      7  0.63 0.71    0.45"

# example:
# BH degree of 14 is for 7 ingoing ties + 7 outgoing ties

# choose degree and betweenness to investigate
# 1) easy to explain the meaning of these and visually understand on
# 2) they show two different aspects of centrality and aren't uber correlated like the other measures are
cor(degree,between)  # 0.7967577
cor(degree,eigen)    # 0.9830744


###################################################################################
# relationship between complexity and centrality
###################################################################################

summary(lm(d$ncats ~ d$degree))
"            Estimate Std. Error t value  Pr(>|t|)    
(Intercept)  0.70779    0.28486   2.485    0.0323 *  
d$degree     0.20563    0.03278   6.272 0.0000923 ***"

# interpretation: one more degree buys you 0.20563 more category systems

summary(lm(d$ncats ~ d$between))
"            Estimate Std. Error t value Pr(>|t|)
(Intercept)   1.6917     0.2562   6.603 6.06e-05
d$between     0.0827     0.0223   3.708  0.00405"

# interpretation: being on one more betweenness path buys you 0.0827 more category systems


###################################################################################
# plot
###################################################################################

require(ggplot2)

n <- 0.2
nudges <- c(rep(n,11),-n)  # because YC is on top of RB
pdf(file = "./Plots/complexity_degree.pdf", width = 5, height = 3)
par(mfrow=c(1,1), mar=c(1,1,1,1))
ggplot(d,aes(degree, ncats)) +
    geom_smooth(method='lm', size = 0.5, color="#63666A") +
    geom_text(aes(label=ss), position = position_nudge(y = nudges)) +
    geom_point() +
    xlab("degree centrality") + ylab("number of categories")
dev.off()

# these points are really jammed up together, fix manually
labs <- c("BH","BY","","","ME","OC","","","","","","")
nudges <- c(rep(n,11),-n)
pdf(file = "./Plots/complexity_between.pdf", width = 5, height = 3)
par(mfrow=c(1,1), mar=c(1,1,1,1))
ggplot(d,aes(between, ncats)) +
    geom_smooth(method='lm', size = 0.5, color="#63666A") +
    geom_text(aes(label=labs), position = position_nudge(y = nudges)) +
    geom_point() +
    # add in a row from lowest to highest betweenness
    annotate("text", x=6.5, y=2, label= "RG RB YC SC CF") +
    annotate("text", x=3, y=1, label= "DH WB WE") +
    xlab("betweenness centrality") + ylab("number of categories")
dev.off()

# By default, the envelope is the 95% confidence level interval for predictions from a linear model ("lm")


###################################################################################
# END
###################################################################################

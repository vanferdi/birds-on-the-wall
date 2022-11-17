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
1  BH     2     14  1.00 0.99   21.02
2  BY     3     13  0.94 1.00   18.15
3  CF     1     12  0.96 0.93    4.78
4  DH     3      2  0.11 0.13    0.00
5  ME     3     12  0.94 0.93   12.73
6  OC     7      9  0.55 0.44   20.65
7  RB     1      8  0.71 0.69    0.50
8  RG     1      4  0.30 0.29    0.00
9  SC     2      8  0.62 0.59    3.33
10 WB     1      0  0.00 0.00    0.00
11 WE     1      4  0.30 0.32    0.00
12 YC     2      8  0.71 0.68    0.83"


###################################################################################
# relationship between complexity and centrality
###################################################################################

summary(lm(d$ncats ~ d$degree))
"            Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.53840    1.03899   1.481    0.169
d$degree     0.09084    0.11604   0.783    0.452"

# interpretation: one more degree buys you 0.09084 more category systems - not significant

summary(lm(d$ncats ~ d$between))
"            Estimate Std. Error t value Pr(>|t|)  
(Intercept)  1.36277    0.50320   2.708   0.0220 *
d$between    0.12985    0.04668   2.782   0.0194 *"

# interpretation: being on one more betweenness paths buys you 0.12985 more category systems


###################################################################################
# plot
###################################################################################

require(ggplot2)

n <- 0.4
nudges <- c(rep(n,6),-n,n,-n,n,-n,n)  
pdf(file = "./Plots/complexity_degree.pdf", width = 5, height = 3)
par(mfrow=c(1,1), mar=c(1,1,1,1))
ggplot(d,aes(degree, ncats)) +
    geom_smooth(method='lm', size = 0.5, color="#63666A") +
    geom_text(aes(label=ss), position = position_nudge(y = nudges)) +
    geom_point() +
    xlab("degree centrality") + ylab("")
dev.off()

# these points are really jammed up together, fix manually
labs <- c("BH","BY","CF","DH","ME","OC","","","SC","","","YC")
pdf(file = "./Plots/complexity_between.pdf", width = 5, height = 3)
par(mfrow=c(1,1), mar=c(1,1,1,1))
ggplot(d,aes(between, ncats)) +
    geom_smooth(method='lm', size = 0.5, color="#63666A") +
    geom_text(aes(label=labs), position = position_nudge(y = n)) +
    geom_point() +
    # add in a row from lowest to highest betweenness
    annotate("text", x=0.3, y=1+n, label= "RG RB") +
    annotate("text", x=0.4, y=1-n, label= "WB WE") +
    xlab("betweenness centrality") + ylab("")
dev.off()

# By default, the envelope is the 95% confidence level interval for predictions from a linear model ("lm")


###################################################################################
# END
###################################################################################

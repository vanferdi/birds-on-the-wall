
require(igraph)

circle_layout <- matrix(0,12,2)
circle_layout[ssi("ME"),] <- c(1.000, 0.000) # this is the far right node
circle_layout[ssi("SC"),] <- c(0.866, 0.500) # then it goes up and around counterclockwise
circle_layout[ssi("BY"),] <- c(0.500, 0.866)
circle_layout[ssi("RG"),] <- c(0.000, 1.000)
circle_layout[ssi("DH"),] <- c(-0.500, 0.866)
circle_layout[ssi("WB"),] <- c(-0.866, 0.500)
circle_layout[ssi("WE"),] <- c(-1.000, 0.000)
circle_layout[ssi("OC"),] <- c(-0.866, -0.500)
circle_layout[ssi("BH"),] <- c(-0.500, -0.866)
circle_layout[ssi("RB"),] <- c(0.000, -1.000)
circle_layout[ssi("CF"),] <- c(0.500, -0.866)
circle_layout[ssi("YC"),] <- c(0.866, -0.500)

save(circle_layout, file = "./Rdata/circle_layout.RData")

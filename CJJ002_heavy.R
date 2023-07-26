# down no use
library(dplyr)
library(tidyverse)
DN_heavy <- read.delim("e:/P23042711/DN-B/heavy_parse-select.tsv", 
                       stringsAsFactors = FALSE)
DP_heavy <- read.delim("e:/P23042711/DP-B/heavy_parse-select.tsv", 
                       stringsAsFactors = FALSE)
NB_heavy <- read.delim("e:/P23042711/NB-B/heavy_parse-select.tsv", 
                       stringsAsFactors = FALSE)
PC_heavy <- read.delim("e:/P23042711/PC-B/heavy_parse-select.tsv", 
                       stringsAsFactors = FALSE)
library(alakazam)
library(dowser)
library(dplyr)
library(scales)
library(shazam)
DN_heavy_dist_ham <- distToNearest(DN_heavy , 
                                   cellIdColumn="cell_id", locusColumn="locus", 
                                   VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)

library(ggplot2)
p1 <- ggplot(subset(DN_heavy_dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)

output <- findThreshold(DN_heavy_dist_ham$dist_nearest, method="density")
threshold <- output@threshold
plot(output, title="Density Method")
print(output)


DP_heavy_dist_ham <- distToNearest(DP_heavy , 
                                   cellIdColumn="cell_id", locusColumn="locus", 
                                   VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)

library(ggplot2)
p1 <- ggplot(subset(DP_heavy_dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)

output <- findThreshold(DP_heavy_dist_ham$dist_nearest, method="density")
threshold <- output@threshold
plot(output, title="Density Method")
print(output)


NB_heavy_dist_ham <- distToNearest(NB_heavy , 
                                   cellIdColumn="cell_id", locusColumn="locus", 
                                   VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)

library(ggplot2)
p1 <- ggplot(subset(NB_heavy_dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)

output <- findThreshold(NB_heavy_dist_ham$dist_nearest, method="density")
threshold <- output@threshold
plot(output, title="Density Method")
print(output)


PC_heavy_dist_ham <- distToNearest(PC_heavy , 
                                   cellIdColumn="cell_id", locusColumn="locus", 
                                   VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)

p1 <- ggplot(subset(PC_heavy_dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)

output <- findThreshold(PC_heavy_dist_ham$dist_nearest, method="density")
threshold <- output@threshold
plot(output, title="Density Method")
print(output)
# up no use

###
# library(dplyr)
# library(tidyverse)
# DN_heavy_new <- read.delim("e:/P23042711/DN-B/20230706/Result/05.match/heavy_parse-select.tsv", 
#                        stringsAsFactors = FALSE)
# DP_heavy_new <- read.delim("e:/P23042711/DP-B/20230706/Result/05.match/heavy_parse-select.tsv", 
#                        stringsAsFactors = FALSE)
# NB_heavy_new <- read.delim("e:/P23042711/NB-B/20230706/Result/05.match/heavy_parse-select.tsv", 
#                        stringsAsFactors = FALSE)
# PC_heavy_new <- read.delim("e:/P23042711/PC-B/20230706/Result/05.match/heavy_parse-select.tsv", 
#                        stringsAsFactors = FALSE)
# library(alakazam)
# library(dowser)
# library(dplyr)
# library(scales)
# library(shazam)
# DN_heavy_new_dist_ham <- distToNearest(DN_heavy_new , 
#                                        cellIdColumn="cell_id", locusColumn="locus", 
#                                        VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)
# 
# p1 <- ggplot(subset(DN_heavy_new_dist_ham, !is.na(dist_nearest)),
#              aes(x=dist_nearest)) + 
#   theme_bw() + 
#   xlab("Hamming distance") + 
#   ylab("Count") +
#   scale_x_continuous(breaks=seq(0, 1, 0.1)) +
#   geom_histogram(color="white", binwidth=0.02) +
#   geom_vline(xintercept=0.12, color="firebrick", linetype=2)
# plot(p1)
# 
# output <- findThreshold(DN_heavy_new_dist_ham$dist_nearest, method="density")
# threshold <- output@threshold
# plot(output, title="Density Method")
# print(output)
# 
# 
# DP_heavy_new_dist_ham <- distToNearest(DP_heavy_new , 
#                                        cellIdColumn="cell_id", locusColumn="locus", 
#                                        VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)
# 
# p1 <- ggplot(subset(DP_heavy_new_dist_ham, !is.na(dist_nearest)),
#              aes(x=dist_nearest)) + 
#   theme_bw() + 
#   xlab("Hamming distance") + 
#   ylab("Count") +
#   scale_x_continuous(breaks=seq(0, 1, 0.1)) +
#   geom_histogram(color="white", binwidth=0.02) +
#   geom_vline(xintercept=0.12, color="firebrick", linetype=2)
# plot(p1)
# 
# output <- findThreshold(DP_heavy_new_dist_ham$dist_nearest, method="density")
# threshold <- output@threshold
# plot(output, title="Density Method")
# print(output)
# 
# 
# NB_heavy_new_dist_ham <- distToNearest(NB_heavy_new , 
#                                        cellIdColumn="cell_id", locusColumn="locus", 
#                                        VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)
# 
# p1 <- ggplot(subset(NB_heavy_new_dist_ham, !is.na(dist_nearest)),
#              aes(x=dist_nearest)) + 
#   theme_bw() + 
#   xlab("Hamming distance") + 
#   ylab("Count") +
#   scale_x_continuous(breaks=seq(0, 1, 0.1)) +
#   geom_histogram(color="white", binwidth=0.02) +
#   geom_vline(xintercept=0.12, color="firebrick", linetype=2)
# plot(p1)
# 
# output <- findThreshold(NB_heavy_new_dist_ham$dist_nearest, method="density")
# threshold <- output@threshold
# plot(output, title="Density Method")
# print(output)
# 
# 
# PC_heavy_new_dist_ham <- distToNearest(PC_heavy_new , 
#                                        cellIdColumn="cell_id", locusColumn="locus", 
#                                        VJthenLen=FALSE, onlyHeavy=FALSE, nproc=40)
# 
# p1 <- ggplot(subset(PC_heavy_new_dist_ham, !is.na(dist_nearest)),
#              aes(x=dist_nearest)) + 
#   theme_bw() + 
#   xlab("Hamming distance") + 
#   ylab("Count") +
#   scale_x_continuous(breaks=seq(0, 1, 0.1)) +
#   geom_histogram(color="white", binwidth=0.02) +
#   geom_vline(xintercept=0.12, color="firebrick", linetype=2)
# plot(p1)
# 
# output <- findThreshold(PC_heavy_new_dist_ham$dist_nearest, method="density")
# threshold <- output@threshold
# plot(output, title="Density Method")
# print(output)


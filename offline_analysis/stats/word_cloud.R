
## install packages using the following
install.packages("devtools") # to install from github
library(devtools)
install_github("ropensci/rorcid") # to access orcid through R
install.packages(c("tm","wordcloud")) # to draw word clouds

## make a word cloud
make.cloud(titles)



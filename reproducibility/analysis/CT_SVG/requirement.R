
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("dmcable/spacexr")

library(Matrix)
library(doParallel)

devtools::install_github("pekjoonwu/CELINA")


devtools::install_git("https://github.com/FIRST-Data-Lab/BPST.git")
devtools::install_github("shanyu-stat/spVC")
devtools::install_github("FIRST-Data-Lab/Triangulation")

devtools::install_github("jingeyu/CTSV")


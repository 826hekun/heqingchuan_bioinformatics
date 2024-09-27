rm(list = ls())
rdatas = dir(pattern = ".Rdata$",recursive = T)
rdatas = rdatas[-c(9,10,17,18)]
#file.remove(rdatas)

rmds = dir(pattern = ".Rmd$",recursive = T)
rmds
for(i in 15:22){
  #i = 20
  rmds = dir(pattern = ".Rmd$",recursive = T)
  rmarkdown::render(rmds[[i]])
  rmds = dir(pattern = ".Rmd$",recursive = T)
}

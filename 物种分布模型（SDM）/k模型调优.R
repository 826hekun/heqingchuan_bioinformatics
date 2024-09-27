#工作路径E:/2023maxentrain/data/Model/

if(!require(devtools)){
  install.packages("devtools")
}
if(!require(kuenm)){
  devtools::install_github("marlonecobos/kuenm")
}

library(kuenm)
setwd("E:/2023maxentrain/data/Model")
# arguments
occs <- read.csv("E:/2023maxentrain/data/Model/abc.csv")
occs <- data.frame(Species = "abc", occs)
train_prop <- 0.75
method = "random"

#running，划分数据
data_split <- kuenm_occsplit(occ = occs, train.proportion = train_prop, method = method, save = TRUE, name = "abc")

### preparing arguments(change "YOUR/DIRECTORY" by your pertinent directory)
occ_joint <- "abc_joint.csv"
occ_tra <- "abc_train.csv"
M_var_dir <- "Environ_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_models"
reg_mult <- seq(0.1, 4, 0.1)
f_clas <- "all" #总共五个
args <- NULL
maxent_path <- "E:/2023maxentrain/data/Model/" #maxent所在位置，最好拷贝到工作目录下
wait <- FALSE
run <- TRUE


## running candidate models
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal, out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas,
          args = args, maxent.path = maxent_path, wait = wait, run = run)

#Candidate model evaluation and selection
### check the functions help to understand arguments

## preparing arguments (change "YOUR/DIRECTORY" by your pertinent directory)
occ_test <- "abc_test.csv"
out_eval <- "Calibration_Results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
parallel_proc <- FALSE  #10

cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra,
                        occ.test = occ_test, batch = batch_cal, out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations, kept = kept, selection = selection, parallel.proc = parallel_proc)



# 设置本地工作目录为输出文件夹
output_path <- "/chenlab/heyuanyu/MegaSDMkun/"
setwd(output_path)

# 开始日志记录
log_file <- file("/chenlab/heyuanyu/Rscripts/megaSDM_kun_log.txt", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

output_path <- "/chenlab/heyuanyu/MegaSDMkun/"
setwd(output_path)
 library(megaSDM)
 input_TA <- list.files(
   "/chenlab/heyuanyu/BioClimTIF/BioClim_CN/", # 使用 system.file 函数来获取 megaSDM 包中 extdata 文件夹下的 trainingarea 子文件夹的路径。
   pattern = ".tif$", # 搜索以 .tif 结尾的文件。
   full.names = TRUE # 返回完整的文件路径名。
   )
 # If you have your own data, replace the system.file command with
 # a pathway to the training area files.

# -----------------------------------------------------------------------------
 envoutput <- "TestRun"

# -----------------------------------------------------------------------------
 # Here we define the extent of the training and study regions in c(xmin, xmax, ymin, ymax) form.
 TSEnv <- TrainStudyEnv(input_TA = input_TA,
                        output = envoutput,
                        clipTrain = NA,
                        clipStudy = NA)

# -----------------------------------------------------------------------------
 Env2050_ssp245 <- list.files("/chenlab/heyuanyu/BioClimTIF/BCC_future_CN/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2041-2060/",
                           pattern = ".tif$",
                           full.names = TRUE)
 Env2070_ssp245 <- list.files("/chenlab/heyuanyu/BioClimTIF/BCC_future_CN/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2061-2080/",
                           pattern = ".tif$",
                           full.names = TRUE)
 Envssp245 <- list(Env2050_ssp245, Env2070_ssp245)
 
 # The "time_periods" argument must contain the current time (the time of the training
 # and study rasters) first and then the time periods for the forecast/hindcast.
 
 PredictEnv(studylayers = TSEnv$study,
            futurelayers = Envssp245,
            time_periods = c(2010, 2050, 2070),
            output = envoutput,
            scenario_name = "ssp245")
 
 #Repeat with a different climate scenario (ssp585):
 Env2050_ssp585 <- list.files("/chenlab/heyuanyu/BioClimTIF/BCC_future_CN/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp585_2041-2060/",
                           pattern = ".tif$",
                           full.names = TRUE)
 Env2070_ssp585 <- list.files("/chenlab/heyuanyu/BioClimTIF/BCC_future_CN/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp585_2061-2080/",
                           pattern = ".tif$",
                           full.names = TRUE)
 
 Envssp585 <- list(Env2050_ssp585, Env2070_ssp585)
 PredictEnv(studylayers = TSEnv$study,
            futurelayers = Envssp585,
            time_periods = c(2010, 2050, 2070),
            output = envoutput,
            scenario_name = "ssp585")
 occ_output <- "occurrences"
# Occurrences <- OccurrenceCollection(spplist = spplist,
#                                     output = occ_output,
#                                     trainingarea = extent_occ)
 
 # 注意：当使用R Markdown运行此操作时，可能会出现“不完整的最终行...”
 # 警告。然而，它们似乎不会影响出现点的总数或身份，
 # 当在控制台上运行代码时，警告不会出现
# 因为一个物种被重命名，重命名物种列表以反映分类学变化
# spplist <- Occurrences$Scientific.Name
 spplist <- "abc_CN"
 
 #使用R合并下载数据为一个物种文件abc_CN.csv
 # -----------------------------------------------------------------------------
 # 首先，获取出现文件列表
 occlist <- list.files(occ_output, pattern = "abc_CN.csv", full.names = TRUE)
 
 OccurrenceManagement(occlist = occlist,
                      output = occ_output,
                      envextract = TRUE,
                      envsample = TRUE,
                      nbins = 25,
                      envdata = TSEnv$training)
 
 # -----------------------------------------------------------------------------
 # 再次获取出现文件列表，即使它们以前被写在同一文件夹中。这确保了出现文件
 # 正确格式化。
 occlist <- list.files(occ_output, pattern = "abc_CN.csv", full.names = TRUE)
 # 打印背景缓冲区(.shp)的位置（如果不存在则创建）
 buff_output <- "TestRun/buffers"
 # 为每个物种生成缓冲区。
 BackgroundBuffers(occlist = occlist,
                   envdata = TSEnv$training,
                   buff_output,
                   ncores = 2)
 # -----------------------------------------------------------------------------
 # 设置背景点生成的参数
 # （多少个点，以及空间约束程度）
 # 应为每个物种生成多少个背景点？
 nbg <- 1000
 # 多少比例的背景点应该从缓冲区内采样？
 spatial_weights <- 0.5
 # 背景点应该是环境子采样（Varela）还是
 # 随机分布（random）？
 sampleMethod <- "Varela"
 # 因为我们希望有部分空间约束（50%的点在缓冲区内），我们必须制作
 # 用于生成背景点的缓冲文件列表。在此示例中，
 # 这些文件是由BackgroundBuffers函数创建的，但也可以在megaSDM之外生成并带入这里。
 bufflist <- list.files(buff_output, pattern = ".shp$", full.names = TRUE)
 # 定义背景点将打印到的位置(.shp)
 # （如果不存在，则创建此目录）
 bg_output <- "TestRun/backgrounds"
 BackgroundPoints(spplist = spplist,
                  envdata = TSEnv$training,
                  output = bg_output,
                  nbg = nbg,
                  spatial_weights = spatial_weights,
                  buffers = bufflist,
                  method = sampleMethod,
                  ncores = 2)
 # -----------------------------------------------------------------------------
 # 定义每个物种要保留的环境变量列表
 # 在这个例子中，我们希望所有物种都有相同的环境变量。
 envvar <- rep("hfp,elev,bio_09,slope,bio_14,ocd,bio_02,aspect,phh2o,nitro,bio_03,bio_15,bio_13", length = length(occlist))
 # 定义背景点文件列表
 # （由BackgroundPoints函数创建或单独生成）
 bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)
 # 在这个例子中，megaSDM覆盖了出现点和背景点，
 # 如果需要，可以放在不同的文件夹中。
 VariableEnv(occlist = occlist,
             bglist = bglist,
             env_vars = envvar,
             occ_output = occ_output,
             bg_output = bg_output)
 # -----------------------------------------------------------------------------
 # 首先，定义所有背景和出现点文件列表
 occlist <- list.files(occ_output, pattern = "abc_CN.csv", full.names = TRUE)
 bglist <- list.files(bg_output, pattern = ".csv", full.names = TRUE)
 # 定义MaxEnt模型运行结果将打印到的位置（作为.lambdas文件）
model_output <- "TestRun/models"
# “nrep”设置为4，意味着MaxEnt算法将运行4次，使用不同的
# 出现点子集以更好地表示栖息地适宜性。
 MaxEntModel( # 运行 MaxEnt 模型。
   occlist = occlist, # 发生数据文件列表。
   bglist = bglist, # 背景点文件列表。
   model_output = model_output, # MaxEnt 模型结果输出文件夹的名称。
   ncores = 2, # 并行处理的核心数量。
   nrep = 4, # MaxEnt 算法运行次数。
   alloutputs = TRUE, # 是否输出所有中间结果。
   reptype = "Subsample",#交叉验证的类型（"Crossvalidate", "Bootstrap", "Subsample"；
   test_percent = 20,#整数，介于 0 和 100 之间：用于交叉验证的点的百分比（例如，测试 AUC 验证等）。默认值为 20。
   features = c("linear", "quadratic", "product", "threshold", "hinge"),#可选。MaxEnt 用于建模物种-环境关系的特征向量。可选特征包括 "linear", "quadratic", "product", "threshold", "hinge"。如果发生点较少，不建议使用 "hinge" 特征。默认值为所有特征类。
   testsamples = FALSE,
   regularization = 1#数值。正则化参数（惩罚复杂模型）。更高的正则化意味着更倾向于简单模型。默认值为 1。
 )
 # -----------------------------------------------------------------------------
 # 首先，创建一个分析中使用的时期和气候情景列表
 # （从模型训练的年份开始）
 time_periods <- c(2010,2050,2070)
 scenarios <- c("ssp245", "ssp585")
 # 定义当前研究区域栅格所在的目录
 # （由TrainStudyEnv函数生成或从单独位置带入）
 study_dir <- "TestRun/studyarea"
 # 定义未来研究区域栅格所在的目录
 # （由PredictEnv函数生成或从单独位置带入）
 # 定义预测气候层的目录列表，
 # 按不同气候情景和年份分开：
 # list(c(Scenario1Year1, Scenario1Year2),
 #      c(Scenario2Year1, Scenario2Year2))
 predictdir <- list(c("TestRun/ssp245/2050",
                      "TestRun/ssp245/2070"),
                    c("TestRun/ssp585/2050",
                      "TestRun/ssp585/2070"))
 # 定义结果将打印到的位置。
 # 对于这个例子，我们将定义一个工作目录内的新文件夹，
 # 专门用于模型投影和分析。
 result_dir <- "Results"
 # 其他选项也可用（参阅文档页面）
 MaxEntProj( # 运行 MaxEnt 模型预测。
   input = model_output, # MaxEnt 模型输出文件夹的名称。
   time_periods = time_periods, # 时间周期列表。
   scenarios = scenarios, # 气候情景列表。
   study_dir = study_dir, # 当前研究区域栅格数据所在的目录。
   predict_dirs = predictdir, # 不同气候情景和年份下的预测气候层目录列表。
   output = result_dir, # 模型预测结果输出文件夹的名称。
   aucval = 0.5, # AUC 值阈值。要大于0.7
   ncores = 2 # 并行处理的核心数量。
 )
 # -----------------------------------------------------------------------------
 # 时间地图将被写入在"result_dir"中提供的目录
   result_dir <- "Results"
   createTimeMaps(result_dir = result_dir,
                  time_periods = time_periods,
                  scenarios = scenarios,
                  dispersal = FALSE,
                  ncores = 1)
#-----------------------------------------------------------------------------
   additionalStats(result_dir = result_dir,
                   time_periods = time_periods,
                   scenarios = scenarios,
                   dispersal = FALSE,
                   ncores = 1)
 # -----------------------------------------------------------------------------
   # Add in dispersal data (normally you would read a .csv file with two columns, but in this example
   # the data is just added in by hand here).
   dispersaldata <- data.frame(Species = spplist, Rate = c(6.92))
 
   dispersalRate(result_dir = result_dir,
                 dispersaldata = dispersaldata,
                 time_periods = time_periods,
                 scenarios = scenarios,
                 ncores = 1)
 
   # Repeat the time map and additional stats steps for the dispersal constrained data.
   # Set dispersal = TRUE this time.
   createTimeMaps(result_dir,
                  time_periods,
                  scenarios,
                  dispersal = TRUE,
                  dispersaldata = dispersaldata,
                  ncores = 1)
 
   # The additional stats function will compare the species ranges between the
   # dispersal-constrained and the regular data.
   additionalStats(result_dir,
                   time_periods,
                   scenarios,
                   dispersal = TRUE,
                   dispersaldata = dispersaldata,
                   ncores = 1)

 # -----------------------------------------------------------------------------
   createRichnessMaps(result_dir = result_dir,   #taken from previous steps
                      time_periods = time_periods, #taken from previous steps
                      scenarios = scenarios, #taken from previous steps
                      dispersal = TRUE,
                      taxonlist = FALSE)
   # A list of the higher taxa of each species (e.g., family) can be provided (taxonlist) to create
   # separate richness maps for each higher taxon.

# 结束日志记录
sink(type = "output")
sink(type = "message")
close(log_file)

 
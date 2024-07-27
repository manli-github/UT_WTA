pacman::p_load(tidyverse,linelist,,ggpubr,flextable,readtext,mlogit,ggplot2)

################## Import data and statified sampling ################## 
# Import the data that has been spatially sampled
mydata <- import('~/mydata.csv')

# Create a balanced panel data containing all years and uid
balanced_id <- expand.grid(
  uid = unique(mydata_s$uid),
  Year = unique(mydata_s$Year)
)
balanced_id <- balanced_id %>% arrange(uid,Year)

# Stratified Sampling
set.seed(40)
stratified_id <- balanced_id %>%
  group_by(uid) %>%
  sample_n(1) %>%
  ungroup() 
sample_data <- merge(stratified_id, mydata, by=c('uid','Year'), all.x=T) %>% drop_na()

# Factoring Choice and Choice_lag to facilitate discrete-choice estimation
sample_data$Choice <- factor(sample_data$Choice, 
                             levels = c("Alfalfa", "Fallow", "Hay", "OtherCrop", "OtherGrain", "Wheat"))
sample_data$Choice_lag <- factor(sample_data$Choice_lag, 
                                 levels = c("Alfalfa", "Fallow", "Hay", "OtherCrop", "OtherGrain", "Wheat"))

# Randomly sampling into train and test data
sample_indices <- sample(1:nrow(sample_data), 0.5 * nrow(sample_data))
train_data <- sample_data[sample_indices, ]
test_data <- sample_data[-sample_indices, ] 

################## Estimate discrete choice model using mlogit function ################## 
# function to convert train_data/test_data from the regular structure to a long structure in mlogit format
convert_to_DCE_data <- function(input) {
  train_data1 <- input %>% 
    mutate(Choice_lagAlfalfa = as.numeric(Choice_lag == 'Alfalfa')) %>%
    mutate(Choice_lagFallow = as.numeric(Choice_lag == 'Fallow')) %>% 
    mutate(Choice_lagHay = as.numeric(Choice_lag == 'Hay')) %>%   
    mutate(Choice_lagOtherCrop = as.numeric(Choice_lag == 'OtherCrop')) %>% 
    mutate(Choice_lagOtherGrain = as.numeric(Choice_lag == 'OtherGrain')) %>% 
    mutate(Choice_lagWheat = as.numeric(Choice_lag == 'Wheat')) %>% 
    mutate(basin_BearRiver = as.numeric(basin == 'BEAR RIVER')) %>% 
    mutate(basin_WeberRiver = as.numeric(basin == 'WEBER RIVER')) %>% 
    mutate(basin_JordanRiver = as.numeric(basin == 'JORDAN RIVER')) %>% 
    mutate(basin_NonGSL = as.numeric(basin == 'NON GSL')) %>% 
    dplyr::select(uid,Choice,rent_Fallow,rent_Alfalfa,rent_Hay,rent_OtherCrop,rent_OtherGrain,rent_Wheat,
                  Choice_lagAlfalfa,Choice_lagFallow,Choice_lagHay,Choice_lagOtherCrop,Choice_lagOtherGrain,Choice_lagWheat,
                  basin_BearRiver,basin_WeberRiver,basin_JordanRiver,basin_NonGSL)

  # Transform regular 'wide' data to 'long' data to run mlogit()
  train_data_long <- mlogit.data(train_data1,shape='wide',choice='Choice') 
  
  # Export and import the mlogit data to further process 'rent' variable
  write.csv(train_data_long,'~/tmp.csv',row.names=FALSE)
  train_data_long1 <- import('~/tmp.csv') %>% 
    mutate(rent = case_when(alt == 'Alfalfa' ~ rent_Alfalfa,
                            alt == 'Fallow' ~ rent_Fallow,
                            alt == 'Hay' ~ rent_Hay,
                            alt == 'OtherCrop' ~ rent_OtherCrop,
                            alt == 'OtherGrain' ~ rent_OtherGrain,
                            alt == 'Wheat' ~ rent_Wheat
    ))
  
  # After importing and further processing, the data is no long mlogit.data; re-transform to mlogit data    
  output <- mlogit.data(data = train_data_long1, choice = 'Choice', shape = "long", id.var = c("chid"),alt.var = c("alt"))
  return(output)
}
DCE_data_train <- convert_to_DCE_data(train_data_rent)
DCE_data_test <- convert_to_DCE_data(test_data_rent)

# Estimation
mod <- mlogit(Choice ~ rent | -1 + Choice_lagFallow + Choice_lagHay + Choice_lagOtherCrop + Choice_lagOtherGrain  + Choice_lagWheat
                 + basin_BearRiver + basin_WeberRiver + basin_JordanRiver + basin_NonGSL, data = DCE_data_train, reflevel = 'Fallow')
summary(mod)

# Calculate pseudo McFadden's R2
mod0 <- mlogit(Choice ~ 1, data = DCE_data_train, reflevel = 'Fallow')
get_mcfadden_pseudo_R2 <- function(model, model0, data) {
  mod.loglik <- model$logLik
  mod0.loglik <- model0$logLik
  return(as.numeric(1 - mod.loglik / mod0.loglik))
}
get_mcfadden_pseudo_R2(mod,mod0,DCE_data) %>% round(.,3)

# Obtain elasticity estimate
elas_all <- effects(mod, covariate = 'rent', type = "ar") # Given one percent change in rent, how much the probability would change in absolute value
elas_own <- as.data.frame(diag(elas)) # Own-elasticity (Table S2)

# Bootstrapping
uid <- train_data %>% dplyr::select(uid) %>% arrange(uid) 
matrix_to_vector <- function(matrix){
  vector <- rep(0,nrow(matrix)*ncol(matrix))
  for (i in 1:nrow(matrix)) {
    for (j in 1:ncol(matrix)) {
      index <- (i - 1) * ncol(matrix) + j
      if (index <= length(vector)) {
        vector[index] <- matrix[i, j]
      }
    }
  }
  return(vector)
}
n_bootstrap <- 1000
beta_bootstrap <- matrix(0, nrow = n_bootstrap, ncol = 46)
elas_bootstrap <- matrix(0, nrow = n_bootstrap, ncol = 30)
set.seed(123)
for (b in 1:n_bootstrap) {
  uid_resample <- mosaic::resample(uid,replace = TRUE) %>% arrange(uid) %>% mutate(pid=row_number()) %>% dplyr::select(uid,pid)
  tmp <- merge(uid_resample,DCE_data,by='uid',all.x = TRUE) %>% drop_na(pid) %>% mutate(chid=pid)
  resampled_DCE_data <- mlogit.data(data = tmp, choice = 'Choice', shape = "long", id.var = c("chid"),alt.var = c("alt"))
  mod_boot <- mlogit(Choice ~ rent | -1 + Choice_lagFallow + Choice_lagHay + Choice_lagOtherCrop + Choice_lagOtherGrain  + Choice_lagWheat
                + basin_BearRiver + basin_WeberRiver + basin_JordanRiver + basin_NonGSL, data = resampled_DCE_data, reflevel = 'Fallow')
  beta <- t(mod_boot$coefficients)
  elas <- effects(mod_boot, covariate = 'rent', type = "ar")[-1,] 
  elas_v <- matrix_to_vector(elas)
  # Store the data frame output in the combined data frame
  beta_bootstrap[b, ] <- beta
  elas_bootstrap[b,] <- elas_v
}

beta_bootstrap <- data.frame(beta_bootstrap) 
elas_bootstrap <- data.frame(me_bootstrap) 
write.csv(beta_bootstrap,'~/beta_bootstrap.csv',row.names=FALSE)
write.csv(elas_bootstrap,'~/elas_bootstrap.csv',row.names=FALSE)

elas_own_bootstrap <- elas_bootstrap[,c(2,9,16,23,30)]
elas_own_se <- sapply(elas_own_bootstrap,sd) # bootstrap standard errors for own elasticity (Table S2)

# Access Model performance and create the output table (Table S3)
console_summary_text <- function(model, train, test, test_dce){
  
  # Function to calculate accuracy
  calculate_accuracy <- function(confusion_matrix) {
    sum(diag(confusion_matrix)) / sum(confusion_matrix)
  }
  # Function to calculate precision
  calculate_precision <- function(confusion_matrix) {
    diag(confusion_matrix) / rowSums(confusion_matrix)
  }
  # Function to calculate recall
  calculate_recall <- function(confusion_matrix) {
    diag(confusion_matrix) / colSums(confusion_matrix)
  }
  
  # train_data
  prob_train <- model$probabilities
  pred_train <- colnames(prob_train)[apply(prob_train, 1, which.max)]
  confusion_matrix_train <- table(Actual = train$Choice, Predicted = pred_train)
  
  # Calculate metrics for trained data
  accuracy_train <- calculate_accuracy(confusion_matrix_train)
  precision_train <- calculate_precision(confusion_matrix_train)
  recall_train <- calculate_recall(confusion_matrix_train)
  
  cat("\nconfusion_matrix:\n")
  print(confusion_matrix_train)
  cat("\naccuracy:\n")
  print(accuracy_train)
  cat("\nprecision:\n")
  print(precision_train)
  cat("\nrecall:\n")
  print(recall_train)
  
  # test_data
  prob_test <- predict(model, newdata = test_dce)
  pred_test <- colnames(prob_test)[apply(prob_test, 1, which.max)]
  confusion_matrix_test <- table(Actual = test$Choice, Predicted = pred_test)
  
  # Calculate metrics for trained data
  accuracy_test <- calculate_accuracy(confusion_matrix_test)
  precision_test <- calculate_precision(confusion_matrix_test)
  recall_test <- calculate_recall(confusion_matrix_test)
  
  cat("\nconfusion_matrix:\n")
  print(confusion_matrix_test)
  cat("\naccuracy:\n")
  print(accuracy_test)
  cat("\nprecision:\n")
  print(precision_test)
  cat("\nrecall:\n")
  print(recall_test)
}
console_summary_text(model = mod, train = train_data, test = test_data, test_dce=DCE_data_test)

############ Application: switching alfalfa to less water-intensive use using 2022 data ############
# Import 2022 data for alfalfa (30 m)
alfalfa_2022 <- import('~/alfalfa_2022.csv')
alfalfa_2022 <- alfalfa_2022 %>% 
  mutate(Choice_lagFallow = as.numeric(Choice == 'Fallow')) %>% 
  mutate(Choice_lagHay = as.numeric(Choice == 'Hay')) %>%   
  mutate(Choice_lagOtherCrop = as.numeric(Choice == 'OtherCrop')) %>% 
  mutate(Choice_lagOtherGrain = as.numeric(Choice == 'OtherGrain')) %>% 
  mutate(Choice_lagWheat = as.numeric(Choice == 'Wheat')) %>% 
  mutate(basin_BearRiver = as.numeric(basin == 'BEAR RIVER')) %>% 
  mutate(basin_WeberRiver = as.numeric(basin == 'WEBER RIVER')) %>% 
  mutate(basin_JordanRiver = as.numeric(basin == 'JORDAN RIVER')) %>% 
  mutate(basin_NonGSL = as.numeric(basin == 'NON GSL')) %>% 
  dplyr::select(pid,basin,county,irri,Choice,
                rent_Fallow,rent_Alfalfa,rent_Hay,rent_OtherCrop,rent_OtherGrain,rent_Wheat,
                Choice_lagFallow,Choice_lagHay,Choice_lagOtherCrop,Choice_lagOtherGrain,Choice_lagWheat,
                basin_BearRiver,basin_WeberRiver,basin_JordanRiver,basin_NonGSL,X,Y)

### Calculate the WTAs, costs, and bootstraping CI ### 
# wta and cost
x_rent <- as.matrix(alfalfa_2022[,c('rent_Alfalfa','rent_Hay','rent_OtherCrop','rent_OtherGrain','rent_Wheat')])
x_othr <- as.matrix(alfalfa_2022[,c('Choice_lagFallow','Choice_lagHay','Choice_lagOtherCrop','Choice_lagOtherGrain','Choice_lagWheat','basin_BearRiver','basin_WeberRiver','basin_JordanRiver','basin_NonGSL')])
beta.v <- as.matrix(mod$coefficients)
# Function to extract coefficients from the model
coef <- function(beta,k){
  alpha <- rep(beta[1],k)
  gamma <- beta[2:nrow(beta)]
  m <- length(gamma)/k
  matrix <- matrix(0, nrow = m, ncol = k)
  # Fill the matrix row by row
  for (i in 1:m) {
    for (j in 1:k) {
      # Calculate index for vector
      index <- (i - 1) * k + j
      # Fill the matrix
      if (index <= length(gamma)) {
        matrix[i, j] <- gamma[index]
      }
    }
  }
  beta.m <- rbind(alpha,matrix)
  rownames(beta.m)<-NULL 
  return(beta.m)
}
# Function to calculate wta, water-saving potential, and unit water-saving cost
wta <- function(beta.v,x_rent,drip,sprinkler,flood){
  coef_values <- coef(beta.v,5)
  lin_preds <- x_rent * coef_values[1,] + x_othr %*% coef_values[-1,]
  lin_preds <- cbind(rep(0,nrow(lin_preds)),lin_preds)
  colnames(lin_preds) <- c('Fallow','Alfalfa','Hay','OtherCrop','Grain','Wheat')
  wta <- cbind(alfalfa_2022,lin_preds) %>% 
    filter(basin_BearRiver+basin_WeberRiver+basin_JordanRiver==1) %>% 
    mutate(AppEff = case_when(irri == 'Drip' ~ drip, irri == 'Sprinkler' ~ sprinkler, irri == 'Flood' ~ flood)) %>% 
    drop_na(AppEff) %>% 
    mutate(wta_Fallow = 1000*(Alfalfa-Fallow)*2.47105/mod$coefficients[1]) %>% # convert $1,000/acre to $/acre, and further convert to $/ha (1 ha = 2.47105 acres)
    mutate(wta_Grain = 1000*(Alfalfa-Grain)*2.47105/mod$coefficients[1]) %>% 
    mutate(wta_Hay = 1000*(Alfalfa-Hay)*2.47105/mod$coefficients[1]) %>% 
    mutate(water_Fallow = 28.73*22.86/AppEff) %>% # 1 inch = 0.0254 meter, 1 m2 inch = 0.0254 m3, 22.86 m3/pixel (22.86=.0254*900), effective rainfall is considered (80%)
    mutate(water_Grain = (28.73-19.21)*22.86/AppEff) %>% 
    mutate(water_Hay = (28.73-20.81)*22.86/AppEff) %>% 
    mutate(cost_Fallow = wta_Fallow*.09/water_Fallow) %>% # $/m3, 1 ha = 10000 m2, 1 pixel = 900 m2, .09 = 900 m2/10000 m2
    mutate(cost_Grain = wta_Grain*.09/water_Grain) %>% 
    mutate(cost_Hay = wta_Hay*.09/water_Hay) %>%  
    dplyr::select(pid,basin,county,irri,AppEff,Choice,wta_Fallow,wta_Grain,wta_Hay,
                  cost_Fallow,cost_Grain,cost_Hay,water_Fallow,water_Grain,water_Hay,X,Y) 
  return(wta)
}

alfalfa_wta <- wta(beta.v,x_rent,drip=.90,sprinkler=.80,flood=.76) # Used to produce Figs 1 and 2

# Check if conversion to fallow is the most cost-effective option
indicator <- ifelse(alfalfa_wta[,'cost_Fallow'] != pmin(alfalfa_wta[,'cost_Fallow'], alfalfa_wta[,'cost_Grain'], alfalfa_wta[,'cost_Hay']), 1, 0)
sum(indicator) # For each pixel, the cost of converting to fallow is the lowest

# Create upper and lower bounds of 95% CI for WTA and cost using bootstrapped sample
beta_bootstrap <- import('~/beta_bootstrap.csv')
beta_bootstrap <- as.matrix(beta_bootstrap)

cost_fallow <- matrix(0, nrow = nrow(alfalfa_wta), ncol = 1000)
cost_grain <- matrix(0, nrow = nrow(alfalfa_wta), ncol = 1000)
cost_hay <- matrix(0, nrow = nrow(alfalfa_wta), ncol = 1000)
for (j in 1:1000){
  tmp <- wta(matrix(beta_bootstrap[j,]),x_rent,drip=.90,sprinkler=.80,flood=.76)
  cost_fallow[,j] <- tmp$cost_Fallow
  cost_grain[,j] <- tmp$cost_Grain
  cost_hay[,j] <- tmp$cost_Hay
}
cost_fallow_bca <- matrix(0, nrow = nrow(alfalfa_wta), ncol = 2)
cost_grain_bca <- matrix(0, nrow = nrow(alfalfa_wta), ncol = 2)
cost_hay_bca <- matrix(0, nrow = nrow(alfalfa_wta), ncol = 2)
for (i in 1:nrow(alfalfa_wta)){
  cost_fallow_bca[i,] <- coxed::bca(cost_fallow[i,], conf.level = .95)
  cost_grain_bca[i,] <- coxed::bca(cost_grain[i,], conf.level = .95)
  cost_hay_bca[i,] <- coxed::bca(cost_hay[i,], conf.level = .95)
}
colnames(cost_fallow_bca) <- c('cost_Fallow_lower','cost_Fallow_upper')
colnames(cost_grain_bca) <- c('cost_Grain_lower','cost_Grain_upper')
colnames(cost_hay_bca) <- c('cost_Hay_lower','cost_Hay_upper')
alfalfa_wta <- cbind(alfalfa_wta,cost_fallow_bca,cost_grain_bca,cost_hay_bca) %>% 
  mutate(wta_Fallow_lower = cost_Fallow_lower*water_Fallow/.09,wta_Fallow_upper = cost_Fallow_upper*water_Fallow/.09,
         wta_Grain_lower = cost_Grain_lower*water_Grain/.09,wta_Grain_upper = cost_Grain_upper*water_Grain/.09,
         wta_Hay_lower = cost_Hay_lower*water_Hay/.09,wta_Hay_upper = cost_Hay_upper*water_Hay/.09) %>% 
  dplyr::select(pid,basin,county,irri,AppEff,Choice,X,Y,wta_Fallow,wta_Grain,wta_Hay,cost_Fallow,cost_Grain,cost_Hay,
                water_Fallow,water_Grain,water_Hay,wta_Fallow_lower,wta_Fallow_upper,wta_Grain_lower,wta_Grain_upper,wta_Hay_lower,wta_Hay_upper,
                cost_Fallow_lower,cost_Fallow_upper,cost_Grain_lower,cost_Grain_upper,cost_Hay_lower,cost_Hay_upper)

# Data frame 'alfalfa_wta' has 28 columns and 752,212 rows; export the data in csv format
write.csv(alfalfa_wta,'~/alfalfa_wta.csv',row.names=FALSE)

### Calculate marginal costs and draw marginal cost curves (i.e., water supply curves in Fig. 3) ### 
# Marginal costs
mc <- function(wta){
  mc_list <- list()
  
  mc_Fallow <- wta %>% 
    arrange(cost_Fallow,water_Fallow) %>% 
    mutate(group = 'Fallow') %>% 
    mutate(water=water_Fallow/1000000) %>% # million m3
    rename(cost=cost_Fallow,cost_lower=cost_Fallow_lower,cost_upper=cost_Fallow_upper) %>% #$/m3
    dplyr::select(basin,county,group,water,cost,cost_lower,cost_upper) 
  
  mc_Grain <-  wta %>% 
    arrange(cost_Grain,water_Grain) %>% 
    mutate(group = 'Spring Grains') %>% 
    mutate(water=water_Grain/1000000) %>% 
    rename(cost=cost_Grain,cost_lower=cost_Grain_lower,cost_upper=cost_Grain_upper) %>% 
    dplyr::select(basin,county,group,water,cost,cost_lower,cost_upper) 
  
  mc_Hay <- wta %>% 
    arrange(cost_Hay,water_Hay) %>% 
    mutate(group = 'Hay') %>% 
    mutate(water=water_Hay/1000000) %>% 
    rename(cost=cost_Hay,cost_lower=cost_Hay_lower,cost_upper=cost_Hay_upper) %>% 
    dplyr::select(basin,county,group,water,cost,cost_lower,cost_upper) 
  
  mc_list[[1]] <- rbind(subset(mc_Fallow,basin=='BEAR RIVER'),
                        subset(mc_Grain,basin=='BEAR RIVER'),subset(mc_Hay,basin=='BEAR RIVER')) %>% 
    group_by(group) %>% 
    mutate(water_cum=cumsum(water))
  
  mc_list[[2]] <- rbind(subset(mc_Fallow,basin=='WEBER RIVER'),
                        subset(mc_Grain,basin=='WEBER RIVER'),subset(mc_Hay,basin=='WEBER RIVER')) %>% 
    group_by(group) %>% 
    mutate(water_cum=cumsum(water))
  
  mc_list[[3]] <- rbind(subset(mc_Fallow,basin=='JORDAN RIVER'),
                        subset(mc_Grain,basin=='JORDAN RIVER'),subset(mc_Hay,basin=='JORDAN RIVER')) %>% 
    group_by(group) %>% 
    mutate(water_cum=cumsum(water))
  
  return(mc_list) 
}
mc_curve <- function(mc){
  mc_Bear <- mc[[1]]
  mc_Weber <- mc[[2]]
  mc_Jordan <- mc[[3]]
  y_limits <- c(min(mc_Bear$cost,mc_Weber$cost,mc_Jordan$cost), max(mc_Bear$cost,mc_Weber$cost,mc_Jordan$cost))
  p1 <- ggplot(mc_Bear, aes(x = water_cum, y = cost, color = group)) +
    geom_ribbon(aes(ymin = cost_lower, ymax = cost_upper), fill = "grey70", alpha = 0.5, lty = 'blank') + # Add confidence interval
    geom_line(lwd = .5) +
    labs(x = 'Conserved water (million cubic meter)', y = "Marginal cost ($/cubic meter)", 
         title = 'Bear') +
    scale_y_continuous(limits = y_limits, breaks = seq(0, y_limits[2], by = .1)) + 
    theme_minimal()
  p2 <- ggplot(mc_Weber, aes(x = water_cum, y = cost, color = group)) +
    geom_ribbon(aes(ymin = cost_lower, ymax = cost_upper), fill = "grey70", alpha = 0.5, lty = 'blank') + # Add confidence interval
    geom_line(lwd = .5) +
    labs(x = 'Conserved water (million cubic meter)', y = "", 
         title = "Weber") +
    scale_y_continuous(limits = y_limits, breaks = seq(0, y_limits[2], by = .1)) + 
    theme_minimal()
  p3 <- ggplot(mc_Jordan, aes(x = water_cum, y = cost, color = group)) +
    geom_ribbon(aes(ymin = cost_lower, ymax = cost_upper), fill = "grey70", alpha = 0.5, lty = 'blank') + # Add confidence interval
    geom_line(lwd = .5) +
    labs(x = 'Conserved water (million cubic meter)', y = "", 
         title = "Jordan") +
    scale_y_continuous(limits = y_limits, breaks = seq(0, y_limits[2], by = .1)) + 
    theme_minimal()
  my_list <- list(p1,p2,p3)
  return(my_list)
}
alfalfa_mc <- mc(alfalfa_wta)
outplots <- mc_curve(alfalfa_mc)
ggarrange(outplots[[1]],outplots[[2]],outplots[[3]],nrow=1,ncol=3,common.legend = TRUE, legend = "bottom")

########### Sensitivity Test ###########
### Sensitivity to changes in crop net returns ###
drops <- c('cost_Fallow_lower','cost_Grain_lower','cost_Hay_lower','cost_Fallow_upper','cost_Grain_upper','cost_Hay_upper')
# Alfalfa net returns
x_rent_low <- x_rent
x_rent_low[,1] <- x_rent[,1]*.9
cost_alfa_low <- wta(beta.v,x_rent_low,drip=.90,sprinkler=.80,flood=.76) %>% 
  rename(cost_Fallow_lower=cost_Fallow,cost_Grain_lower=cost_Grain,cost_Hay_lower=cost_Hay) %>% 
  dplyr::select(pid,cost_Fallow_lower,cost_Grain_lower,cost_Hay_lower)
x_rent_high <- x_rent
x_rent_high[,1] <- x_rent[,1]*1.1
cost_alfa_high <- wta(beta.v,x_rent_high,drip=.90,sprinkler=.80,flood=.76) %>% 
  rename(cost_Fallow_upper=cost_Fallow,cost_Grain_upper=cost_Grain,cost_Hay_upper=cost_Hay) %>% 
  dplyr::select(pid,cost_Fallow_upper,cost_Grain_upper,cost_Hay_upper)
alfalfa_wta_drop <- alfalfa_wta[, !(names(alfalfa_wta) %in% drops)]
alfalfa_wta_alfa <- Reduce(function(x, y) merge(x, y, by = 'pid', all = TRUE), list(alfalfa_wta_drop,cost_alfa_low,cost_alfa_high))
alfalfa_mc_alfa <- mc(alfalfa_wta_alfa)

# Spring grain net returns 
x_rent_low <- x_rent
x_rent_low[,4] <- x_rent[,4]*1.1
cost_grain_low <- wta(beta.v,x_rent_low,drip=.90,sprinkler=.80,flood=.76) %>% 
  rename(cost_Fallow_lower=cost_Fallow,cost_Grain_lower=cost_Grain,cost_Hay_lower=cost_Hay) %>% 
  dplyr::select(pid,cost_Fallow_lower,cost_Grain_lower,cost_Hay_lower)
x_rent_high <- x_rent
x_rent_high[,4] <- x_rent[,4]*.9
cost_grain_high <- wta(beta.v,x_rent_high,drip=.90,sprinkler=.80,flood=.76) %>% 
  rename(cost_Fallow_upper=cost_Fallow,cost_Grain_upper=cost_Grain,cost_Hay_upper=cost_Hay) %>% 
  dplyr::select(pid,cost_Fallow_upper,cost_Grain_upper,cost_Hay_upper)
alfalfa_wta_grain <- Reduce(function(x, y) merge(x, y, by = 'pid', all = TRUE), list(alfalfa_wta_drop,cost_grain_low,cost_grain_high))
alfalfa_mc_grain <- mc(alfalfa_wta_grain)

# Hay net returns
x_rent_low <- x_rent
x_rent_low[,2] <- x_rent[,2]*1.1
cost_hay_low <- wta(beta.v,x_rent_low,drip=.90,sprinkler=.80,flood=.76) %>% 
  rename(cost_Fallow_lower=cost_Fallow,cost_Grain_lower=cost_Grain,cost_Hay_lower=cost_Hay) %>% 
  dplyr::select(pid,cost_Fallow_lower,cost_Grain_lower,cost_Hay_lower)
x_rent_high <- x_rent
x_rent_high[,2] <- x_rent[,2]*.9
cost_hay_high <- wta(beta.v,x_rent_high,drip=.90,sprinkler=.80,flood=.76) %>% 
  rename(cost_Fallow_upper=cost_Fallow,cost_Grain_upper=cost_Grain,cost_Hay_upper=cost_Hay) %>% 
  dplyr::select(pid,cost_Fallow_upper,cost_Grain_upper,cost_Hay_upper)
alfalfa_wta_hay <- Reduce(function(x, y) merge(x, y, by = 'pid', all = TRUE), list(alfalfa_wta_drop,cost_hay_low,cost_hay_high))
alfalfa_mc_hay <- mc(alfalfa_wta_hay)

### Sensitivity to assumptions of net irrigation requirement (NIR) ###
wta_nir <- function(beta.v,x_rent){
  coef_values <- coef(beta.v,5)
  lin_preds <- x_rent * coef_values[1,] + x_othr %*% coef_values[-1,]
  lin_preds <- cbind(rep(0,nrow(lin_preds)),lin_preds)
  colnames(lin_preds) <- c('Fallow','Alfalfa','Hay','OtherCrop','Grain','Wheat')
  wta <- cbind(alfalfa_rent,lin_preds) %>% 
    filter(basin_BearRiver+basin_WeberRiver+basin_JordanRiver==1) %>% 
    mutate(AppEff = case_when(irri == 'Drip' ~ .90, irri == 'Sprinkler' ~ .80, irri == 'Flood' ~ .76)) %>% 
    drop_na(AppEff) %>% 
    mutate(wta_Fallow = 1000*(Alfalfa-Fallow)*2.47105/mod_c2$coefficients[1]) %>% # convert $1,000/acre to $/acre, and further convert to $/ha (1 ha = 2.47105 acres)
    mutate(wta_Grain = 1000*(Alfalfa-Grain)*2.47105/mod_c2$coefficients[1]) %>% 
    mutate(wta_Hay = 1000*(Alfalfa-Hay)*2.47105/mod_c2$coefficients[1]) %>% 
    mutate(water_Fallow = 28.73*22.86/AppEff) %>% # 1 inch = 0.0254 meter, 1 m2 inch = 0.0254 m3, 22.86 m3/pixel (22.86=.0254*900), effective rainfall is considered (80%)
    mutate(water_Grain = (28.73-19.21)*22.86/AppEff) %>% 
    mutate(water_Hay = (28.73-20.81)*22.86/AppEff) %>% 
    mutate(cost_Fallow = wta_Fallow*.09/water_Fallow) %>% # $/m3, 1 ha = 10000 m2, 1 pixel = 900 m2, .09 = 900 m2/10000 m2
    mutate(cost_Grain = wta_Grain*.09/water_Grain) %>% 
    mutate(cost_Hay = wta_Hay*.09/water_Hay) %>% 
    mutate(water_Fallow_low = 28.73*22.86/AppEff*.9) %>% # 1 inch = 0.0254 meter, 1 m2 inch = 0.0254 m3, 22.86 m3/pixel (22.86=.0254*900), effective rainfall is considered (80%)
    mutate(water_Grain_low = (28.73-19.21)*22.86/AppEff*.9) %>% 
    mutate(water_Hay_low = (28.73-20.81)*22.86/AppEff*.9) %>% 
    mutate(cost_Fallow_low = wta_Fallow*.09/water_Fallow_low) %>% # $/m3, 1 ha = 10000 m2, 1 pixel = 900 m2, .09 = 900 m2/10000 m2
    mutate(cost_Grain_low = wta_Grain*.09/water_Grain_low) %>% 
    mutate(cost_Hay_low = wta_Hay*.09/water_Hay_low) %>%     
    mutate(water_Fallow_high = 28.73*22.86/AppEff*1.1) %>% # 1 inch = 0.0254 meter, 1 m2 inch = 0.0254 m3, 22.86 m3/pixel (22.86=.0254*900), effective rainfall is considered (80%)
    mutate(water_Grain_high = (28.73-19.21)*22.86/AppEff*1.1) %>% 
    mutate(water_Hay_high = (28.73-20.81)*22.86/AppEff*1.1) %>% 
    mutate(cost_Fallow_high = wta_Fallow*.09/water_Fallow_high) %>% # $/m3, 1 ha = 10000 m2, 1 pixel = 900 m2, .09 = 900 m2/10000 m2
    mutate(cost_Grain_high = wta_Grain*.09/water_Grain_high) %>% 
    mutate(cost_Hay_high = wta_Hay*.09/water_Hay_high) %>%  
    dplyr::select(pid,basin,county,irri,AppEff,Choice,X,Y,wta_Fallow,wta_Grain,wta_Hay,
                  cost_Fallow,cost_Grain,cost_Hay,water_Fallow,water_Grain,water_Hay,
                  cost_Fallow_low,cost_Grain_low,cost_Hay_low,water_Fallow_low,water_Grain_low,water_Hay_low,
                  cost_Fallow_high,cost_Grain_high,cost_Hay_high,water_Fallow_high,water_Grain_high,water_Hay_high) 
  return(wta)
}
mc_nir <- function(wta){
  mc_list <- list()
  
  mc_Fallow <- wta %>% 
    arrange(cost_Fallow,water_Fallow) %>% 
    mutate(group = 'Fallow') %>% 
    mutate(scenario = 'Baseline') %>% 
    mutate(water=water_Fallow/1000000) %>% # million m3
    rename(cost=cost_Fallow) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Fallow_low <- wta %>% 
    arrange(cost_Fallow_low,water_Fallow_low) %>% 
    mutate(group = 'Fallow') %>% 
    mutate(scenario = 'NIR decreases by 10%') %>% 
    mutate(water=water_Fallow_low/1000000) %>% # million m3
    rename(cost=cost_Fallow_low) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Fallow_high <- wta %>% 
    arrange(cost_Fallow_high,water_Fallow_high) %>% 
    mutate(group = 'Fallow') %>% 
    mutate(scenario = 'NIR increases by 10%') %>% 
    mutate(water=water_Fallow_high/1000000) %>% # million m3
    rename(cost=cost_Fallow_high) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Grain <- wta %>% 
    arrange(cost_Grain,water_Grain) %>% 
    mutate(group = 'Spring Grains') %>% 
    mutate(scenario = 'Baseline') %>% 
    mutate(water=water_Grain/1000000) %>% # million m3
    rename(cost=cost_Grain) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Grain_low <- wta %>% 
    arrange(cost_Grain_low,water_Grain_low) %>% 
    mutate(group = 'Spring Grains') %>% 
    mutate(scenario = 'NIR decreases by 10%') %>% 
    mutate(water=water_Grain_low/1000000) %>% # million m3
    rename(cost=cost_Grain_low) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Grain_high <- wta %>% 
    arrange(cost_Grain_high,water_Grain_high) %>% 
    mutate(group = 'Spring Grains') %>% 
    mutate(scenario = 'NIR increases by 10%') %>% 
    mutate(water=water_Grain_high/1000000) %>% # million m3
    rename(cost=cost_Grain_high) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Hay <- wta %>% 
    arrange(cost_Hay,water_Hay) %>% 
    mutate(group = 'Hay') %>% 
    mutate(scenario = 'Baseline') %>% 
    mutate(water=water_Hay/1000000) %>% # million m3
    rename(cost=cost_Hay) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Hay_low <- wta %>% 
    arrange(cost_Hay_low,water_Hay_low) %>% 
    mutate(group = 'Hay') %>% 
    mutate(scenario = 'NIR decreases by 10%') %>% 
    mutate(water=water_Hay_low/1000000) %>% # million m3
    rename(cost=cost_Hay_low) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Hay_high <- wta %>% 
    arrange(cost_Hay_high,water_Hay_high) %>% 
    mutate(group = 'Hay') %>%  
    mutate(scenario = 'NIR increases by 10%') %>% 
    mutate(water=water_Hay_high/1000000) %>% # million m3
    rename(cost=cost_Hay_high) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_list[[1]] <- rbind(subset(mc_Fallow,basin=='BEAR RIVER'),subset(mc_Fallow_low,basin=='BEAR RIVER'),subset(mc_Fallow_high,basin=='BEAR RIVER'),
                        subset(mc_Grain,basin=='BEAR RIVER'),subset(mc_Grain_low,basin=='BEAR RIVER'),subset(mc_Grain_high,basin=='BEAR RIVER'),
                        subset(mc_Hay,basin=='BEAR RIVER'),subset(mc_Hay_low,basin=='BEAR RIVER'),subset(mc_Hay_high,basin=='BEAR RIVER')) %>% 
    group_by(group,scenario) %>% 
    mutate(water_cum=cumsum(water))
  
  mc_list[[2]] <- rbind(subset(mc_Fallow,basin=='WEBER RIVER'),subset(mc_Fallow_low,basin=='WEBER RIVER'),subset(mc_Fallow_high,basin=='WEBER RIVER'),
                        subset(mc_Grain,basin=='WEBER RIVER'),subset(mc_Grain_low,basin=='WEBER RIVER'),subset(mc_Grain_high,basin=='WEBER RIVER'),
                        subset(mc_Hay,basin=='WEBER RIVER'),subset(mc_Hay_low,basin=='WEBER RIVER'),subset(mc_Hay_high,basin=='WEBER RIVER')) %>% 
    group_by(group,scenario) %>% 
    mutate(water_cum=cumsum(water))
  
  mc_list[[3]] <- rbind(subset(mc_Fallow,basin=='JORDAN RIVER'),subset(mc_Fallow_low,basin=='JORDAN RIVER'),subset(mc_Fallow_high,basin=='JORDAN RIVER'),
                        subset(mc_Grain,basin=='JORDAN RIVER'),subset(mc_Grain_low,basin=='JORDAN RIVER'),subset(mc_Grain_high,basin=='JORDAN RIVER'),
                        subset(mc_Hay,basin=='JORDAN RIVER'),subset(mc_Hay_low,basin=='JORDAN RIVER'),subset(mc_Hay_high,basin=='JORDAN RIVER')) %>% 
    group_by(group,scenario) %>% 
    mutate(water_cum=cumsum(water))
  
  return(mc_list) 
}
alfalfa_wta_nir <- wta_nir(beta.v,x_rent)
alfalfa_mc_nir <- mc_nir(alfalfa_wta_nir)

### Sensitivity to assumptions of irrigation efficiency (Ea) ###
wta_irri <- function(beta.v,x_rent){
  coef_values <- coef(beta.v,5)
  lin_preds <- x_rent * coef_values[1,] + x_othr %*% coef_values[-1,]
  lin_preds <- cbind(rep(0,nrow(lin_preds)),lin_preds)
  colnames(lin_preds) <- c('Fallow','Alfalfa','Hay','OtherCrop','Grain','Wheat')
  wta <- cbind(alfalfa_rent,lin_preds) %>% 
    filter(basin_BearRiver+basin_WeberRiver+basin_JordanRiver==1) %>% 
    mutate(AppEff = case_when(irri == 'Drip' ~ .90, irri == 'Sprinkler' ~ .80, irri == 'Flood' ~ .76)) %>% 
    drop_na(AppEff) %>% 
    mutate(AppEff_l = case_when(irri == 'Drip' ~ .85, irri == 'Sprinkler' ~ .75, irri == 'Flood' ~ .71)) %>% 
    mutate(AppEff_h = case_when(irri == 'Drip' ~ .95, irri == 'Sprinkler' ~ .85, irri == 'Flood' ~ .81)) %>% 
    mutate(wta_Fallow = 1000*(Alfalfa-Fallow)*2.47105/mod_c2$coefficients[1]) %>% # convert $1,000/acre to $/acre, and further convert to $/ha (1 ha = 2.47105 acres)
    mutate(wta_Grain = 1000*(Alfalfa-Grain)*2.47105/mod_c2$coefficients[1]) %>% 
    mutate(wta_Hay = 1000*(Alfalfa-Hay)*2.47105/mod_c2$coefficients[1]) %>% 
    mutate(water_Fallow = 28.73*22.86/AppEff) %>% # 1 inch = 0.0254 meter, 1 m2 inch = 0.0254 m3, 22.86 m3/pixel (22.86=.0254*900), effective rainfall is considered (80%)
    mutate(water_Grain = (28.73-19.21)*22.86/AppEff) %>% 
    mutate(water_Hay = (28.73-20.81)*22.86/AppEff) %>% 
    mutate(cost_Fallow = wta_Fallow*.09/water_Fallow) %>%   # $/m3, 1 ha = 10000 m2, 1 pixel = 900 m2, .09 = 900 m2/10000 m2
    mutate(cost_Grain = wta_Grain*.09/water_Grain) %>% 
    mutate(cost_Hay = wta_Hay*.09/water_Hay) %>% 
    mutate(water_Fallow_low = 28.73*22.86/AppEff_l) %>% # 1 inch = 0.0254 meter, 1 m2 inch = 0.0254 m3, 22.86 m3/pixel (22.86=.0254*900), effective rainfall is considered (80%)
    mutate(water_Grain_low = (28.73-19.21)*22.86/AppEff_l) %>% 
    mutate(water_Hay_low = (28.73-20.81)*22.86/AppEff_l) %>% 
    mutate(cost_Fallow_low = wta_Fallow*.09/water_Fallow_low) %>%   # $/m3, 1 ha = 10000 m2, 1 pixel = 900 m2, .09 = 900 m2/10000 m2
    mutate(cost_Grain_low = wta_Grain*.09/water_Grain_low) %>% 
    mutate(cost_Hay_low = wta_Hay*.09/water_Hay_low) %>%     
    mutate(water_Fallow_high = 28.73*22.86/AppEff_h) %>% # 1 inch = 0.0254 meter, 1 m2 inch = 0.0254 m3, 22.86 m3/pixel (22.86=.0254*900), effective rainfall is considered (80%)
    mutate(water_Grain_high = (28.73-19.21)*22.86/AppEff_h) %>% 
    mutate(water_Hay_high = (28.73-20.81)*22.86/AppEff_h) %>% 
    mutate(cost_Fallow_high = wta_Fallow*.09/water_Fallow_high) %>%   # $/m3, 1 ha = 10000 m2, 1 pixel = 900 m2, .09 = 900 m2/10000 m2
    mutate(cost_Grain_high = wta_Grain*.09/water_Grain_high) %>% 
    mutate(cost_Hay_high = wta_Hay*.09/water_Hay_high) %>%  
    dplyr::select(pid,basin,county,irri,AppEff,Choice,X,Y,wta_Fallow,wta_Grain,wta_Hay,
                  cost_Fallow,cost_Grain,cost_Hay,water_Fallow,water_Grain,water_Hay,
                  cost_Fallow_low,cost_Grain_low,cost_Hay_low,water_Fallow_low,water_Grain_low,water_Hay_low,
                  cost_Fallow_high,cost_Grain_high,cost_Hay_high,water_Fallow_high,water_Grain_high,water_Hay_high) 
  return(wta)
}
alfalfa_wta_irri <- wta_irri(beta.v,x_rent)
mc_irri <- function(wta){
  mc_list <- list()
  
  mc_Fallow <- wta %>% 
    arrange(cost_Fallow,water_Fallow) %>% 
    mutate(group = 'Fallow') %>% 
    mutate(scenario = 'Baseline') %>% 
    mutate(water=water_Fallow/1000000) %>% # million m3
    rename(cost=cost_Fallow) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Fallow_low <- wta %>% 
    arrange(cost_Fallow_low,water_Fallow_low) %>% 
    mutate(group = 'Fallow') %>% 
    mutate(scenario = 'Ea decreases by 5 pp') %>% 
    mutate(water=water_Fallow_low/1000000) %>% # million m3
    rename(cost=cost_Fallow_low) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Fallow_high <- wta %>% 
    arrange(cost_Fallow_high,water_Fallow_high) %>% 
    mutate(group = 'Fallow') %>% 
    mutate(scenario = 'Ea increases by 5 pp') %>% 
    mutate(water=water_Fallow_high/1000000) %>% # million m3
    rename(cost=cost_Fallow_high) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Grain <- wta %>% 
    arrange(cost_Grain,water_Grain) %>% 
    mutate(group = 'Spring Grains') %>% 
    mutate(scenario = 'Baseline') %>% 
    mutate(water=water_Grain/1000000) %>% # million m3
    rename(cost=cost_Grain) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Grain_low <- wta %>% 
    arrange(cost_Grain_low,water_Grain_low) %>% 
    mutate(group = 'Spring Grains') %>% 
    mutate(scenario = 'Ea decreases by 5 pp') %>% 
    mutate(water=water_Grain_low/1000000) %>% # million m3
    rename(cost=cost_Grain_low) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Grain_high <- wta %>% 
    arrange(cost_Grain_high,water_Grain_high) %>% 
    mutate(group = 'Spring Grains') %>% 
    mutate(scenario = 'Ea increases by 5 pp') %>% 
    mutate(water=water_Grain_high/1000000) %>% # million m3
    rename(cost=cost_Grain_high) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Hay <- wta %>% 
    arrange(cost_Hay,water_Hay) %>% 
    mutate(group = 'Hay') %>% 
    mutate(scenario = 'Baseline') %>% 
    mutate(water=water_Hay/1000000) %>% # million m3
    rename(cost=cost_Hay) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Hay_low <- wta %>% 
    arrange(cost_Hay_low,water_Hay_low) %>% 
    mutate(group = 'Hay') %>% 
    mutate(scenario = 'Ea decreases by 5 pp') %>% 
    mutate(water=water_Hay_low/1000000) %>% # million m3
    rename(cost=cost_Hay_low) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_Hay_high <- wta %>% 
    arrange(cost_Hay_high,water_Hay_high) %>% 
    mutate(group = 'Hay') %>%  
    mutate(scenario = 'Ea increases by 5 pp') %>% 
    mutate(water=water_Hay_high/1000000) %>% # million m3
    rename(cost=cost_Hay_high) %>% 
    dplyr::select(basin,county,group,scenario,water,cost) 
  
  mc_list[[1]] <- rbind(subset(mc_Fallow,basin=='BEAR RIVER'),subset(mc_Fallow_low,basin=='BEAR RIVER'),subset(mc_Fallow_high,basin=='BEAR RIVER'),
                        subset(mc_Grain,basin=='BEAR RIVER'),subset(mc_Grain_low,basin=='BEAR RIVER'),subset(mc_Grain_high,basin=='BEAR RIVER'),
                        subset(mc_Hay,basin=='BEAR RIVER'),subset(mc_Hay_low,basin=='BEAR RIVER'),subset(mc_Hay_high,basin=='BEAR RIVER')) %>% 
    group_by(group,scenario) %>% 
    mutate(water_cum=cumsum(water))
  
  mc_list[[2]] <- rbind(subset(mc_Fallow,basin=='WEBER RIVER'),subset(mc_Fallow_low,basin=='WEBER RIVER'),subset(mc_Fallow_high,basin=='WEBER RIVER'),
                        subset(mc_Grain,basin=='WEBER RIVER'),subset(mc_Grain_low,basin=='WEBER RIVER'),subset(mc_Grain_high,basin=='WEBER RIVER'),
                        subset(mc_Hay,basin=='WEBER RIVER'),subset(mc_Hay_low,basin=='WEBER RIVER'),subset(mc_Hay_high,basin=='WEBER RIVER')) %>% 
    group_by(group,scenario) %>% 
    mutate(water_cum=cumsum(water))
  
  mc_list[[3]] <- rbind(subset(mc_Fallow,basin=='JORDAN RIVER'),subset(mc_Fallow_low,basin=='JORDAN RIVER'),subset(mc_Fallow_high,basin=='JORDAN RIVER'),
                        subset(mc_Grain,basin=='JORDAN RIVER'),subset(mc_Grain_low,basin=='JORDAN RIVER'),subset(mc_Grain_high,basin=='JORDAN RIVER'),
                        subset(mc_Hay,basin=='JORDAN RIVER'),subset(mc_Hay_low,basin=='JORDAN RIVER'),subset(mc_Hay_high,basin=='JORDAN RIVER')) %>% 
    group_by(group,scenario) %>% 
    mutate(water_cum=cumsum(water))
  
  return(mc_list) 
}
alfalfa_mc_irri <- mc_irri(alfalfa_wta_irri)

### Summary table (Table S5) ###
# Base model
mc_Bear <- alfalfa_mc[[1]]
mc_Weber <- alfalfa_mc[[2]]
mc_Jordan <- alfalfa_mc[[3]]
mc_combined <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc = cost*water, tc_lower = cost_lower*water, tc_upper = cost_upper*water)  # million dollar (= $/m3*million m3) 
tab0 <- mc_combined %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc,tc_lower,tc_upper), sum, na.rm=T) %>% 
  mutate(ac=tc/water,ac_lower=tc_lower/water,ac_upper=tc_upper/water) # $/m3

# Panel A: Crop net economic returns
mc_Bear <- alfalfa_mc_alfa[[1]]
mc_Weber <- alfalfa_mc_alfa[[2]]
mc_Jordan <- alfalfa_mc_alfa[[3]]
mc_combined_alfa <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc_alfa_low = cost_lower*water, tc_alfa_high = cost_upper*water)  # million dollar (= $/m3*million m3) 
tab61<- mc_combined_alfa %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc_alfa_low,tc_alfa_high), sum, na.rm=T) %>% 
  mutate(ac_alfa_low=tc_alfa_low/water,ac_alfa_high=tc_alfa_high/water) %>% dplyr::select(-water)

mc_Bear <- alfalfa_mc_grain[[1]]
mc_Weber <- alfalfa_mc_grain[[2]]
mc_Jordan <- alfalfa_mc_grain[[3]]
mc_combined_grain <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc_grain_low = cost_lower*water, tc_grain_high = cost_upper*water)  # million dollar (= $/m3*million m3) 
tab2 <- mc_combined_grain %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc_grain_low,tc_grain_high), sum, na.rm=T) %>% 
  mutate(ac_grain_low=tc_grain_low/water,ac_grain_high=tc_grain_high/water) %>% dplyr::select(-water)

mc_Bear <- alfalfa_mc_hay[[1]]
mc_Weber <- alfalfa_mc_hay[[2]]
mc_Jordan <- alfalfa_mc_hay[[3]]
mc_combined_hay <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc_hay_low = cost_lower*water, tc_hay_high = cost_upper*water)  # million dollar (= $/m3*million m3) 
tab3 <- mc_combined_hay %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc_hay_low,tc_hay_high), sum, na.rm=T) %>% 
  mutate(ac_hay_low=tc_hay_low/water,ac_hay_high=tc_hay_high/water) %>% dplyr::select(-water)

# Panel B: Crop water use
mc_Bear <- alfalfa_mc_nir[[1]] %>% filter(scenario == 'NIR decreases by 10%')
mc_Weber <- alfalfa_mc_nir[[2]] %>% filter(scenario == 'NIR decreases by 10%')
mc_Jordan <- alfalfa_mc_nir[[3]] %>% filter(scenario == 'NIR decreases by 10%')
mc_combined_nir_low <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc_nir_low = cost*water) # million dollar (= $/m3*million m3) 
tab4 <- mc_combined_nir_low %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc_nir_low), sum, na.rm=T) %>% 
  mutate(ac_nir_low=tc_nir_low/water) %>% 
  rename(water_nir_low=water)

mc_Bear <- alfalfa_mc_nir[[1]] %>% filter(scenario == 'NIR increases by 10%')
mc_Weber <- alfalfa_mc_nir[[2]] %>% filter(scenario == 'NIR increases by 10%')
mc_Jordan <- alfalfa_mc_nir[[3]] %>% filter(scenario == 'NIR increases by 10%')
mc_combined_nir_high <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc_nir_high = cost*water) # million dollar (= $/m3*million m3) 
tab5 <- mc_combined_nir_high %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc_nir_high), sum, na.rm=T) %>% 
  mutate(ac_nir_high=tc_nir_high/water) %>% 
  rename(water_nir_high=water)

mc_Bear <- alfalfa_mc_irri[[1]] %>% filter(scenario == 'Ea decreases by 5 pp')
mc_Weber <- alfalfa_mc_irri[[2]] %>% filter(scenario == 'Ea decreases by 5 pp')
mc_Jordan <- alfalfa_mc_irri[[3]] %>% filter(scenario == 'Ea decreases by 5 pp')
mc_combined_irri_low <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc_irri_low = cost*water) # million dollar (= $/m3*million m3) 
tab6 <- mc_combined_irri_low %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc_irri_low), sum, na.rm=T) %>% 
  mutate(ac_irri_low=tc_irri_low/water) %>% 
  rename(water_irri_low=water)

mc_Bear <- alfalfa_mc_irri[[1]] %>% filter(scenario == 'Ea increases by 5 pp')
mc_Weber <- alfalfa_mc_irri[[2]] %>% filter(scenario == 'Ea increases by 5 pp')
mc_Jordan <- alfalfa_mc_irri[[3]] %>% filter(scenario == 'Ea increases by 5 pp')
mc_combined_irri_high <- rbind(mc_Bear,mc_Weber,mc_Jordan) %>% 
  mutate(tc_irri_high = cost*water) # million dollar (= $/m3*million m3) 
tab7 <- mc_combined_irri_high %>% 
  group_by(group,basin) %>% 
  summarise_at(vars(water,tc_irri_high), sum, na.rm=T) %>% 
  mutate(ac_irri_high=tc_irri_high/water) %>% 
  rename(water_irri_high=water)

tab_combined <- Reduce(function(x, y) merge(x, y, by = c("group", "basin"), all = TRUE), list(tab0, tab1, tab2, tab3, tab4, tab5, tab6, tab7))

########### Water conservation policy design ###########
# Histogram plot of estimated WTA for fallow alfalfa by watershed (Fig. S3)
basin1 <- as.matrix(subset(alfalfa_wta,basin=='BEAR RIVER')[,'wta_Fallow'] )
basin2 <- as.matrix(subset(alfalfa_wta,basin=='WEBER RIVER')[,'wta_Fallow'])
basin3 <- as.matrix(subset(alfalfa_wta,basin=='JORDAN RIVER')[,'wta_Fallow'])
hist(basin1,prob=T,border=F,col='darkblue', density=80, xlim=c(3200,6700), breaks = 100, xlab='WTA ($/ha)', main = '')
hist(basin2,prob=T,border=F,col='darkred', density=80, xlim=c(3200,6700), breaks = 100, add=T)
hist(basin3,prob=T,border=F,col='darkgreen', density=80, xlim=c(3200,6700), breaks = 100, add=T)
legend('topleft',c('BEAR','WEBER','JORDAN'), fill = c('darkblue','darkred','darkgreen'), bty = 'n', border = NA)

# Watershed-level payment system (used to produce Fig 4)
cutoff_wta_basin <- data.frame(
  county = c('BOX ELDER','CACHE','RICH','WEBER','DAVIS','MORGAN','SUMMIT','SALT LAKE','UTAH','WASATCH','JUAB'),
  cutoff_wta = c(rep(6093,7),rep(4331,4)),
  cutoff_wta_upper = c(rep(6711,7),rep(4567,4))
)

alfalfa_wta_basin <- merge(alfalfa_wta,cutoff_wta_basin,by = 'county',all.x=T) %>% 
  mutate(enroll = as.numeric(wta_Fallow<=cutoff_wta), # Baseline scenario
         enroll_lower = as.numeric(wta_Fallow_lower<=cutoff_wta), # Lower bound of 95% CI in the baseline scenario
         enroll_upper = as.numeric(wta_Fallow_upper<=cutoff_wta)) %>% # Upper bound of 95% CI in the baseline scenario
  mutate(enroll_a = as.numeric(wta_Fallow<=cutoff_wta_upper), # Conservative scenario
         enroll_a_lower = as.numeric(wta_Fallow_lower<=cutoff_wta_upper), # Lower bound of 95% CI in the conservative scenario
         enroll_a_upper = as.numeric(wta_Fallow_upper<=cutoff_wta_upper)) %>% # Upper bound of 95% CI in the conservative scenario
  dplyr::select(X,Y,basin,county,wta_Fallow,enroll,enroll_lower,enroll_upper,enroll_a,enroll_a_lower,enroll_a_upper)
write.csv(alfalfa_wta_basin,'~/alfalfa_wta_basin.csv',row.names=FALSE)

# County-level payment system (used to produce Fig 5)
cutoff_wta_county <- data.frame(
  county = c('BOX ELDER','CACHE','RICH','WEBER','DAVIS','MORGAN','SUMMIT','SALT LAKE','UTAH','WASATCH','JUAB'),
  cutoff_wta = c(6093,5573,5284,rep(6093,2),6051,5987,4331,3976,3621,3587),
  cutoff_wta_upper = c(6711,5767,5477,rep(6711,2),6669,6597,4567,4163,3785,3744)
)
alfalfa_wta_county <- merge(alfalfa_wta,cutoff_wta_county,by='county',all.x=T) %>% 
  mutate(enroll = as.numeric(wta_Fallow<=cutoff_wta), # Baseline scenario
         enroll_lower = as.numeric(wta_Fallow_lower<=cutoff_wta), # Lower bound of 95% CI in the baseline scenario
         enroll_upper = as.numeric(wta_Fallow_upper<=cutoff_wta)) %>% # Upper bound of 95% CI in the baseline scenario
  mutate(enroll_a = as.numeric(wta_Fallow<=cutoff_wta_upper), # Conservative scenario
         enroll_a_lower = as.numeric(wta_Fallow_lower<=cutoff_wta_upper), # Lower bound of 95% CI in the conservative scenario
         enroll_a_upper = as.numeric(wta_Fallow_upper<=cutoff_wta_upper)) %>% # Upper bound of 95% CI in the conservative scenario
  dplyr::select(X,Y,basin,county,wta_Fallow,enroll,enroll_lower,enroll_upper,enroll_a,enroll_a_lower,enroll_a_upper)
write.csv(alfalfa_wta_county,'~/alfalfa_wta_county.csv',row.names=FALSE)
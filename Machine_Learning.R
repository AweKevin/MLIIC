
# BlackBoost --------------------------------------------------------------

BlackBoost.mod <- function(train, minsplit = 10,
                           minbucket = 3, maxdepth = 5, mstop = 1000, Seed = 123456) {
  fit <- blackboost(Surv(time, status) ~ ., train,
    family = CoxPH(),
    control = boost_control(mstop = mstop),
    tree_controls = partykit::ctree_control(
      teststat = "quadratic",
      testtype = "Bonferroni",
      mincriterion = 0.95,
      minsplit = minsplit,
      minbucket = minbucket,
      maxdepth = maxdepth,
      saveinfo = FALSE
    )
  )
  set.seed(Seed)
  cvm <- cvrisk(fit, papply = lapply, folds = cv(model.weights(fit), type = "kfold"))
  fit[mstop(cvm), return = F]
  return(fit)
}


# GlmBoost ----------------------------------------------------------------

GlmBoost.mod <- function(train, mstop = 1000, Seed = 123456) {
  fit <- glmboost(Surv(time, status) ~ ., train,
    family = CoxPH(),
    control = boost_control(mstop = mstop),
    center = F
  )
  set.seed(Seed)
  cvm <- cvrisk(fit, papply = lapply, folds = cv(model.weights(fit), type = "kfold"))
  fit[mstop(cvm), return = F]
  return(fit)
}


# Akritas -----------------------------------------------------------------

Akritas.mod <- function(train) {
  fit <- akritas(Surv(time, status) ~ ., train,
    lambda = 0.5, # Bandwidth parameter for uniform smoothing kernel in nearest neighbours estimation.
    # The default value of 0.5 is arbitrary and should be chosen by the user.
    # when lambda = 1, identical to Kaplan-Meier
    reverse = F
  )
  return(fit)
}



# CForest -----------------------------------------------------------------

CForest.mod <- function(train, mstop = 1000, minsplit = 10,
                        minbucket = 3, maxdepth = 5, Seed = 123456) {
  set.seed(Seed)
  fit <- cforest(Surv(time, status) ~ ., train,
    ntree = 1000,
    control = partykit::ctree_control(
      teststat = "quadratic",
      testtype = "Bonferroni",
      mincriterion = 0.95,
      minsplit = minsplit,
      minbucket = minbucket,
      maxdepth = maxdepth,
      saveinfo = FALSE
    )
  )
  return(fit)
}


# CTree -----------------------------------------------------------------

CTree.mod <- function(train, minsplit = 10, minbucket = 3, maxdepth = 5) {
  fit <- ctree(Surv(time, status) ~ ., train,
    control = partykit::ctree_control(
      teststat = "quadratic",
      testtype = "Bonferroni",
      mincriterion = 0.95,
      minsplit = minsplit,
      minbucket = minbucket,
      maxdepth = maxdepth,
      saveinfo = FALSE
    )
  )
  return(fit)
}



# CoxBoost --------------------------------------------------------------

CoxBoost.mod <- function(train, Seed = 123456) {
  time <- train$time
  status <- train$status
  x <- train %>%
    dplyr::select(-c(1:2)) %>%
    as.matrix()
  # determine penalty parameter
  set.seed(Seed)
  optim.CoxBoost <- optimCoxBoostPenalty(
    time = time, status = status, x = x,
    trace = T, start.penalty = 100
  )
  # Fit with obtained penalty parameter and optimal number of boosting
  # steps obtained by cross-validation
  mod.CoxBoost <- CoxBoost(
    time = time, status = status, x = x,
    stepno = optim.CoxBoost$cv.res$optimal.step,
    penalty = optim.CoxBoost$penalty
  )
  return(mod.CoxBoost)
}


# StepwiseCox -----------------------------------------------------------------

StepwiseCox.mod <- function(train) {
  fit <- step(coxph(Surv(time, status) ~ ., train), direction = "backward", trace = F)
  return(fit)
}


# CoxPH -----------------------------------------------------------------

CoxPH.mod <- function(train) {
  fit <- coxph(Surv(time, status) ~ ., train)
  return(fit)
}


# GBM --------------------------------------------------------------

GBM.mod <- function(train, Seed = 123456, maxdepth = 5, minbucket = 3, train.fraction = 0.5) {
  set.seed(Seed)
  mod.GBM <- gbm(Surv(time, status) ~ .,
    distribution = "coxph",
    data = train,
    n.trees = 1000,
    interaction.depth = maxdepth,
    n.minobsinnode = minbucket,
    shrinkage = 0.05,
    bag.fraction = 0.5,
    train.fraction = train.fraction,
    cv.folds = 5, n.cores = 6
  )
  # which.min(mod.GBM$cv.error)
  return(mod.GBM)
}



# SurvReg -----------------------------------------------------------------

SurvReg.mod <- function(train) {
  fit <- survreg(Surv(time, status) ~ ., train)
  return(fit)
}



# ObliqueRSF -----------------------------------------------------------------

ObliqueRSF.mod <- function(train, Seed = 123456) {
  fit <- ORSF(train,
    alpha = 0.5, ntree = 1000, time = "time", status = "status",
    max_pval_to_split_node = 0.05, use.cv = T, verbose = F, random_seed = Seed
  )
  return(fit)
}


# Ranger --------------------------------------------------------------

Ranger.mod <- function(train, Seed = 123456, splitrule = "logrank", # "logrank", "extratrees", "C" or "maxstat" with default "logrank".
                       minbucket = 3, maxdepth = 5) {
  set.seed(Seed)
  mod.Ranger <- ranger(Surv(time, status) ~ ., train,
    num.trees = 1000, min.node.size = minbucket,
    max.depth = maxdepth, seed = Seed, splitrule = splitrule
  )
  return(mod.Ranger)
}


# RSF --------------------------------------------------------------


RSF.mod <- function(train) {
  mod.RSF <- rfsrc(Surv(time, status) ~ .,
    data = train,
    ntree = 1000,
    splitrule = "logrank",
    importance = T,
    proximity = T,
    forest = T,
    seed = 123456
  )
  # mod.RSF <- rfsrc(Surv(time, status) ~ .,
  #                  data = train,
  #                  ntree = which.min(mod.RSF$err.rate),
  #                  splitrule = "logrank",
  #                  importance = T,
  #                  proximity = T,
  #                  forest = T,
  #                  seed = 123456
  # )
  return(mod.RSF)
}


# Rpart --------------------------------------------------------------

Rpart.mod <- function(train, minbucket = 3, maxdepth = 5, minsplit = 10) {
  mod.Rpart <- rpart(Surv(time, status) ~ ., train,
    method = "exp",
    control = rpart.control(
      minsplit = minsplit, minbucket = minbucket,
      maxdepth = maxdepth
    )
  )
  return(mod.Rpart)
}



# SurvivalSVM --------------------------------------------------------------

SurvivalSVM.mod <- function(train) {
  mod.survivalSVM <- survivalsvm(Surv(time, status) ~ .,
    data = train,
    gamma.mu = 0.5,
    type = "regression",
    opt.meth = "ipop", kernel = "add_kernel"
  )
  return(mod.survivalSVM)
}



# Ridge --------------------------------------------------------------

Ridge.mod <- function(train, Seed = 123456) {
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(Seed)
  cv.fit <- cv.glmnet(x2, y2,
    # nfolds = nrow(x2), # 余一交叉验证
    nfolds = 10,
    family = "cox", # cox
    grouped = FALSE, # 在余一交叉验证的时候设置为FALSE
    alpha = 0, # 原文参数（alpha=0为岭回归；=1为lasso；=0~1为弹性网络）
    type.measure = "mse"
  )
  # plot(cv.fit)
  fit <- glmnet(x2, y2, family = "cox", alpha = 0)
  # plot(fit)
  mod.Ridge <- list(cv.fit, fit)
  return(mod.Ridge)
}



# Enet --------------------------------------------------------------

Enet.mod <- function(train, Seed = 123456, alpha = 0.5) {
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(Seed)
  cv.fit <- cv.glmnet(x2, y2,
    # nfolds = nrow(x2), # 余一交叉验证
    nfolds = 10,
    family = "cox", # cox
    grouped = FALSE, # 在余一交叉验证的时候设置为FALSE
    alpha = alpha, # 原文参数（alpha=0为岭回归；=1为lasso；=0~1为弹性网络）
    type.measure = "mse"
  )
  # plot(cv.fit)
  fit <- glmnet(x2, y2, family = "cox", alpha = alpha)
  # plot(fit)
  mod.Enet <- list(cv.fit, fit)
  return(mod.Enet)
}


# Lasso --------------------------------------------------------------

Lasso.mod <- function(train, Seed = 123456) {
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(Seed)
  cv.fit <- cv.glmnet(x2, y2,
    # nfolds = nrow(x2), # 余一交叉验证
    nfolds = 10,
    family = "cox", # cox
    grouped = FALSE, # 在余一交叉验证的时候设置为FALSE
    alpha = 1, # 原文参数（alpha=0为岭回归；=1为lasso；=0~1为弹性网络）
    type.measure = "mse"
  )
  # plot(cv.fit)
  fit <- glmnet(x2, y2, family = "cox", alpha = 1)
  # plot(fit)
  mod.Lasso <- list(cv.fit, fit)
  return(mod.Lasso)
}


# PlsRcox --------------------------------------------------------------

PlsRcox.mod <- function(train, Seed = 123456) {
  time <- train$time
  status <- train$status
  x <- train %>%
    dplyr::select(-c(1:2)) %>%
    as.matrix()
  gene.names <- colnames(x)
  set.seed(Seed)
  plsRcox.cv <- cv.plsRcox(data = list(x = x, time = time, status = status), nfold = 10, nt = 10)
  nt <- which.max(plsRcox.cv$cv.error5[-1])
  mod.plsRcox <- plsRcox(x, time = time, event = status, nt = nt, sparse = F)
  mod.plsRcox <- list(mod.plsRcox, gene.names)
  return(mod.plsRcox)
}


# SuperPC --------------------------------------------------------------

SuperPC.mod <- function(train, Seed = 123456) {
  x2 <- train %>%
    dplyr::select(-c(1:2)) %>%
    t()
  y2 <- train$time
  censoring.status <- train$status
  featurenames <- rownames(x2)
  data.train <- list(
    x = x2,
    y = y2,
    censoring.status = censoring.status,
    featurenames = featurenames
  )
  obj.SuperPC <- superpc.train(data.train, type = "survival")
  set.seed(Seed)
  cv.SuperPC <- superpc.cv(obj.SuperPC,
                           data.train,
                           n.components = 1,
                           n.threshold = 20,
                           min.features = 0,
                           max.features = length(featurenames),
                           n.fold = 10
  )
  sel <- cv.SuperPC$scor %>% which.max()
  threshold <- cv.SuperPC$thresholds[sel]
  mod.SuperPC <- list(obj.SuperPC, threshold, data.train)
  return(mod.SuperPC)
}


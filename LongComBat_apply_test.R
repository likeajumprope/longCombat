## =====================================================================
## Longitudinal ComBat: train on a training set, apply to a test set
## ---------------------------------------------------------------------
## 1. Simulate longitudinal multi-site data with KNOWN site effects
## 2. Split by subject into training and test sets
## 3. Fit longCombat on the training set and store its parameters
## 4. Apply the stored parameters to the test set
## 5. Check: (a) applying to the training set reproduces longCombat
##           (b) site effects are removed in the test set
##           (c) harmonised test data match the simulated "truth"
##
## Requires: lme4, longCombat  (remotes::install_github("jcbeer/longCombat"))
## =====================================================================

install.packages("Matrix")
install.packages("lme4")

library(lme4)
library(Matrix)

install.packages("devtools") 
devtools::load_all("/Users/johannabayer/Dropbox/ENIGMA_site_effects/Site_effects/longCombat")


## ---------------------------------------------------------------------
## 1. Simulate data
##    y_ijv(t) = alpha_v + beta_age*age + beta_sex*sex + eta_j
##               + sigma_v*gamma_iv + sigma_v*delta_iv*eps
##    The "truth" is the same data without site effects (gamma = 0, delta = 1).
## ---------------------------------------------------------------------
simulate_long_sites <- function(n_subj = 400, n_sites = 4, n_features = 10,
                                visits = 2:4, seed = 1) {
  set.seed(seed)
  sites <- paste0("site", seq_len(n_sites))
  
  # subject-level variables
  subj <- data.frame(
    subid    = sprintf("s%04d", seq_len(n_subj)),
    site     = sample(sites, n_subj, replace = TRUE),
    sex      = rbinom(n_subj, 1, 0.5),
    age0     = runif(n_subj, 9, 14),
    n_visits = sample(visits, n_subj, replace = TRUE)
  )
  
  # long format: one row per scan
  dat <- subj[rep(seq_len(n_subj), subj$n_visits), ]
  dat$time <- ave(seq_len(nrow(dat)), dat$subid, FUN = seq_along) - 1   # 0, 1, 2, ...
  dat$age  <- dat$age0 + 1.5 * dat$time + runif(nrow(dat), -0.2, 0.2)
  dat$site <- factor(dat$site, levels = sites)
  rownames(dat) <- NULL
  
  # true parameters per feature
  alpha    <- runif(n_features, 2.5, 3.0)            # intercept (e.g. thickness, mm)
  beta_age <- runif(n_features, -0.04, -0.01)        # per year
  beta_sex <- runif(n_features, -0.05, 0.05)
  sigma    <- runif(n_features, 0.05, 0.08)          # within-subject noise SD
  tau      <- runif(n_features, 0.10, 0.20)          # between-subject SD
  
  # true site effects (site1 = no shift, others shifted and rescaled)
  gamma <- matrix(rnorm(n_sites * n_features, 0, 0.8), n_sites, n_features)
  delta <- matrix(runif(n_sites * n_features, 0.7, 1.6), n_sites, n_features)
  rownames(gamma) <- rownames(delta) <- sites
  
  s_idx <- as.integer(dat$site)
  j_idx <- match(dat$subid, subj$subid)
  
  for (v in seq_len(n_features)) {
    eta  <- rnorm(n_subj, 0, tau[v])[j_idx]
    eps  <- rnorm(nrow(dat), 0, 1)
    bio  <- alpha[v] + beta_age[v] * dat$age + beta_sex[v] * dat$sex + eta
    dat[[paste0("f", v)]]        <- bio + sigma[v] * gamma[s_idx, v] +
      sigma[v] * delta[s_idx, v] * eps        # observed
    dat[[paste0("f", v, "_true")]] <- bio + sigma[v] * eps                 # no site effect
  }
  attr(dat, "truth") <- list(gamma = gamma, delta = delta)
  dat
}

data <- simulate_long_sites()

## Make a training and test split
split_by_subject <- function(dat, prop_train = 0.7, seed = 2) {
  set.seed(seed)
  subj_site <- unique(dat[, c("subid", "site")])
  train_ids <- unlist(lapply(split(subj_site$subid, subj_site$site),
                             function(ids) sample(ids, round(prop_train * length(ids)))))
  list(train = dat[dat$subid %in% train_ids, ],
       test  = dat[!dat$subid %in% train_ids, ])
}

sp    <- split_by_subject(data)    # dat = your long-format data frame
train <- sp$train
test  <- sp$test

# Run local long ComBat
feats <- paste0("f", 1:10)            # measured features only, not the _true columns
train$site <- factor(train$site)      # batch variable must be a factor

res <- longCombat(
  idvar    = "subid",                 # subject ID
  timevar  = "time",                  # visit / time variable (used to label the output)
  batchvar = "site",                  # scanner / site
  features = feats,
  formula  = "age + sex",             # fixed effects to protect; do NOT include site or random effects
  ranef    = "(1|subid)",             # random subject intercept
  data     = train[, c("subid", "time", "site", "age", "sex", feats)]
)

# Use the  parameters on the test set
# See first what we have
res$fixef$f1        # intercept, age, sex and site coefficients for feature f1
res$varcorr$f1      # subject variance and residual variance for f1
res$gammastarhat    # site shifts (sites x features)
res$delta2starhat
mod   <- res$model

## ---- Step 0: choose the feature and get the training parameters ----
f <- "f1" #adjust feature, this can be a loop

fe      <- res$fixef[[f]]                 # fixed effects from training
vc      <- res$varcorr[[f]]               # variance components from training
sigma   <- vc$sdcor[vc$grp == "Residual"] # residual SD - how much do subejct differ on average between each other
tau2    <- vc$vcov[vc$grp == "subid"]     # subject variance - how much do subjects differ on average between each other
gstar   <- res$gammastarhat[, f]          # ComBat site shift (per site)
dstar   <- sqrt(res$delta2starhat[, f])   # ComBat site scaling (per site, as SD)

## ---- Step 1: site effects on the raw scale (as longCombat uses them) ----
sites   <- rownames(res$gammastarhat)     # site1 is the reference level
n_site  <- table(factor(train$site, levels = sites))
b       <- fe[paste0("site", sites[-1])]  # lmer site coefficients (vs. site1)
g1      <- -sum(n_site[-1] * b) / sum(n_site)
site_eff <- setNames(c(g1, b + g1), sites)

## ---- Step 2: expected value m for each test scan (no site) ----
alpha <- fe["(Intercept)"] - site_eff[sites[1]]          # grand mean intercept
m     <- alpha + fe["age"] * test$age + fe["sex"] * test$sex

## ---- Step 3: residual after removing the site shift ----
y      <- test[[f]]
s      <- as.character(test$site)
e      <- y - m - site_eff[s]

## ---- Step 4: each test subject's own level (shrunken average of e) ----
n_j    <- table(test$subid) # how many scans has each subject
sum_e  <- tapply(e, test$subid, sum) #sum e per subject
eta_j  <- sum_e / (n_j + sigma^2 / tau2)   # one value per test subject, this shinks the scan below the average
eta    <- as.numeric(eta_j[as.character(test$subid)])   # back to one per scan

## ---- Step 5: standardise, remove ComBat site effects, transform back ----
z      <- (y - m - eta) / sigma
z_star <- (z - gstar[s]) / dstar[s]
y_star <- sigma * z_star + m + eta

test[[paste0(f, ".combat")]] <- as.numeric(y_star)

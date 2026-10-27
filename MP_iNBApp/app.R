library(shiny)
library(tidyverse)

predict3.1_addMP <- function(year=2010, age=69, screen=0.204, size=11, grade=2, nodes=15,
                             er=0, her2=NA, ki67=0, pr=0, mp=2, probmppos = 0, ######## here mp and probmppos are added
                             horm=0, traz=1, bis=0, radio=1, heart.gy=0, gen=0,
                             smoker=0.15, censor=10) {
  
  #Logical function to determine where a number is an integer or not
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  censor <- censor
  # Vector of follow-up years
  if (is.wholenumber(censor)) { time <- 1:censor } else {time <- c(1:floor(censor), censor)}
  
  # Year of diagnosis set to 2017 if diagnosed after 2017
  year <- ifelse(year>2017, 2017, year)
  
  # Age at diagnosis rescaled
  age <- age - 24
  
  # Set up model coefficients
  ag1_er0 <- 1.7560679383953317
  ag2_er0 <- 4.555274634880487
  sz1_er0 <- 0.7443277265210039
  nd1_er0 <- 0.6308590364310884
  gr1_er0 <- 0.3462287839000938
  sc1_er0 <- -0.211096932659867
  yr1_er0 <- -0.0455046104161727
  ag1_er1 <- 0.1958986203588129
  ag2_er1 <- 2.929125051895787
  sz1_er1 <- 2.27417420673498
  nd1_er1 <- 0.6724828209061207
  gr1_er1 <- 0.7050339629739382
  sc1_er1 <- -0.320393633530154
  yr1_er1 <- -0.0482513707127767
  ag_ot_ea.beta.1 <- 4.211443507724736
  ag_ot_ea.beta.2 <- -31.4120281536527
  yr_ot_ea.beta <- -0.0211864621118188
  h0_br_i <- -3.01531397651658
  h0_br_t1 <- -0.575538052404425
  h0_br_t2 <- -0.102843931717805
  h1_br_i <- -2.31935100392816
  h1_br_t1 <- -3.62253833264139
  h1_br_t2 <- -0.542240944493945
  h_ot_i <- -4.84565475828399
  h_ot_t1 <- 1.341348310005262
  h_ot_t2 <- 0.49539394353046057
  
  # Generate the coefficients depending on values of input parameters
  age.mfp.1   <- ifelse(er==1, ((age/100)^-0.5), (age/100))
  age.beta.1  <- ifelse(er==1, ag1_er1, ag1_er0)
  age.mfp.2   <- ifelse(er==1, ((age/100)^2), (age/100) * log(age/100))
  age.beta.2  <- ifelse(er==1, ag2_er1, ag2_er0)
  size.mfp  <- ifelse(er==1, 1-exp(-size/20), log(size))
  size.beta <- ifelse(er==1, sz1_er1, sz1_er0)
  nodes.mfp   <- ifelse(er==1, log(nodes + 1), log(nodes + 1))
  nodes.beta  <- ifelse(er==1, nd1_er1, nd1_er0)
  grade.beta  <- ifelse(er==1, gr1_er1, gr1_er0)
  screen.beta <- ifelse(er==1, sc1_er1, sc1_er0)
  year.br.beta  <- ifelse(er==1, yr1_er1, yr1_er0)
  ## HER2 time varying effect for ER positive disease
  her2.beta   <- case_when(her2==1 & er==1 ~ c(0.608, .532, .457, .382, .307, .231, .156, .081, .006, 0, 0, 0, 0, 0, 0),
                           her2==0 & er==1 ~ c(-.053, -.046, -.040, -.033, -.027, -.020, -.014, -.007, 0, 0, 0, 0, 0, 0, 0),
                           her2==1 & er==0 ~ 0.2316,
                           her2==0 & er==0 ~ -0.08589,
                           TRUE ~ 0)
  ki67.beta   <- case_when(ki67==1 & er==1 ~ 0.14904,
                           ki67==0 & er==1 ~ -0.11333, 
                           TRUE ~0 ) 
  pr.beta    <- case_when(pr==1 & er==1 ~ -0.0619,
                          pr==0 & er==1 ~ 0.2624,
                          pr==1 & er==0 ~ -0.2231,
                          pr==0 & er==0 ~ 0.0296, 
                          TRUE ~0)
  #########################################
  #########################################
  # Manually added to determine the mp.beta
  mp.beta <- ifelse(mp == 2, 0, 
                    ifelse(mp == 1, log(2.1)*(1-probmppos), 
                           log(2.1)*((1-probmppos)-1)))
  #########################################
  #########################################
  
  ## --- ACM log relative hazard by smoking status ------------------------------------------------
  smoker.beta <- ifelse(smoker==0, log(1/1.025), log(1.219512)) # Assume smoking RR of 2 and 10% population smokers
  
  ## ----baseline_adjust-----------------------------------------------------
  c.other  <- 1.2  # RH non-breast mortality from Kerr et al (2022)
  r.prop   <- 0.69 # Proportion of population receiving radiotherapy
  r.breast <- 0.82 # Relative hazard breast specific mortality from Darby et al
  r.other  <- 1.02 # Relative hazard other disease mortality per Gy heart dose from Taylor et al (2017)
  # RR 1.15 for average heart dose of 6 Gy.
  # r.base.br  <- log(1/((1-r.prop) + r.prop*r.breast))
  # r.base.oth <- log(1/((1-r.prop) + r.prop*(r.other^2))) # Assume 2 Gy average heart dose
  
  ## ----prognostic indices -----
  # Other mortality prognostic index (mi)
  c.oth <- ifelse(gen==0, 0, log(c.other)) # winton fix, for easier comparisons with cljs code
  # c.oth <- log(c.other)
  r.oth <- ifelse(radio==0, 0, log(r.other)*heart.gy)
  mi <-  ag_ot_ea.beta.1*((age/100)^3) + ag_ot_ea.beta.2*((age/100)^3*log(age/100))  + yr_ot_ea.beta*(year-2000) +
    smoker.beta - 0.062
  
  # Breast cancer mortality prognostic index (pi)
  pi.adj <- ifelse(er==1, 0.7, 0.55)
  pi <- age.beta.1*age.mfp.1 + age.beta.2*age.mfp.2 + size.beta*size.mfp +
    nodes.beta*nodes.mfp + grade.beta*grade +
    screen.beta*screen + year.br.beta*(year-2000) +
    her2.beta + ki67.beta + pr.beta + pi.adj + mp.beta ########################### MP: ADDED HERE
  pi <- pi[1:length(time)]
  
  # Treatment coefficients
  c     <- ifelse(gen== 0, 0, ifelse(gen== 2, -.248, -.446))
  h <- if(horm==0) {rep(0, length(time))
  } else { if(horm==5 & er==1) {rep(-.3857, length(time))
  } else { if(horm==10 & er==1 & length(time)<=10) {rep(-.3857, length(time))
  } else{c(rep(-.3857, 10), rep(-.6867, length(time)-10))}
  }
  }
  r     <- ifelse(radio==1, log(r.breast), 0)
  t     <- rep(ifelse(!is.na(her2) & her2==1 & traz==1, -.3567, 0), length(time)) # ADDED MANUALLY: !is.na(her2) & 
  b     <- rep(ifelse(bis==1, -0.198, 0), length(time))
  
  # Treatment combined 
  rx <- tibble(s=rep(0, length(time)), rx = r + h + c + t + b)
  pi.rx <- rx + pi
  mi.rx <- tibble(s = rep(mi, length(time)),
                  mi.rx = rep(mi + c.oth + r.oth, length(time))) #
  
  cols <- ncol(rx)
  
  # Non breast cancer mortality
  # Generate cumulative baseline other mortality
  base.m.cum.oth <- exp(h_ot_i+ h_ot_t1*log(time/10) + h_ot_t2*(time/10)) 
  
  # Generate annual baseline other mortality
  base.m.oth <- base.m.cum.oth
  for (i in 2:length(time)) {
    base.m.oth[i] <- base.m.cum.oth[i] - base.m.cum.oth[i-1] }
  
  # Generate the annual other cancer specific mortality rate
  m.oth <- base.m.oth*exp(mi.rx)
  
  # Calculate the cumulative other cancer mortality rate
  m.cum.oth <- apply(m.oth, 2, cumsum)
  
  # Calculate the cumulative breast cancer survival
  s.cum.oth <- exp(- m.cum.oth)
  m.cum.oth <- 1 - s.cum.oth
  
  m.oth <- m.cum.oth
  for (j in 1:cols) {
    for (i in 2:length(time)) {
      m.oth[i,j] <- m.cum.oth[i,j] - m.cum.oth[i-1,j]
    }
  }
  
  ## ------------------------------------------------------------------------
  # Generate cumulative baseline breast mortality
  base.m.cum.br <- if (er==1) { 
    exp(h1_br_i + h1_br_t1*((time/10)^-0.5) + h1_br_t2*((time/10)^-0.5*log(time/10)))
  } else {exp(h0_br_i + h0_br_t1*((time/10)^-1) + h0_br_t2*((time/10)^-1*log(time/10)))
  }
  
  # Generate annual baseline breast mortality
  base.m.br <- base.m.cum.br
  for (i in 2:length(time)) {
    base.m.br[i] <- base.m.cum.br[i] - base.m.cum.br[i-1] }
  
  # Generate the annual breast cancer specific mortality rate
  m.br <- base.m.br*exp(pi.rx)
  
  # Calculate the cumulative breast cancer mortality rate
  m.cum.br <- apply(m.br, 2, cumsum)
  
  # Calculate the cumulative breast cancer survival
  s.cum.br <- exp(- m.cum.br)
  m.cum.br <- 1- s.cum.br
  
  m.br <- m.cum.br
  for (j in 1:cols) {
    for (i in 2:length(time)) {
      m.br[i,j] <- m.cum.br[i,j] - m.cum.br[i-1,j]
    }
  }
  
  # Cumulative all cause mortality conditional on surviving breast and all cause mortality
  m.cum.all <- 1 - s.cum.oth*s.cum.br
  s.cum.all <- 100 - 100*m.cum.all
  
  # Annual all cause mortality
  m.all <- m.cum.all
  for (j in 1:cols) {
    for (i in 2:length(time)) {
      m.all[i,j] <- m.cum.all[i,j] - m.cum.all[i-1,j]
    }
  }
  
  ## ------------------------------------------------------------------------
  # Proportion of all cause mortality that is breast cancer
  prop.br      <- m.br/(m.br + m.oth)
  pred.m.br    <- prop.br*m.all
  pred.cum.br  <- apply(pred.m.br, 2, cumsum)
  pred.m.oth   <- m.all - pred.m.br
  pred.cum.oth <- apply(pred.m.oth, 2, cumsum)
  pred.cum.all <- pred.cum.br + pred.cum.oth
  
  surv <- as_tibble(100*(1-pred.cum.all))
  
  return(surv)
}
apply_predict <- function(row) {
  # Call PREDICTv3 using values from the row, while keeping fixed values
  result <- predict3.1_addMP(
    year = 2017,
    age = as.numeric(row["age"]),
    smoker = as.numeric(row["smoker"]),
    er = as.numeric(row["er"]),
    pr = ifelse(as.numeric(row["pr"]) == 2, NA, as.numeric(row["pr"])),
    her2 = ifelse(as.numeric(row["her2"]) == 2, NA, as.numeric(row["her2"])),
    ki67 = ifelse(as.numeric(row["ki67"]) == 2, NA, as.numeric(row["ki67"])),
    mp = 2, # meaning unknown
    size = as.numeric(row["size"]),
    grade = as.numeric(row["grade"]),
    screen = as.numeric(row["screen"]),
    nodes = as.numeric(row["nodes"]),
    radio = as.numeric(row["radio"]),
    heart.gy = ifelse(as.numeric(row["radio"]) == 1, 1, 0),
    horm = as.numeric(row["horm"]),
    gen = 0,
    traz = 0, # No impact
    bis = 0,  # No impact
    probmppos = 0, ##### MANUALLY ADDED
    censor = 10  # 10-yr OS
  )
  not = result$rx[10]/100
  result <- predict3.1_addMP(
    year = 2017,
    age = as.numeric(row["age"]),
    smoker = as.numeric(row["smoker"]),
    er = as.numeric(row["er"]),
    pr = ifelse(as.numeric(row["pr"]) == 2, NA, as.numeric(row["pr"])),
    her2 = ifelse(as.numeric(row["her2"]) == 2, NA, as.numeric(row["her2"])),
    ki67 = ifelse(as.numeric(row["ki67"]) == 2, NA, as.numeric(row["ki67"])),
    mp = 2,
    size = as.numeric(row["size"]),
    grade = as.numeric(row["grade"]),
    screen = as.numeric(row["screen"]),
    nodes = as.numeric(row["nodes"]),
    radio = as.numeric(row["radio"]),
    heart.gy = ifelse(as.numeric(row["radio"]) == 1, 1, 0),
    horm = as.numeric(row["horm"]),
    gen = 3,
    traz = 0,
    bis = 0,  # No info
    probmppos = 0, ##### MANUALLY ADDED
    censor = 10  # 10-yr OS
  )
  refscore = s = result$s[10]/100
  t = result$rx[10]/100
  
  probmp1 = as.numeric(linkdf$probs[which.min(abs(linkdf$risks - (1-refscore)))])
  
  no_mp <- c(surv_notreat = s,
             surv_treat_nochemo = not,
             surv_treat_chemo   = t,
             chemo_benefit = (t - not))
  chemo_benefit = (t - not)
  
  # Call PREDICTv3 using values from the row, while keeping fixed values
  result <- predict3.1_addMP(
    year = 2017,
    age = as.numeric(row["age"]),
    smoker = as.numeric(row["smoker"]),
    er = as.numeric(row["er"]),
    pr = ifelse(as.numeric(row["pr"]) == 2, NA, as.numeric(row["pr"])),
    her2 = ifelse(as.numeric(row["her2"]) == 2, NA, as.numeric(row["her2"])),
    ki67 = ifelse(as.numeric(row["ki67"]) == 2, NA, as.numeric(row["ki67"])),
    mp = 0,
    size = as.numeric(row["size"]),
    grade = as.numeric(row["grade"]),
    screen = as.numeric(row["screen"]),
    nodes = as.numeric(row["nodes"]),
    radio = as.numeric(row["radio"]),
    heart.gy = ifelse(as.numeric(row["radio"]) == 1, 1, 0),
    horm = as.numeric(row["horm"]),
    gen = 0,
    traz = 0,
    bis = 0,  # No info
    probmppos = probmp1, ##### MANUALLY ADDED
    censor = 10  # 10-yr OS
  )
  not = result$rx[10]/100
  result <- predict3.1_addMP(
    year = 2017,
    age = as.numeric(row["age"]),
    smoker = as.numeric(row["smoker"]),
    er = as.numeric(row["er"]),
    pr = ifelse(as.numeric(row["pr"]) == 2, NA, as.numeric(row["pr"])),
    her2 = ifelse(as.numeric(row["her2"]) == 2, NA, as.numeric(row["her2"])),
    ki67 = ifelse(as.numeric(row["ki67"]) == 2, NA, as.numeric(row["ki67"])),
    mp = 0,
    size = as.numeric(row["size"]),
    grade = as.numeric(row["grade"]),
    screen = as.numeric(row["screen"]),
    nodes = as.numeric(row["nodes"]),
    radio = as.numeric(row["radio"]),
    heart.gy = ifelse(as.numeric(row["radio"]) == 1, 1, 0),
    horm = as.numeric(row["horm"]),
    gen = 3,
    traz = 0,
    bis = 0,  # No info
    probmppos = probmp1, ##### MANUALLY ADDED
    censor = 10  # 10-yr OS
  )
  s = result$s[10]/100
  t = result$rx[10]/100
  
  mp0 <- c(surv_notreatmp0 = s,
           surv_treat_nochemomp0 = not,
           surv_treat_chemomp0   = t,
           chemo_benefitmp0 = (t - not))
  chemo_benefitmp0 = (t - not)
  
  # Call PREDICTv3 using values from the row, while keeping fixed values
  result <- predict3.1_addMP(
    year = 2017,
    age = as.numeric(row["age"]),
    smoker = as.numeric(row["smoker"]),
    er = as.numeric(row["er"]),
    pr = ifelse(as.numeric(row["pr"]) == 2, NA, as.numeric(row["pr"])),
    her2 = ifelse(as.numeric(row["her2"]) == 2, NA, as.numeric(row["her2"])),
    ki67 = ifelse(as.numeric(row["ki67"]) == 2, NA, as.numeric(row["ki67"])),
    mp = 1,
    size = as.numeric(row["size"]),
    grade = as.numeric(row["grade"]),
    screen = as.numeric(row["screen"]),
    nodes = as.numeric(row["nodes"]),
    radio = as.numeric(row["radio"]),
    heart.gy = ifelse(as.numeric(row["radio"]) == 1, 1, 0),
    horm = as.numeric(row["horm"]),
    gen = 0,
    traz = 0,
    bis = 0,  # No info
    probmppos = probmp1, ##### MANUALLY ADDED
    censor = 10  # 10-yr OS
  )
  not = result$rx[10]/100
  result <- predict3.1_addMP(
    year = 2017,
    age = as.numeric(row["age"]),
    smoker = as.numeric(row["smoker"]),
    er = as.numeric(row["er"]),
    pr = ifelse(as.numeric(row["pr"]) == 2, NA, as.numeric(row["pr"])),
    her2 = ifelse(as.numeric(row["her2"]) == 2, NA, as.numeric(row["her2"])),
    ki67 = ifelse(as.numeric(row["ki67"]) == 2, NA, as.numeric(row["ki67"])),
    mp = 1,
    size = as.numeric(row["size"]),
    grade = as.numeric(row["grade"]),
    screen = as.numeric(row["screen"]),
    nodes = as.numeric(row["nodes"]),
    radio = as.numeric(row["radio"]),
    heart.gy = ifelse(as.numeric(row["radio"]) == 1, 1, 0),
    horm = as.numeric(row["horm"]),
    gen = 3,
    traz = 0,
    bis = 0,  # No info
    probmppos = probmp1, ##### MANUALLY ADDED
    censor = 10  # 10-yr OS
  )
  s = result$s[10]/100
  t = result$rx[10]/100
  
  mp1 <- c(surv_notreatmp1 = s,
           surv_treat_nochemomp1 = not,
           surv_treat_chemomp1   = t,
           chemo_benefitmp1 = (t - not))
  chemo_benefitmp1 = (t - not)
  
  t = as.numeric(row["t_threshold"])
  inbmpclin <- ifelse(chemo_benefitmp1 > t, (chemo_benefitmp1 - t) * probmp1, 0) +
    ifelse(chemo_benefitmp0 > t, (chemo_benefitmp0 - t) * (1-probmp1),0)
  inbclin <- ifelse(chemo_benefit > t, (chemo_benefitmp1 - t) * probmp1, 0) +
    ifelse(chemo_benefit > t, (chemo_benefitmp0 - t) * (1-probmp1),0)
  inbmp <- inbmpclin - inbclin
  
  fin <- c(probmp1 = probmp1,
           inbclin = inbclin,
           inbmpclin   = inbmpclin,
           inbmp = inbmp)
  
  return(c(no_mp, mp0, mp1, fin))
}
probss <- c(0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0003, 0.0005, 0.0008, 0.0012, 0.0018, 0.0024, 0.0033, 0.0043, 0.0055, 0.0069, 0.0084, 0.0102, 0.0121, 0.0142, 0.0166, 0.0190, 0.0217, 0.0245, 0.0275, 0.0307, 0.0339, 0.0374, 0.0409, 0.0445, 0.0483, 0.0521, 0.0561, 0.0601, 0.0642, 0.0684, 0.0726, 0.0769, 0.0812, 0.0856, 0.0900, 0.0945, 0.0990, 0.1035, 0.1081, 0.1126, 0.1172, 0.1218, 0.1264, 0.1310, 0.1357, 0.1403, 0.1449, 0.1495, 0.1542, 0.1588, 0.1634, 0.1680, 0.1726, 0.1772, 0.1818, 0.1864, 0.1909, 0.1955, 0.2000, 0.2045, 0.2090, 0.2135, 0.2180, 0.2224, 0.2269, 0.2313, 0.2357, 0.2401, 0.2445, 0.2488, 0.2531, 0.2574, 0.2617, 0.2660, 0.2702, 0.2744, 0.2787, 0.2828, 0.2870, 0.2911, 0.2953, 0.2994, 0.3034, 0.3075, 0.3115, 0.3155, 0.3195, 0.3235, 0.3274, 0.3314, 0.3353, 0.3392, 0.3430, 0.3469, 0.3507, 0.3545, 0.3583, 0.3620, 0.3658, 0.3695, 0.3732, 0.3769, 0.3805, 0.3842, 0.3878, 0.3914, 0.3950, 0.3985, 0.4020, 0.4056, 0.4091, 0.4125, 0.4160, 0.4194, 0.4228, 0.4262, 0.4296, 0.4330, 0.4363, 0.4396, 0.4430, 0.4462, 0.4495, 0.4528, 0.4560, 0.4592, 0.4624, 0.4656, 0.4687, 0.4719, 0.4750, 0.4781, 0.4812, 0.4842, 0.4873, 0.4903, 0.4933, 0.4963, 0.4993, 0.5023, 0.5052, 0.5082, 0.5111, 0.5140, 0.5169, 0.5197, 0.5226, 0.5254, 0.5283, 0.5311, 0.5338, 0.5366, 0.5394, 0.5421, 0.5448, 0.5475, 0.5502, 0.5529, 0.5556, 0.5582, 0.5609, 0.5635, 0.5661, 0.5687, 0.5713, 0.5738, 0.5764, 0.5789, 0.5814, 0.5839, 0.5864, 0.5889, 0.5914, 0.5938, 0.5963, 0.5987, 0.6011, 0.6035, 0.6059, 0.6082, 0.6106, 0.6129, 0.6153, 0.6176, 0.6199, 0.6222, 0.6245, 0.6267, 0.6290, 0.6312, 0.6334, 0.6357, 0.6379, 0.6400, 0.6422, 0.6444, 0.6465, 0.6487, 0.6508, 0.6529, 0.6550, 0.6571, 0.6592, 0.6613, 0.6634, 0.6654, 0.6674, 0.6695, 0.6715, 0.6735, 0.6755, 0.6775, 0.6794, 0.6814, 0.6833, 0.6853, 0.6872, 0.6891, 0.6910, 0.6929, 0.6948, 0.6967, 0.6985, 0.7004, 0.7022, 0.7041, 0.7059, 0.7077, 0.7095, 0.7113, 0.7131, 0.7149, 0.7166, 0.7184, 0.7201, 0.7219, 0.7236, 0.7253, 0.7270, 0.7287, 0.7304, 0.7321, 0.7337, 0.7354, 0.7370, 0.7387, 0.7403, 0.7419, 0.7435, 0.7451, 0.7467, 0.7483, 0.7499, 0.7515, 0.7530, 0.7546, 0.7561, 0.7576, 0.7592, 0.7607, 0.7622, 0.7637, 0.7652, 0.7667, 0.7681, 0.7696, 0.7710, 0.7725, 0.7739, 0.7754, 0.7768, 0.7782, 0.7796, 0.7810, 0.7824, 0.7838, 0.7852, 0.7866, 0.7879, 0.7893, 0.7906, 0.7920, 0.7933, 0.7946, 0.7959, 0.7972, 0.7985, 0.7998, 0.8011, 0.8024, 0.8037, 0.8049, 0.8062, 0.8074, 0.8087, 0.8099, 0.8111, 0.8124, 0.8136, 0.8148, 0.8160, 0.8172, 0.8184, 0.8196, 0.8207, 0.8219, 0.8231, 0.8242, 0.8254, 0.8265, 0.8276, 0.8288, 0.8299, 0.8310, 0.8321, 0.8332, 0.8343, 0.8354, 0.8365, 0.8376, 0.8386, 0.8397, 0.8408, 0.8418, 0.8429, 0.8439, 
            0.8449, 0.8460, 0.8470, 0.8480, 0.8490, 0.8500, 0.8510, 0.8520, 0.8530, 0.8540, 0.8550, 0.8559, 0.8569, 0.8578, 0.8588, 0.8598, 0.8607, 0.8616, 0.8626, 0.8635, 0.8644, 0.8653, 0.8662, 0.8671, 0.8681, 0.8689, 0.8698, 0.8707, 0.8716, 0.8725, 0.8733, 0.8742, 0.8750, 0.8759, 0.8767, 0.8776, 0.8784, 0.8793, 0.8801, 0.8809, 0.8817, 0.8825, 0.8833, 0.8841, 0.8849, 0.8857, 0.8865, 0.8873, 0.8881, 0.8889, 0.8896, 0.8904, 0.8912, 0.8919, 0.8927, 0.8934, 0.8941, 0.8949, 0.8956, 0.8964, 0.8971, 0.8978, 0.8985, 0.8992, 0.8999, 0.9006, 0.9013, 0.9020, 0.9027, 0.9034, 0.9041, 0.9048, 0.9054, 0.9061, 0.9068, 0.9074, 0.9081, 0.9088, 0.9094, 0.9100, 0.9107, 0.9113, 0.9119, 0.9126, 0.9132, 0.9138, 0.9145, 0.9151, 0.9157, 0.9163, 0.9169, 0.9175, 0.9181, 0.9187, 0.9192, 0.9198, 0.9204, 0.9210, 0.9215, 0.9221, 0.9227, 0.9232, 0.9238, 0.9244, 0.9249, 0.9255, 0.9260, 0.9266, 0.9271, 0.9277, 0.9282, 0.9287, 0.9292, 0.9298, 0.9303, 0.9308, 0.9313, 0.9318, 0.9323, 0.9328, 0.9333, 0.9338, 0.9343, 0.9348, 0.9352, 0.9357, 0.9362, 0.9367, 0.9372, 0.9377, 0.9381, 0.9386, 0.9390, 0.9395, 0.9399, 0.9404, 0.9408, 0.9413, 0.9417, 0.9422, 0.9426, 0.9431, 0.9435, 0.9439, 0.9443, 0.9448, 0.9452, 0.9456, 0.9460, 0.9464, 0.9468, 0.9473, 0.9477, 0.9481, 0.9485, 0.9488, 0.9493, 0.9496, 0.9501, 0.9504, 0.9508, 0.9512, 0.9516, 0.9519, 0.9523, 0.9527, 0.9530, 0.9535, 0.9538, 0.9541, 0.9545, 0.9549, 0.9552, 0.9555, 0.9559, 0.9563, 0.9566, 0.9570, 0.9572, 0.9576, 0.9579, 0.9583, 0.9586, 0.9589, 0.9593, 0.9596, 0.9599, 0.9602, 0.9606, 0.9609, 0.9612, 0.9615, 0.9618, 0.9621, 0.9625, 0.9627, 0.9630, 0.9633, 0.9637, 0.9639, 0.9642, 0.9645, 0.9648, 0.9651, 0.9654, 0.9656, 0.9659, 0.9663, 0.9664, 0.9667, 0.9671, 0.9674, 0.9676, 0.9678, 0.9681, 0.9684, 0.9686, 0.9689, 0.9692, 0.9694, 0.9697, 0.9699, 0.9702, 0.9705, 0.9707, 0.9708, 0.9711, 0.9714, 0.9716, 0.9718, 0.9721, 0.9724, 0.9726, 0.9729, 0.9730, 0.9733, 0.9735, 0.9736, 0.9740, 0.9742, 0.9744, 0.9747, 0.9748, 0.9751, 0.9753, 0.9755, 0.9757, 0.9759, 0.9761, 0.9763, 0.9766, 0.9767, 0.9770, 0.9770, 0.9773, 0.9775, 0.9777, 0.9779, 0.9781, 0.9783, 0.9784, 0.9786, 0.9788, 0.9789, 0.9792, 0.9794, 0.9796, 0.9798, 0.9798, 0.9801, 0.9803, 0.9805, 0.9808, 0.9808, 0.9810, 0.9812, 0.9813, 0.9815, 0.9816, 0.9819, 0.9820, 0.9821, 0.9823, 0.9824, 0.9826, 0.9828, 0.9831, 0.9831, 0.9832, 0.9834, 0.9836, 0.9837, 0.9837, 0.9842, 0.9842, 0.9843, 0.9845, 0.9846, 0.9848, 0.9848, 0.9850, 0.9851, 0.9853, 0.9854, 0.9856, 0.9856, 0.9858, 0.9860, 0.9860, 0.9863, 0.9863, 0.9866, 0.9866, 0.9867, 0.9868, 0.9869, 0.9869, 0.9873, 0.9873, 0.9873, 0.9878, 0.9878, 0.9878, 0.9878, 0.9883, 0.9883, 0.9883, 0.9883, 0.9886, 0.9887, 0.9887,
            0.9890, 0.9890, 0.9891, 0.9891, 0.9891, 0.9891, 0.9891, 0.9901, 0.9901, 0.9901, 0.9901, 0.9901, 0.9901, 0.9904, 0.9904, 0.9904, 0.9906, 0.9906, 0.9907, 0.9907, 0.9907, 0.9911, 0.9911, 0.9911, 0.9911, 0.9915, 0.9915, 0.9915, 0.9915, 0.9915, 0.9920, 0.9920, 0.9920, 0.9920, 0.9920, 0.9923, 0.9923, 0.9923, 0.9923, 0.9923, 0.9923, 0.9923, 0.9923, 0.9923, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9935, 0.9940, 0.9940, 0.9940, 0.9940, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9941, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960, 0.9960)
riskss <- c(0.0086, 0.0086, 0.0086, 0.0086, 0.0086, 0.0086, 0.0086, 0.0086, 0.0086, 0.0089, 0.0100, 0.0110, 0.0120, 0.0130, 0.0140, 0.0150, 0.0160, 0.0170, 0.0180, 0.0190, 0.0200, 0.0210, 0.0220, 0.0230, 0.0240, 0.0250, 0.0260, 0.0270, 0.0280, 0.0290, 0.0300, 0.0310, 0.0320, 0.0330, 0.0340, 0.0350, 0.0360, 0.0370, 0.0380, 0.0390, 0.0400, 0.0410, 0.0420, 0.0430, 0.0440, 0.0450, 0.0460, 0.0470, 0.0480, 0.0490, 0.0500, 0.0510, 0.0520, 0.0530, 0.0540, 0.0550, 0.0560, 0.0570, 0.0580, 0.0590, 0.0600, 0.0610, 0.0620, 0.0630, 0.0640, 0.0650, 0.0660, 0.0670, 0.0680, 0.0690, 0.0700, 0.0710, 0.0720, 0.0730, 0.0740, 0.0750, 0.0760, 0.0770, 0.0780, 0.0790, 0.0800, 0.0810, 0.0820, 0.0830, 0.0840, 0.0850, 0.0860, 0.0870, 0.0880, 0.0890, 0.0900, 0.0910, 0.0920, 0.0930, 0.0940, 0.0950, 0.0960, 0.0970, 0.0980, 0.0990, 0.1000, 0.1010, 0.1020, 0.1030, 0.1040, 0.1050, 0.1060, 0.1070, 0.1080, 0.1090, 0.1100, 0.1110, 0.1120, 0.1130, 0.1140, 0.1150, 0.1160, 0.1170, 0.1180, 0.1190, 0.1200, 0.1210, 0.1220, 0.1230, 0.1240, 0.1250, 0.1260, 0.1270, 0.1280, 0.1290, 0.1300, 0.1310, 0.1320, 0.1330, 0.1340, 0.1350, 0.1360, 0.1370, 0.1380, 0.1390, 0.1400, 0.1410, 0.1420, 0.1430, 0.1440, 0.1450, 0.1460, 0.1470, 0.1480, 0.1490, 0.1500, 0.1510, 0.1520, 0.1530, 0.1540, 0.1550, 0.1560, 0.1570, 0.1580, 0.1590, 0.1600, 0.1610, 0.1620, 0.1630, 0.1640, 0.1650, 0.1660, 0.1670, 0.1680, 0.1690, 0.1700, 0.1710, 0.1720, 0.1730, 0.1740, 0.1750, 0.1760, 0.1770, 0.1780, 0.1790, 0.1800, 0.1810, 0.1820, 0.1830, 0.1840, 0.1850, 0.1860, 0.1870, 0.1880, 0.1890, 0.1900, 0.1910, 0.1920, 0.1930, 0.1940, 0.1950, 0.1960, 0.1970, 0.1980, 0.1990, 0.2000, 0.2010, 0.2020, 0.2030, 0.2040, 0.2050, 0.2060, 0.2070, 0.2080, 0.2090, 0.2100, 0.2110, 0.2120, 0.2130, 0.2140, 0.2150, 0.2160, 0.2170, 0.2180, 0.2190, 0.2200, 0.2210, 0.2220, 0.2230, 0.2240, 0.2250, 0.2260, 0.2270, 0.2280, 0.2290, 0.2300, 0.2310, 0.2320, 0.2330, 0.2340, 0.2350, 0.2360, 0.2370, 0.2380, 0.2390, 0.2400, 0.2410, 0.2420, 0.2430, 0.2440, 0.2450, 0.2460, 0.2470, 0.2480, 0.2490, 0.2500, 0.2510, 0.2520, 0.2530, 0.2540, 0.2550, 0.2560, 0.2570, 0.2580, 0.2590, 0.2600, 0.2610, 0.2620, 0.2630, 0.2640, 0.2650, 0.2660, 0.2670, 0.2680, 0.2690, 0.2700, 0.2710, 0.2720, 0.2730, 0.2740, 0.2750, 0.2760, 0.2770, 0.2780, 0.2790, 0.2800, 0.2810, 0.2820, 0.2830, 0.2840, 0.2850, 0.2860, 0.2870, 0.2880, 0.2890, 0.2900, 0.2910, 0.2920, 0.2930, 0.2940, 0.2950, 0.2960, 0.2970, 0.2980, 0.2990, 0.3000, 0.3010, 0.3020, 0.3030, 0.3040, 0.3050, 0.3060, 0.3070, 0.3080, 0.3090, 0.3100, 0.3110, 0.3120, 0.3130, 0.3140, 0.3150, 0.3160, 0.3170, 0.3180, 0.3190, 0.3200, 0.3210, 0.3220, 0.3230, 0.3240, 0.3250, 0.3260, 0.3270, 0.3280, 0.3290, 0.3300, 0.3310, 0.3320, 0.3330, 0.3340, 
            0.3350, 0.3360, 0.3370, 0.3380, 0.3390, 0.3400, 0.3410, 0.3420, 0.3430, 0.3440, 0.3450, 0.3460, 0.3470, 0.3480, 0.3490, 0.3500, 0.3510, 0.3520, 0.3530, 0.3540, 0.3550, 0.3560, 0.3570, 0.3580, 0.3590, 0.3600, 0.3610, 0.3620, 0.3630, 0.3640, 0.3650, 0.3660, 0.3670, 0.3680, 0.3690, 0.3700, 0.3710, 0.3720, 0.3730, 0.3740, 0.3750, 0.3760, 0.3770, 0.3780, 0.3790, 0.3800, 0.3810, 0.3820, 0.3830, 0.3840, 0.3850, 0.3860, 0.3870, 0.3880, 0.3890, 0.3900, 0.3910, 0.3920, 0.3930, 0.3940, 0.3950, 0.3960, 0.3970, 0.3980, 0.3990, 0.4000, 0.4010, 0.4020, 0.4030, 0.4040, 0.4050, 0.4060, 0.4070, 0.4080, 0.4090, 0.4100, 0.4110, 0.4120, 0.4130, 0.4140, 0.4150, 0.4160, 0.4170, 0.4180, 0.4190, 0.4200, 0.4210, 0.4220, 0.4230, 0.4240, 0.4250, 0.4260, 0.4270, 0.4280, 0.4290, 0.4300, 0.4310, 0.4320, 0.4329, 0.4340, 0.4350, 0.4360, 0.4370, 0.4380, 0.4390, 0.4400, 0.4410, 0.4420, 0.4430, 0.4441, 0.4450, 0.4460, 0.4470, 0.4480, 0.4490, 0.4500, 0.4510, 0.4520, 0.4530, 0.4540, 0.4550, 0.4561, 0.4570, 0.4580, 0.4590, 0.4600, 0.4610, 0.4620, 0.4630, 0.4640, 0.4650, 0.4660, 0.4670, 0.4679, 0.4690, 0.4699, 0.4710, 0.4720, 0.4730, 0.4740, 0.4750, 0.4760, 0.4770, 0.4780, 0.4790, 0.4800, 0.4810, 0.4820, 0.4830, 0.4840, 0.4850, 0.4860, 0.4870, 0.4880, 0.4890, 0.4899, 0.4911, 0.4920, 0.4931, 0.4940, 0.4950, 0.4960, 0.4970, 0.4980, 0.4990, 0.5000, 0.5009, 0.5021, 0.5031, 0.5039, 0.5051, 0.5061, 0.5069, 0.5079, 0.5090, 0.5100, 0.5110, 0.5120, 0.5129, 0.5139, 0.5149, 0.5160, 0.5170, 0.5178, 0.5190, 0.5201, 0.5210, 0.5220, 0.5230, 0.5240, 0.5249, 0.5259, 0.5271, 0.5280, 0.5291, 0.5300, 0.5308, 0.5320, 0.5331, 0.5340, 0.5350, 0.5360, 0.5370, 0.5380, 0.5390, 0.5400, 0.5409, 0.5422, 0.5428, 0.5440, 0.5451, 0.5463, 0.5470, 0.5480, 0.5490, 0.5501, 0.5508, 0.5520, 0.5530, 0.5539, 0.5551, 0.5560, 0.5571, 0.5582, 0.5590, 0.5596, 0.5610, 0.5621, 0.5629, 0.5637, 0.5651, 0.5662, 0.5670, 0.5686, 0.5686, 0.5700, 0.5710, 0.5716, 0.5731, 0.5741, 0.5750, 0.5762, 0.5769, 0.5781, 0.5791, 0.5800, 0.5809, 0.5819, 0.5829, 0.5838, 0.5855, 0.5861, 0.5875, 0.5876, 0.5891, 0.5898, 0.5910, 0.5920, 0.5931, 0.5941, 0.5945, 0.5959, 0.5969, 0.5975, 0.5991, 0.5999, 0.6009, 0.6023, 0.6023, 0.6040, 0.6050, 0.6060, 0.6079, 0.6079, 0.6089, 0.6100, 0.6106, 0.6120, 0.6128, 0.6147, 0.6150, 0.6158, 0.6169, 0.6178, 0.6189, 0.6199, 0.6217, 0.6219, 0.6227, 0.6241, 0.6253, 0.6258, 0.6258, 0.6292, 0.6292, 0.6299, 0.6316, 0.6319, 0.6331, 0.6331, 0.6351, 0.6358, 0.6367, 0.6377, 0.6395, 0.6395, 0.6409, 0.6420, 0.6425, 0.6443, 0.6447, 0.6470, 0.6470, 0.6479, 0.6487, 0.6494, 0.6494, 0.6527, 0.6527, 0.6527, 0.6565, 0.6565, 0.6565, 0.6565, 0.6608, 0.6608, 0.6608, 0.6613, 0.6637, 0.6642, 0.6642, 
            0.6668, 0.6668, 0.6685, 0.6685, 0.6685, 0.6685, 0.6685, 0.6774, 0.6774, 0.6774, 0.6774, 0.6774, 0.6774, 0.6805, 0.6805, 0.6805, 0.6830, 0.6830, 0.6841, 0.6841, 0.6841, 0.6884, 0.6884, 0.6884, 0.6884, 0.6929, 0.6929, 0.6929, 0.6929, 0.6929, 0.6980, 0.6980, 0.6980, 0.6980, 0.6980, 0.7025, 0.7025, 0.7025, 0.7025, 0.7025, 0.7025, 0.7025, 0.7025, 0.7025, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7174, 0.7247, 0.7247, 0.7247, 0.7247, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7264, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603, 0.7603)            
linkdf <- data.frame(risks = riskss,probs = probss)            

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      
      .panel-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(350px, 1fr));
        gap: 20px;
        align-items: start;
      }

      .panel-box {
        width: 100%;
      }

      .panel-title {
        text-align: center;
        margin-bottom: 10px;
      }

      /* Header styling */
      .app-header {
        text-align: center;
        padding: 18px 10px 10px 10px;
        border-bottom: 1px solid #ddd;
        margin-bottom: 15px;
      }

      .app-header h1 {
        margin-bottom: 6px;
        font-weight: 600;
      }

      .app-header h4 {
        margin-top: 0px;
        font-weight: 400;
        color: #555;
      }

      .app-header p {
        font-size: 12px;
        color: #777;
        margin-top: 8px;
      }

      /* Phone layout */
      @media (max-width: 768px) {
        
        .panel-grid {
          grid-template-columns: 1fr;
        }
      }
    "))
  ),
  
  # ---- NEW HEADER (replaces titlePanel) ----
  tags$div(
    class = "app-header",
    
    tags$h1("MammaPrint Decision Support Dashboard"),
    
    tags$h4("Extension to Breast PREDICT to quantify the value of MammaPrint testing"),
    
    tags$p("Decision support tool only — it can be assistive but does not replace clinical judgement"),
    
    tags$details(
      tags$summary("Disclaimer and Technical details"),
      
      tags$p(
        "This tool is an extension of Breast PREDICT and does not replace its use."
      ),
      tags$p(
        "It assumes that MammaPrint updates the breast cancer–specific mortality estimate from Breast PREDICT using a hazard ratio of 2.1."
      ),
      tags$p(
        "The probability of a high-risk MammaPrint result is estimated based on a correlation of 0.45 with clinical risk characteristics."
      ),
      tags$p(
        "These assumptions can be modified in future versions to improve individualized Net Bnefit estimates."
      )
    )
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4("Patient characteristics in Breast PREDICT"),
      
      numericInput("age", "Age at diagnosis (25 to 85)", value = 51, min = 25, max = 85),
      
      selectInput("smoker", "Smoker",
                  c("Never or Ex" = 0, "Yes" = 1)),
      
      selectInput("er", "ER status",
                  c("Negative" = 0, "Positive" = 1)),
      
      selectInput("pr", "Progesterone status",
                  c("Negative" = 0, "Positive" = 1, "Unknown" = 2)),
      
      selectInput("her2", "HER2/ERRB2 status",
                  c("Negative" = 0, "Positive" = 1, "Unknown" = 2)),
      
      selectInput("ki67", "Ki-67 status",
                  c("Negative" = 0, "Positive" = 1, "Unknown" = 2)),
      
      numericInput("size",
                   "Invasive tumour size (mm)",
                   value = 18,
                   min = 1,
                   max = 100,
                   step = 1),
      
      selectInput("grade",
                  "Tumour grade",
                  c(1, 2, 3)),
      
      selectInput("screen", "Detected by",
                  c("Symptoms" = 0, "Screening" = 1, "Unknown" = 0.204)),
      
      numericInput("nodes",
                   "Positive nodes",
                   min = 1,
                   max = 100,
                   value = 1,
                   step = 1),
      
      selectInput("radio", "Radiotherapy",
                  c("No" = 0, "Yes" = 1)),
      
      selectInput("horm",
                  "Hormone therapy",
                  c(0, 5, 10)),
      
      tags$hr(),
      
      h4("Patient specific treatment threshold"),
      
      helpText("Minimum benefit required to accept chemotherapy"),
      
      sliderInput("t_threshold",
                  "Treatment threshold",
                  min = 0,
                  max = 10,
                  value = 3,
                  step = 0.1,
                  post = "%")
    ),
    
    mainPanel(
      h3(textOutput("decision")),
      hr(),
      
      div(
        class = "panel-grid",
        
        div(
          class = "panel-box",
          h4("Chemotherapy benefit estimated by Breast PREDICT",
             class = "panel-title"),
          plotOutput("base_plot", height = "520px")
        ),
        
        div(
          class = "panel-box",
          h4("Probability high risk MammaPrint",
             class = "panel-title"),
          plotOutput("prob", height = "520px")
        ),
        
        div(
          class = "panel-box",
          h4("Updating benefit with MammaPrint",
             class = "panel-title"),
          plotOutput("mp_plot", height = "520px")
        ),
        
        div(
          class = "panel-box",
          h4("The iNB of adding MammaPrint",
             class = "panel-title"),
          plotOutput("delta", height = "520px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ---- reactive model ----
  res <- reactive({
    
    patient <- data.frame(
      age = input$age,
      smoker = as.numeric(input$smoker),
      er = as.numeric(input$er),
      pr = as.numeric(input$pr),
      her2 = as.numeric(input$her2),
      ki67 = as.numeric(input$ki67),
      size = input$size,
      grade = as.numeric(input$grade),
      screen = as.numeric(input$screen),
      nodes = input$nodes,
      radio = as.numeric(input$radio),
      horm = input$horm,
      censor = 10,
      t_threshold = (input$t_threshold/100)
    )
    
    apply_predict(patient)
  })
  
  # ---- decision label ----
  output$decision <- renderText({
    
    r <- res()
    inb_mp <- round(r[["inbmp"]]*100,2)
    
    if (inb_mp > 0) {
      paste0(
        "Consider MammaPrint \u27A1 MammaPrint is expected to improve the treatment decision ",
        "(Net survival Benefit (iNB) = ", round(inb_mp, 4),"%)")
    } else {
      "Consider not doing MammaPrint \u27A1 MammaPrint is not expected to improve the treatment decision (no gain in iNB)"
    }
  })
  
  output$base_plot <- renderPlot({
    
    benefit <- res()[["chemo_benefit"]]
    
    plot(0, 0,
         type = "n",
         xlim = c(0, 1),
         ylim = c(0, 1),
         axes = FALSE,
         xlab = "",
         ylab = "")
    
    # wider rectangle
    rect(0.0, 0.28, 1, 0.72,
         col = "grey95",
         border = "grey40",
         lwd = 3)
    
    # title
    text(0.5, 0.58,
         "Estimated chemotherapy\nbenefit by\nBreast PREDICT",
         cex = 1.6,
         font = 2)
    
    # big number
    text(0.5, 0.40,
         paste0(round(benefit * 100, 1), "%"),
         cex = 3.2,
         font = 2)
  })
  
  # ---- probability ----
  output$prob <- renderPlot({
    p = res()[["probmp1"]]
    plot(1, xlim = c(0, 1), ylim = c(0, 1),
         type = "n", axes = FALSE, xlab = "", ylab = "")
    # background bar
    rect(0, 0, 1, 1, col = "grey90", border = NA)
    # filled probability
    rect(0, 0, 1, p, col = "steelblue", border = NA)
    # left axis (0 to 1 by 0.1)
    axis(2, at = seq(0, 1, by = 0.1), las = 1)
    title(main = paste0(round(p * 100, 1), "% probability of high-risk result"))
  })
  
  # ---- MP effects ----
  output$mp_plot <- renderPlot({
    
    r <- res()
    
    base <- r[["chemo_benefit"]]
    mp0  <- r[["chemo_benefitmp0"]]
    mp1  <- r[["chemo_benefitmp1"]]
    p    <- r[["probmp1"]]
    
    plot(0, 0,
         type = "n",
         xlim = c(0, 1),
         ylim = c(0, 1),
         axes = FALSE,
         xlab = "",
         ylab = "")
    
    # --------------------------
    # Baseline card (top)
    # --------------------------
    rect(0.18, 0.68, 0.82, 0.92,
         col = "grey95",
         border = "grey40",
         lwd = 3)
    
    text(0.5, 0.86,
         "Breast PREDICT",
         cex = 1.6,
         font = 2)
    
    text(0.5, 0.75,
         paste0(round(base * 100, 1), "%"),
         cex = 2.4,
         font = 2)
    
    # --------------------------
    # Arrows
    # --------------------------
    arrows(0.42, 0.68, 0.24, 0.46,
           length = 0.10,
           lwd = 3)
    
    arrows(0.58, 0.68, 0.76, 0.46,
           length = 0.10,
           lwd = 3)
    
    # --------------------------
    # Probability labels
    # --------------------------
    text(0.18, 0.54,
         paste0(
           "Low-risk\n",
           round((1 - p) * 100, 1),
           "%"
         ),
         cex = 1.3,
         font = 2)
    
    text(0.82, 0.54,
         paste0(
           "High-risk\n",
           round(p * 100, 1),
           "%"
         ),
         cex = 1.3,
         font = 2)
    
    # --------------------------
    # Low-risk card
    # --------------------------
    rect(0.04, 0.12, 0.44, 0.38,
         col = "grey95",
         border = "grey40",
         lwd = 3)
    
    text(0.24, 0.31,
         "Low-risk MP",
         cex = 1.3,
         font = 2)
    
    text(0.24, 0.20,
         paste0(round(mp0 * 100, 1), "%"),
         cex = 3,
         font = 2)
    
    # --------------------------
    # High-risk card
    # --------------------------
    rect(0.56, 0.12, 0.96, 0.38,
         col = "grey95",
         border = "grey40",
         lwd = 3)
    
    text(0.76, 0.31,
         "High-risk MP",
         cex = 1.3,
         font = 2)
    
    text(0.76, 0.20,
         paste0(round(mp1 * 100, 1), "%"),
         cex = 3,
         font = 2)
  })
  
  # ---- iNB decomposition (UPDATED CORRECTLY) ----
  output$delta <- renderPlot({
    
    r <- res()
    inb_mp <- round(r[["inbmp"]]*100,2)
    
    plot(1, xlim = c(0, 1), ylim = c(-0.0001, max(1, inb_mp)),
         type = "n", axes = FALSE, xlab = "", ylab = "")
    
    rect(-0.0001, -0.0001, 1, 1, col = "grey90", border = NA)
    
    col <- ifelse(inb_mp > 0, "darkgreen", "red")
    
    rect(-0.0001, -0.0001, 1,  max(inb_mp,0), col = col, border = NA)
    
    axis(2, las = 1)
    abline(h = 0, lty = 2)
    
    title(main = paste0("iNB of MammaPrint = ", round(inb_mp, 4)))
  })
}

source("ui.R")
source("server.R")
source("helpers.R")

shinyApp(ui = ui, server = server)
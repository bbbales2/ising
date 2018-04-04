library(MASS)

source("casm_helper.R")
source("ising_helpers.R")
require(Rcpp)

library(reshape2)
library(tidyverse)
library(ggplot2)
library(rstan)
library(stats)
library(parallel)
library(gtools)
library(GGally)
library(purrr)

path = "/home/bbales2/casm/invent2"
ecis = rep(0, length(getECIs(path)))
N = getSupercell(path)
keep = 2:13
nonZero = c(3, 4, 5, 6, 7, 14, 15, 16, 17, 18)
ecis[nonZero] = c(0.440, 0.280, 0.100, -0.100, -0.133, 0.010, -0.180, 0.170, 0.080, -0.080)
nonZero = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)

paths = c("/home/bbales2/ising/paper_outputs/test1.dat",
          "/home/bbales2/ising/paper_outputs/test2.dat",
          "/home/bbales2/ising/paper_outputs/test3.dat",
          "/home/bbales2/ising/paper_outputs/test4.dat",
          "/home/bbales2/ising/paper_outputs/test5.dat")

env = new.env()
load(paths[1], envir = env)

getData = function(data, vars) {
  data = map(paths, function(path) {
    env = new.env()
    load(path, envir = env)
    out = list()
      
    for(i in 1:length(vars)) {
      if(class(env[[vars[i]]])[[1]] == "list") {
        out[[vars[i]]] = do.call(rbind, env[[vars[i]]]) %>%
          as.tibble %>%
          mutate(which = basename(path),
                 opt = row_number())
      } else {
        out[[vars[i]]] = env[[vars[i]]] %>%
          mutate(which = basename(path))
      }
    }
    
    out
  })
  
  outT = list()
  for(i in 1:length(vars)) {
    outT[[vars[i]]] = map(data, ~ .[[vars[i]]]) %>%
      do.call(rbind, .)
  }
  
  outT
}

smu = function(corrs) {
  res = corrs %>% select(starts_with("corr")) %>% as.matrix
  
  f = matrix(0, nrow = ncol(res), ncol = 1)
  rownames(f) = colnames(res)
  for(i in 1:ncol(res)) {
    f[i, 1] = mean(res[, i])
  }
  
  ## I'm not sure why I need the two transposes on Eg to get things to work, but seems like I do
  list(Eg = t(t(f[keep, 1])))
}

runSimulation = function(g) {
  setECIs(path, g)
  runMC(path)
  corrs = getCorrs(path)
  
  corrs %>%
    group_by(mu) %>%
    filter(mci > 500) %>%
    do(out = smu(.)) %>% pull(out)
}

data = getData(paths, c('opts', 'opts2', 'clexes', 'tclex', 'crs', 'tcr'))

# Compare results of optimization to true ecis
levels = names(data$opts2)[nonZero]
data$opts2 %>%
  #filter(which %in% c("test1.dat", "test2.dat", "test3.dat")) %>%
  select(nonZero, which) %>%
  bind_rows(ecis[nonZero] %>%
        setNames(levels) %>%
        t %>%
        as.tibble %>%
        mutate(which = "TRUTH")) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = levels)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 2, stroke = 2, position = position_dodge(width = 0.75)) +
  theme_set(theme_gray(base_size = 18))

i = 1

runSimDataToDf = function(data) {
  do.call(cbind, map(data, ~ .$Eg)) %>%
    as.tibble %>%
    setNames(getChemicalPotentials(path)) %>%
    mutate(corr = as.factor(row_number())) %>%
    gather(mu, value, -corr) %>%
    mutate(mu = as.numeric(mu))
}

pData = list()
for(i in 1:nrow(data$opts2)) {
  setTemperatureFraction(path, 1.0)
  datt1 = runSimulation(data$opts2[i,] %>% select(starts_with("corr")) %>% unlist())
  setTemperatureFraction(path, 0.5)
  datt2 = runSimulation(data$opts2[i,] %>% select(starts_with("corr")) %>% unlist())

  pData[[i]] = bind_rows(runSimDataToDf(datt1) %>%
                           mutate(temp = 1.0),
                         runSimDataToDf(datt2) %>%
                           mutate(temp = 0.5)) %>%
    mutate(which = data$opts2[i,]$which,
           opt = data$opts2[i,]$opt)
  
  cat("Finished ", i, " out of ", nrow(data$opts2), "\n")
}

pData2 = map(pData, function(data) { data %>% group_by(corr, temp) %>% mutate(mu = getChemicalPotentials(path)) %>% ungroup} )

setTemperatureFraction(path, 1.0)
datt1 = runSimulation(ecis)
setTemperatureFraction(path, 0.5)
datt2 = runSimulation(ecis)

tpData = bind_rows(runSimDataToDf(datt1) %>%
                     mutate(temp = 1.0),
                   runSimDataToDf(datt2) %>%
                     mutate(temp = 0.5)) %>%
  mutate(which = "truth")

lps_filtered = (left_join(pData %>%
                            bind_rows, tpData, by = c("corr", "mu", "temp")) %>%
  group_by(which.x, opt) %>%
  summarize(lp = sum((value.x - value.y)^2)) -> lps) %>%
  top_n(10, -lp) %>%
  ungroup()

pData %>%
  bind_rows %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = which), alpha = 0.25) +
  geom_line(data = tpData, alpha = 0.5) +
  facet_grid(temp ~ corr) +
  theme_set(theme_gray(base_size = 18))

pData2 %>%
  bind_rows %>%
  #filter(which %in% c("test1.dat", "test2.dat", "test3.dat")) %>%
  ggplot(aes(mu, value)) +
  geom_point(aes(colour = which), alpha = 0.25) +
  geom_line(data = tpData, alpha = 0.5) +
  facet_grid(temp ~ corr) +
  theme_set(theme_gray(base_size = 18))

tpData %>%
  filter(corr == 1 & temp == 1.0) %>%
  ggplot(aes(mu, value)) +
  geom_line(alpha = 0.5)# +
#  facet_grid(temp ~ corr)

# Compare results of optimization to true ecis
data$opts2 %>%
  select(nonZero, which) %>%
  mutate(lp = lps$lp, opt = lps$opt) %>%
  group_by(which) %>%
  top_n(10, -lp) %>%
  ungroup() %>%
  bind_rows(ecis[nonZero] %>%
              setNames(levels) %>%
              t %>%
              as.tibble %>%
              mutate(which = "TRUTH")) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = levels)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 2, stroke = 2, position = position_dodge(width = 0.75))

#runSimulation(opts2[[i]])
# getResults(path) %>% mutate(data = FALSE) %>%
#   bind_rows(datr) %>%
#   select(data, param_chem_pota, starts_with("corr")) %>%
#   select(1:2, mixedsort(names(.))[14:18]) %>%
#   gather(corr, avg_value, starts_with("corr")) %>%
#   ggplot(aes(param_chem_pota, avg_value)) +
#   geom_point(aes(group = corr, colour = corr, shape = data), alpha = 0.5, size = 2) +
#   theme_grey(base_size = 18)

# Make convex hull plots
hull = data$tclex %>%
  group_by(comp, which) %>%
  filter(row_number() == which.min(formation_energy)) %>%
  ungroup() %>%
  mutate(reference = TRUE)

hullNames = hull %>%
  pull(configname)

data$clexes %>%
  filter(configname %in% hullNames) %>%
  mutate(reference = FALSE) %>%
  bind_rows(hull) %>%
  ggplot(aes(comp, formation_energy)) +
  geom_point(aes(color = reference), shape = 4, size = 2, stroke = 2, alpha = 0.5) +
  facet_wrap( ~ which) +
  theme_set(theme_gray(base_size = 18))

## Make the cooling run plots
data$crs %>%
  #left_join(lps_filtered, by = c("opt", "which" = "which.x")) %>%
  #drop_na() %>%
  mutate(chem = factor(param_chem_pota, levels = sample(unique(param_chem_pota)))) %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_point(aes(group = chem, colour = chem), alpha = 0.5) +
  geom_point(data = data$tcr, aes(group = param_chem_pota), colour = "black", alpha = 0.5) +
  geom_hline(aes(yintercept = 1.0), color = "black") +
  geom_hline(aes(yintercept = 0.5), color = "black") +
  facet_grid(which ~ .) +
  theme_set(theme_gray(base_size = 18))

crs %>%
  bind_rows %>%
  mutate(which = data$crs$which,
         opt = data$crs$opt) %>%
  #left_join(lps_filtered, by = c("opt", "which" = "which.x")) %>%
  #drop_na() %>%
  mutate(chem = factor(param_chem_pota, levels = sample(unique(param_chem_pota)))) %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_point(aes(group = chem, colour = chem), alpha = 0.5) +
  geom_point(data = data$tcr, aes(group = param_chem_pota), colour = "black", alpha = 0.5) +
  facet_grid(which ~ .)

## Recompute coolingRun
crs = list()
for(j in 1:nrow(data$opts2)) {
  crs[[j]] = coolingRun2(path, data$opts2[j,] %>% select(starts_with("corr")) %>% unlist(), seq(0.1, 1.0, length = 20)) %>%
    mutate(which = "optimization", opt = j)
  cat("Finished cr: ", j, " of ", nrow(data$opts2), "\n")
}


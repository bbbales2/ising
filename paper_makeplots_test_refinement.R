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

paths = c("/home/bbales2/ising/paper_outputs/test1.2.dat",
          "/home/bbales2/ising/paper_outputs/test1.refined.dat",
          "/home/bbales2/ising/paper_outputs/test2.2.dat",
          "/home/bbales2/ising/paper_outputs/test2.refined.dat",
          "/home/bbales2/ising/paper_outputs/test3.2.dat",
          "/home/bbales2/ising/paper_outputs/test3.refined.dat")

env = new.env()
load(paths[3], envir = env)

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

data = getData(paths, c('opts', 'opts2', 'clexes', 'tclex', 'crs', 'tcr'))

# Compare results of optimization to true ecis
levels = names(data$opts2)[nonZero]
whichReplace = list("test1.2.dat" = "test1",
                    "test1.refined.dat" = "test1.refined",
                    "test2.2.dat" = "test2",
                    "test2.refined.dat" = "test2.refined",
                    "test3.2.dat" = "test3",
                    "test3.refined.dat" = "test3.refined",
                    "TRUTH" = "TRUTH")
data$opts2 %>%
  #filter(which %in% c("test1.dat", "test2.dat", "test3.dat")) %>%
  select(nonZero, which) %>%
  bind_rows(ecis[nonZero] %>%
              setNames(levels) %>%
              t %>%
              as.tibble %>%
              mutate(which = "TRUTH")) %>%
  mutate(which = unlist(whichReplace[which])) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = levels)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 2, stroke = 2, position = position_dodge(width = 0.75)) +
  theme_set(theme_gray(base_size = 18))

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
  mutate(chem = factor(param_chem_pota, levels = sample(unique(param_chem_pota)))) %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_point(aes(group = chem, colour = chem), alpha = 0.5) +
  geom_point(data = data$tcr, aes(group = param_chem_pota), colour = "black", alpha = 0.5) +
  geom_hline(aes(yintercept = 1.0), color = "black") +
  geom_hline(aes(yintercept = 0.5), color = "black") +
  facet_grid(which ~ .) +
  theme_set(theme_gray(base_size = 18))

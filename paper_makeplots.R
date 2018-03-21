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

# Compare results of optimization to true ecis
do.call(rbind, opts2) %>% as.tibble %>%
  mutate(which = "optimization") %>%
  bind_rows(ecis %>% setNames(names(opts3[[1]]))) %>%
  mutate(which = replace(which, is.na(which), "truth")) %>%
  select(nonZero, which) %>%
  mutate(opt = row_number()) %>%
  gather(corr, eci, starts_with("corr")) %>%
  mutate(corr = factor(corr, levels = corr %>% unique %>% mixedsort)) %>%
  ggplot(aes(corr, eci)) +
  geom_point(aes(color = which), shape = 4, size = 3, stroke = 1)

setTemperatureFraction(path, 0.5)
i = 1

runSimulation(ecis)
datr = getResults(path) %>% mutate(data = TRUE)
runSimulation(opts2[[i]])
getResults(path) %>% mutate(data = FALSE) %>%
  bind_rows(datr) %>%
  select(data, param_chem_pota, starts_with("corr")) %>%
  select(1:2, mixedsort(names(.))[14:18]) %>%
  gather(corr, avg_value, starts_with("corr")) %>%
  ggplot(aes(param_chem_pota, avg_value)) +
  geom_point(aes(group = corr, colour = corr, shape = data), alpha = 0.5, size = 2) +
  theme_grey(base_size = 18)

# Make convex hull plots
tclex = getClex(path, ecis) %>% mutate(which = "truth")
clexes = list()
for(j in 1:length(opts2)) {
  clexes[[j]] = getClex(path, opts2[[j]]) %>% mutate(which = "optimization", opt = j)
}

hull = tclex %>% group_by(comp) %>% filter(row_number() == which.min(formation_energy))

do.call("bind_rows", clexes) %>%
  filter(configname %in% (hull %>% pull(configname))) %>%
  bind_rows(hull) %>%
  ggplot(aes(comp, formation_energy)) +
  geom_point(aes(color = which), shape = 4)

do.call("bind_rows", clexes) %>%
  group_by(comp, opt) %>%
  mutate(is_minimum = (row_number() == which.min(formation_energy))) %>%
  ungroup() %>%
  filter(configname %in% (hull %>% pull(configname))) %>%
  group_by(comp) %>%
  mutate(n = n()) %>%
  mutate(num_mins = sum(is_minimum)) %>%
  mutate(num_not_mins = n()) %>%
  ungroup() %>%
  gather(which_count, num, num_mins, num_not_mins) %>%
  ggplot(aes(comp, num)) +
  geom_bar(aes(fill = which_count), stat = "identity")

## Make the cooling run plots

crs %>% bind_rows %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_point(aes(colour = param_chem_pota), alpha = 0.1) +
  geom_path(data = tcr, aes(group = param_chem_pota), colour = "red")

crs %>% bind_rows %>%
  mutate(chem = factor(param_chem_pota, levels = sample(unique(param_chem_pota)))) %>%
  ggplot(aes(corr1, Tfrac)) +
  geom_path(aes(colour = chem), alpha = 0.5) +
  geom_path(data = tcr, aes(group = param_chem_pota), colour = "black", alpha = 0.5)

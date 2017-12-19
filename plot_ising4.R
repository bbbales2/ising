library(reshape2)
library(tidyverse)
library(ggplot2)
library(GGally)

opt_df = opts %>% bind_rows %>%
  group_by(which_opt) %>%
  mutate(type = "opt") %>%
  ungroup %>%
  bind_rows(beta %>% as.list %>% as.tibble %>% mutate %>%
              mutate(which_opt = 0, type = "truth", lp = 0, rn = 1)) %>%
  bind_rows(rnorm(5 * 2000, 0.1, 0.25) %>% matrix(ncol = 5) %>% as.tibble %>%
              setNames(names(b)) %>% mutate(which_opt = 0, type = "prior", lp = 0, rn = 1))

opt_plot = pmap(combn(names(b), 2) %>% t %>% as.tibble,
                ~ opt_df %>% select(..1, ..2, which_opt, lp, rn, type) %>%
                  rename(x = !!..1, y = !!..2) %>%
                  mutate(which_x = !!..1,
                         which_y = !!..2)) %>%
  bind_rows

opt_dens = map(names(b) %>% as.list, ~ opt_df %>%
                 rename(x = !!.x) %>%
                 mutate(which = !!.x)) %>%
  bind_rows

opt_dens %>% filter(lp > -5) %>%
  ggplot(aes(x)) +
  geom_density(aes(colour = type, fill = type), alpha = 0.15) +
  facet_grid(. ~ which)

opt_plot %>% filter(type == "opt" & lp > -5 & y < 1.0 & y > -1.0 & x < 1.0 & x > -1.0) %>%
  ggplot(aes(y, x)) +
  geom_density2d(data = opt_plot %>% filter(type == "prior"), bins = 10) +
  geom_point(aes(group = which_opt, colour = lp), size = 0.5) +
  scale_colour_gradient(low = "blue", high = "red") +
  geom_point(data = opt_plot %>% filter(type == "truth"), colour = "black", size = 2.0, shape = 4, stroke = 2) +
  facet_grid(which_x ~ which_y)

opt_samples %>% filter(x == 0.0 & lp > -5.0) %>% gather(coeff, value, starts_with("Q")) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  geom_vline(data = data$Q %>% as.list %>% as.tibble %>% gather(coeff, value),
             aes(xintercept = value), color = "red") +
  facet_grid(. ~ coeff, scales = "free_x")

opt_samples %>% filter(lp > -5.0) %>%
  ggplot(aes(x, X0)) +
  geom_line(aes(group = opt), alpha = 0.25) +
  geom_point(data = list(x = mus, X0 = data$X0, which = "data") %>% as.tibble, color = "red") +
  xlab("mu") + ylab("avg. composition") +
  ggtitle("truth are red dots,\ndistribution of data generated from fits in black/grey\nmean + 95% interval estimates")


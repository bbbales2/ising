library(jsonlite)
library(tidyverse)
library(stringr)

runMC = function(path) {
  system(paste0("(cd ",
                path,
                "; rm -rf conditions.* results.json; /home/bbales2/local/casm/bin/casm monte -s monte.json)"))
}

getResults = function(path) {
  df = fromJSON(paste0(path, "/results.json")) %>% as.tibble
  colnames(df) = gsub('[<>\\(\\),]', '', colnames(df))
  df
}

getCorrs = function(path) {
  files = list.files(path, pattern = "conditions.")
  fileNum = as.numeric(gsub('^conditions.([0123456789]*)$', '\\1', files))
  files = files[order(fileNum)]

  dfs = list()  
  for(i in 1:length(files)) {
    df = fromJSON(paste0(path, "/", files[[i]], "/observations.json.gz")) %>% as.tibble %>%
      select(starts_with("corr"))
    
    colnames(df) = gsub('[\\(\\)]', '', colnames(df))
    colNumbers = order(as.numeric(gsub('^corr([0123456789]*)$', '\\1', colnames(df))))
    dfs[[i]] = df %>% select(noquote(colnames(df)[colNumbers])) %>%
      mutate(t = i, mci = row_number())
  }
  
  dfs %>% bind_rows %>% select(t, everything())
}

getECIs = function(path) {
  fpath = paste0(path, "/cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.default/eci.json")
  a = fromJSON(fpath)
  out = a$cluster_functions$eci
  if(is.null(out)) {
    out = rep(0, nrow(a$cluster_functions))
  }
  out
}

#getECIs("/home/bbales2/casm/cubic_3d/")

setECIs = function(path, ecis) {
  fpath = paste0(path, "/cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.default/eci.json")
  a = fromJSON(fpath)
  a$cluster_functions$eci = ecis
  writeLines(toJSON(a, pretty = TRUE), fpath)
}

#setECIs("/home/bbales2/casm/cubic_3d/", rep(0, 61))



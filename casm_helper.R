library(jsonlite)
library(tidyverse)
library(stringr)
library(parallel)

runMC = function(path) {
  a = fromJSON(paste0(path, "/monte.json"))
  
  iv = a$driver$initial_conditions$param_chem_pot$a
  fv = a$driver$final_conditions$param_chem_pot$a
  sv = a$driver$incremental_conditions$param_chem_pot$a
  
  system(paste0("(cd ",
                path,
                "; rm -rf sim.*)"))
  
  vs = seq(iv, fv, sv)
  
  for(i in 1:length(vs)) {
    system(paste0("(mkdir ", path, "/sim.", i, ")"))
    a$driver$initial_conditions$param_chem_pot$a = vs[i]
    a$driver$final_conditions$param_chem_pot$a = vs[i]
    a$driver$incremental_conditions$param_chem_pot$a = 0.0
    writeLines(toJSON(a, pretty = TRUE, auto_unbox = TRUE), paste0(path, "/sim.", i, "/monte.json"))
  }
  
  out = mclapply(1:length(vs), function(i) {
    system(paste0("(cd ",
                  paste0(path, "/sim.", i),
                  "; /home/bbales2/local/casm/bin/casm monte -s monte.json)"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
  }, mc.cores = 20)

  #system(paste0("(cd ",
  #              path,
  #              "; rm -rf conditions.* results.json; /home/bbales2/local/casm/bin/casm monte -s monte.json)"),
  #       ignore.stdout = TRUE, ignore.stderr = TRUE)
}

getResults = function(path) {
  files = list.files(path, pattern = "sim.")
  fileNum = as.numeric(gsub('^sim.([0123456789]*)$', '\\1', files))
  files = files[order(fileNum)]
  
  dfs = list()
  for(i in 1:length(files)) {
    df = fromJSON(paste0(path, "/", files[i], "/results.json")) %>% as.tibble
    colnames(df) = gsub('[<>\\(\\),]', '', colnames(df))
    dfs[[i]] = df
  }
  
  do.call("bind_rows", dfs)
}

getCorrs = function(path) {
  files = list.files(path, pattern = "sim.")
  fileNum = as.numeric(gsub('^sim.([0123456789]*)$', '\\1', files))
  files = files[order(fileNum)]

  dfs = list()  
  for(i in 1:length(files)) {
    df = fromJSON(paste0(path, "/", files[[i]], "/conditions.0/observations.json.gz")) %>% as.tibble %>%
      select(starts_with("corr"))
    
    colnames(df) = gsub('[\\(\\)]', '', colnames(df))
    colNumbers = order(as.numeric(gsub('^corr([0123456789]*)$', '\\1', colnames(df))))
    dfs[[i]] = df %>% select(noquote(colnames(df)[colNumbers])) %>%
      mutate(mu = i, mci = row_number())
  }
  
  dfs %>% bind_rows %>% select(mu, everything())
}

getECIs = function(path) {
  fpath = paste0(path, "/cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.default/eci.json")
  a = fromJSON(fpath)
  out = a$cluster_functions$eci
  if(is.null(out)) {
    out = rep(0, nrow(a$cluster_functions))
  }
  out[is.na(out)] = 0.0
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

setSupercell = function(path, N) {
  a = fromJSON(paste0(path, "/monte.json"))
  a$supercell[1, 1] = N
  a$supercell[2, 2] = N
  writeLines(toJSON(a, pretty = TRUE, auto_unbox = TRUE), paste0(path, "/monte.json"))
}

getSupercell = function(path) {
  a = fromJSON(paste0(path, "/monte.json"))
  a$supercell[1, 1]
}

getClex = function(path, b = NULL) {
  if(!is.null(b)) {
    setECIs(path, b)
  }
  
  system(paste0("(cd ", path,
                "; /home/bbales2/local/casm/bin/casm query -k 'clex(formation_energy) comp' -j -o clex_query.json)"),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  fromJSON(paste0(path, "/clex_query.json")) %>% as.tibble %>% unnest %>% rename(formation_energy = 'clex(formation_energy)')
}

getChemicalPotentials = function(path) {
  a = fromJSON(paste0(path, "/monte.json"))
  seq(a$driver$initial_conditions$param_chem_pot$a,
      a$driver$final_conditions$param_chem_pot$a,
      a$driver$incremental_conditions$param_chem_pot$a)
}

setTemperatureFraction = function(path, frac) {
  a = fromJSON(paste0(path, "/monte.json"))
  a$driver$initial_conditions$temperature = 11604.97 * frac
  a$driver$final_conditions$temperature = 11604.97 * frac
  a$driver$incremental_conditions$temperature = 0.0
  writeLines(toJSON(a, pretty = TRUE, auto_unbox = TRUE), paste0(path, "/monte.json"))
}

getTemperatureFraction = function(path, frac) {
  a = fromJSON(paste0(path, "/monte.json"))
  a$driver$initial_conditions$temperature / 11604.97
}

coolingRun = function(path, b, fracs, debug = FALSE) {
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
  
  out = map(fracs, function(frac) {
    if(debug) {
      cat("Running simulations for frac: ", frac, "\n");
    }
    
    setTemperatureFraction(path, frac)
    dat = runSimulation(b)
    
    list(corr1 = map(dat, ~ .$Eg[[1]]) %>%
           unlist()) %>%
      as.tibble %>%
      mutate(param_chem_pota = getChemicalPotentials(path),
             Tfrac = frac)
  }) %>% bind_rows
  
  setTemperatureFraction(path, 1.0)
  
  out
}

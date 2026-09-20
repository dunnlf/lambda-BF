### Applying lambda-BF to example data from papers
###
### Data from the following are considered:
###    - Gônet, J. et al. (2023)
###    - Round, E., Dockum, R., Ryder, R. J. (2022)
###    - Arato, J., Fitch, W. T. (2021)
###

library(ape)
library(phytools)
library(lambdaBF)
library(geiger)

# ---- Arato & Fitch ----

# fixing species mismatches due to synonyms
# replace species in data by those found in Ericson birdtree set
data_missing_species <- c('Ceryle_alcyon', 'Chroicocephalus_ridibundus',
                          'Mareca_americana', 'Cecropis_daurica',
                          'Poecile_atricapillus', 'Poecile_palustris')
synonyms <- c('Megaceryle_alcyon', 'Larus_ridibundus', 'Anas_americana',
              'Hirundo_daurica', 'Parus_atricapillus', 'Parus_palustris')

synonym_df <- data.frame(missing_species=data_missing_species, synonym=synonyms)

replace_synonyms <- function(df, synonym_df){
  df_ <- df
  replaced <- 0
  for (i in 1:nrow(df_)){
    species_i <- df_[i,'Species']
    if (species_i %in% synonym_df[,'missing_species']) {
      df_[i, 'Species'] <- synonym_df[synonym_df['missing_species']==species_i, 'synonym']
      replaced <- replaced + 1
    }
  }
  print(paste('Replaced ', replaced, ' species names'))
  return(df_)
}

# acoustic data 
nonpasserine_feat <- read.csv("lambda-BF-example-data/arato-fitch/All Non-Passerines.csv", row.names=1)
passerine_song_feat <- read.csv("lambda-BF-example-data/arato-fitch/All Passerine Song.csv", row.names=1)

nonpasserine_tree <- read.nexus("lambda-BF-example-data/arato-fitch/All_NonPasserines.nex")
passerine_song_tree <- read.nexus("lambda-BF-example-data/arato-fitch/All_Passerine_Song.nex")

passerine_song_feat <- replace_synonyms(passerine_song_feat, synonym_df)
nonpasserine_feat <- replace_synonyms(nonpasserine_feat, synonym_df)

# name checks
check_nonpass <- name.check(nonpasserine_tree[[1]], data.names=nonpasserine_feat[,'Species'])
check_pass <- name.check(passerine_song_tree[[1]], data.names=passerine_song_feat[,'Species'])

pass_species <- passerine_song_feat[,'Species']
nonpass_species <- nonpasserine_feat[,'Species']

# compute BFs for each trait
traits <- c('RMS', 'SpectCentroid', 'SpectralFlux', 'Entropy',
            'Flatness', 'Contrast1', 'Contrast2', 'Contrast3', 'Contrast4')
res_pass <- list()
for (i in seq_along(traits)){
    t <- traits[i]
    trait <- setNames(passerine_song_feat[,t], pass_species)
    res_i <- lambdaBF(passerine_song_tree, trait,a=3, d=5)
    res_pass[[i]] <- res_i
}

res_nonpass <- list()
for (i in seq_along(traits)){
    t <- traits[i]
    trait <- setNames(nonpasserine_feat[,t], nonpass_species)
    res_i <- lambdaBF(nonpasserine_tree, trait,a=3, d=5)
    res_nonpass[[i]] <- res_i
}


# print summaries
bfs_pass <- sapply(res_pass, function(x) x$bf)
bfs_nonpass <- sapply(res_nonpass, function(x) x$bf)

print("Bayes Factors for Passerine Song:")
print(log10(bfs_pass))

print("Bayes Factors for Non-Passerine:")
print(log10(bfs_nonpass))

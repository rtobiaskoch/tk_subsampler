rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_unload()

pkg = c("tidyverse", "seqinr")

pacman::p_load(pkg, character.only = T)

#list all metadata files
fn = list.files(pattern = "metadata.tsv", 
                recursive = T, 
                full.names = T)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#METADATA ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#read files
t = map(fn, ~read_tsv(.,
                      col_types = rep("c", ncol(read_delim(., delim = "\t", n_max = 1))))
        ) 

#data <- read_tsv("myfile.tsv", col_types = rep("c", ncol(read_delim("myfile.tsv", delim = "\t", n_max = 1))))

#create function to merge
to_character <- function(df) {
  df %>% mutate_all(as.character)
}

# apply the function to each data frame in the list
mdata <- map(t, to_character) %>%
  bind_rows()

save(mdata, file = "data_mid/metadata.RData")



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FASTA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#list all metadata files
fn = list.files(pattern = ".fasta$", 
                recursive = T, 
                full.names = T)


#read files
t = map(fn, read.fasta) 

fasta = do.call(c, t)
  
# apply the function to each data frame in the list

save(fasta, file = "data_mid/fasta.RData")





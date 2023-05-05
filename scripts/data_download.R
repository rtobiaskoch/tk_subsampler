rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_unload()

pkg = c("tidyverse", "seqinr")

pacman::p_load(pkg, character.only = T)

#list all metadata files
fn = list.files("data_input",
                pattern = "metadata.tsv", 
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
  bind_rows() %>%
  distinct_all()

save(mdata, file = "data_mid/metadata.RData")

rm(t)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FASTA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#list all metadata files
fn = list.files("data_input",
                pattern = ".fasta$", 
                recursive = T, 
                full.names = T)

sero = str_extract(fn, "(?<=/)[^/]+")




#read files
fasta = map(fn, read.fasta) 
names(fasta) = sero

fasta_all = do.call(c, fasta)

#strains in fasta but not in metadata
x = list(setdiff(names(fasta_all),mdata$strain))
x = data.frame(x)
write.csv(x, "data_output/missing sample from metadata.csv")
  
# apply the function to each data frame in the list

save(fasta, file = "data_mid/fasta.RData")





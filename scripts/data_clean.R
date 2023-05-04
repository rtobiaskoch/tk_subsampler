rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_unload()

pkg = c("tidyverse", "lubridate", "janitor", "seqinr")

pacman::p_load(pkg, character.only = T)

set.seed(420)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# METADATA ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>
load("data_mid/metadata.RData")

caribbean_countries <- c("Anguilla|Antigua-and-Barbuda|Aruba|Bahamas|Barbados|Bonaire, Sint-Eustatius-and-Saba|
                         British-Virgin-Islands|Cayman-Islands|Cuba|Curaçao|Dominica|Dominican-Republic|
                         Grenada|Guadeloupe|Haiti|Jamaica|Martinique|Montserrat|Puerto-Rico|
                         Saint-Barthelemy|Saint-Kitts-and-Nevis|Saint-Lucia|Saint-Martin|Saint-Vincent-and-the-Grenadines|
                         Sint-Maarten|Trinidad-and-Tobago|Turks-and-Caicos-Islands|United-States-Virgin-Islands")


mdata_c = mdata %>%
  mutate(yale = if_else(grepl("^Yale", strain),
                        1,0),
        caribbean = if_else(grepl(caribbean_countries, Location),
                              1,0),
         florida = if_else(grepl("^Florida", Location),
                           1,0)
        ) %>%
  mutate(year = substr(date, 1,4)) %>%
  group_by(Location, year) %>% 
  mutate(n_obs = n()) %>% #get size of group to filter for downsample
  ungroup()
  

t = distinct(mdata_c, Location, caribbean, florida)


keep = mdata_c %>%
  filter(florida == 0 &
         caribbean == 1 | yale == 1)

keep_lt10= mdata_c %>%
  filter(florida == 0 &
         caribbean == 0 &
         yale == 0,
         n_obs <= 10) 

ds = mdata_c %>%
  filter(florida == 0 &
         caribbean == 0 &
         yale == 0,
         n_obs > 10)

#before sampling
t = ds %>% group_by(Location, year) %>%
  dplyr::count()

hist(t$n, breaks = seq(0,500,10))

#downsample 
ds10 = ds %>% 
  group_by(Location, year) %>% 
  sample_n(10) %>% 
  ungroup()

mdata_ds = rbind(keep, keep_lt10, ds10) %>%
  group_by(Location, year) %>% 
  mutate(n_obs_ds = n()) %>% #get size of group to filter for downsample
  ungroup()

t = mdata_ds %>% 
  group_by(Location, year, yale, caribbean) %>%
  dplyr::count() %>%
  mutate(yale = if_else(yale == 1, "yale", "other"),
         caribbean = as.factor(caribbean)) %>%
  arrange(desc(caribbean),Location, year, desc(n))

ggplot(t, aes(x = year, y = n, color = caribbean, size = n)) +
     geom_point(alpha = 0.5) +
     theme(axis.text.x = element_text(angle = 90)) +
     theme_classic() +
     scale_x_discrete(breaks = seq(1900,2020,10)) +
     scale_y_continuous(breaks = seq(0,200,10)) +
     facet_wrap(~yale, ncol = 1)

write_tsv(mdata_ds, "data_output/metadata_ds.tsv")

#END METADATA ####

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FASTA ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

filter_id = mdata_ds$strain

devtools::source_url("https://raw.githubusercontent.com/rtobiaskoch/TK_useful_functions/main/fasta_filter.R")

load("data_mid/fasta.RData")

fasta_ds = filter_fasta(fasta, filter_id)

#find missing
f_nm = names(fasta_ds)

miss_id = setdiff(filter_id, f_nm)

write.fasta(sequences = fasta_ds,
            names = names(fasta_ds),
            file.out = "data_output/fasta_ds.fasta")

#github test

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


#create subsetting variables for easier filtering
mdata_c = mdata %>%
  mutate(yale = if_else(grepl("^Yale", strain),
                        1,0),
        caribbean = if_else(grepl(caribbean_countries, Location),
                              1,0),
        florida = if_else(grepl("^Florida", Location),
                           1,0),
        cmplt_date = if_else(grepl("XX$", date),
                                0,1)
        ) %>%
  mutate(year = substr(date, 1,4)) %>%
  group_by(Location, year) %>% 
  mutate(n_obs = n()) %>% #get size of group to filter for downsample
  group_by(Location, year, cmplt_date) %>%
  mutate(n_obs_date = n()) %>%
  ungroup()

t = mdata_c %>%
  get_dupes(Location, year)

write.csv(mdata_c, "data_mid/mdata_c.csv")

t = distinct(mdata_c, Location, caribbean, florida)

#variables to keep in final dataset
keep = mdata_c %>%
  filter(florida == 0) %>%
  filter(caribbean == 1 | yale == 1)

#variables that aren't florida, caribbean or yale but have less than 10 samples so dont need to be downsampled
keep_lt10 =  mdata_c %>%
  filter(florida == 0) %>%
  filter(caribbean == 0 &
         yale == 0 &
         n_obs <= 10)
  
#variables that aren't florida, caribbean or yale but have more than 10 samples so dont need to be downsampled
ds = mdata_c %>%
  filter(florida == 0) %>%
  filter(caribbean == 0 &
         yale == 0 &
         n_obs > 10)


#before sampling
t = ds %>% group_by(Location, year) %>%
  dplyr::count()

hist(t$n, breaks = seq(0,500,10))

#downsample 
ds10 = ds %>%
  sample_n(nrow(.)) %>% #scramble to create randomization
  group_by(Location, year) %>%
  arrange(desc(cmplt_date)) %>% #prioritize complete date samples
  slice_head(n = 10) %>% #take top 10
  ungroup()

mdata_ds = rbind(keep, keep_lt10, ds10) %>%
  group_by(Location, year) %>% 
  mutate(n_obs_ds = n()) %>% #get size of group to filter for downsample
  ungroup()

#check prioritization of complete dates worked
t = mdata_ds %>%
  group_by(Location, year, cmplt_date, n_obs, n_obs_date) %>%
  dplyr::count() %>%
  ungroup() %>%
  get_dupes(Location, year) %>%
  filter(cmplt_date == 1)


#for plot
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


#split up data again 

sero = unique(mdata_ds$Serotype)

fn = paste0("data_output/", sero, "_metadata.tsv")

df_list = mdata_ds %>%
  group_split(Serotype)

map2(.x = df_list, .y = fn, 
     ~write_tsv(.x, .y))


write_tsv(mdata_ds, "data_output/all_metadata_ds.tsv")

#END METADATA ####

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FASTA ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

filter_id = mdata_ds$strain

devtools::source_url("https://raw.githubusercontent.com/rtobiaskoch/TK_useful_functions/main/fasta_filter.R")

load("data_mid/fasta.RData")

fasta_all = do.call(c, fasta)

fasta_ds = map(fasta, ~filter_fasta(.x, filter_id))

fasta_ds_all = do.call(c, fasta_ds)
length(unique(names(fasta_ds_all))) == length(names(fasta_ds_all)) #should be TRUE



#iteratively write fasta files
sero2 = names(fasta_ds)

fn2 = paste0("data_output/", sero2, "_sequences.fasta")

map2(.x = fasta_ds, .y = fn2, 
     ~write.fasta(sequences = .x, 
                  names = names(.x),
                  file.out = .y)
    )

#if you are in R studio, run the next line
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### Initialize ####

library(tidyverse)
library(magrittr)
library(janitor)

#### User Defined Variables ####
alleledepthFILE = '../leiognathus_leuciscus/filterVCF/lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.100.AD.tsv'
genoFILE = '../leiognathus_leuciscus/filterVCF/lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.100.GT.tsv'
alleledepthFILE = '../leiognathus_leuciscus/filterVCF/lle.B.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr15.5.recode.AD.tsv'
genoFILE = '../leiognathus_leuciscus/filterVCF/lle.B.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr15.5.recode.GT.tsv'
indvPATTERN = "lle_"
modIdPATTERN = "_ll_.*"
het_cutoff_pop1 = 0.1
het_cutoff_pop3 = 0.125

alleledepthFILE = '../atherinomorus_endrachtensis/filterVCF_ceb/Aen.A3.rad.RAW-6-6.Fltr07.19.AD.tsv'
genoFILE = '../atherinomorus_endrachtensis/filterVCF_ceb/Aen.A3.rad.RAW-6-6.Fltr07.19.GT.tsv'
alleledepthFILE = '../atherinomorus_endrachtensis/filterVCF_ceb/Aen.A3.rad.RAW-6-6.Fltr15.9.recode.AD.tsv'
genoFILE = '../atherinomorus_endrachtensis/filterVCF_ceb/Aen.A3.rad.RAW-6-6.Fltr15.9.recode.GT.tsv'
indvPATTERN = "aen_"
modIdPATTERN = "_ae_.*"
het_cutoff_pop1 = 0.2125
het_cutoff_pop3 = 0.275

#### Calculate Variables ####
outDIR=dirname(genoFILE)
outFilePREFIX <-
  basename(genoFILE) %>%
  str_remove("GT.tsv")

#### Define Functions ####

getMode <- function(x) {
  keys <- 0:100 / 100
  keys[which.max(tabulate(match(round(x,
                                      2), 
                                keys)))]
}

# ex usage
# all_data %>%
# pull(allele_balance) %>%
# getMode()

  

#### READ in DATA ####
allele_depths <-
  read_tsv(alleledepthFILE,
           col_types=cols(.default = "c")) %>%
  clean_names() %>%
  pivot_longer(cols=contains(indvPATTERN),
               names_to="id") %>%
  separate(col=value,
           into=(c("num_reads_ref",
                  "num_reads_alt")),
           convert=TRUE) %>%
  mutate(num_reads = num_reads_alt + num_reads_ref)

# geno_likes <-
#   read_tsv(genolikelihoodFILE,
#            col_types=cols(.default = "c")) %>%
#   clean_names() %>%
#   pivot_longer(cols=contains(indvPATTERN),
#                names_to="id") %>%
#   separate(col=value,
#            into=(c("like_homo_ref",
#                    "like_hetero",
#                    "like_homo_alt")),
#            convert=TRUE)

genotypes <-
  read_tsv(genoFILE,
           col_types=cols(.default = "c")) %>%
  clean_names() %>%
  pivot_longer(cols=contains(indvPATTERN),
               names_to="id") %>%
  # mutate(value = str_replace_all(value,
  #                            pattern = "\\.",
  #                            replacement = "NA")) %>%
  separate(col=value,
           into=(c("state_1",
                   "state_2"))) %>%
  mutate(state_1 = na_if(state_1,
                        ""),
         state_2 = na_if(state_2,
                         ""),
         genotype = case_when(state_1 != state_2 ~ "hetero",
                              is.na(state_1) & is.na(state_2) ~ NA_character_,
                              state_1 == ref ~ "homo_ref",
                              state_1 == alt ~ "homo_alt",
                              TRUE ~ "error"))

#### COMBINE DATA ####
num_loci <-
  genotypes %>%
  mutate(chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) %>%
  pull(chrom_pos) %>%
  unique() %>%
  length()

# heterozygosity by individual
genotypes %>%
  group_by(id,
           genotype) %>%
  summarize(num_pos = n()) %>%
  pivot_wider(names_from = "genotype",
              values_from = "num_pos") %>%
  mutate(heterozygosity_obs_ind = hetero/(hetero + homo_ref + homo_alt),
         pop = case_when(heterozygosity_obs_ind >= het_cutoff_pop3 ~ 3,
                         heterozygosity_obs_ind >= het_cutoff_pop1 ~ 1,
                         TRUE ~ 2)) %>%
  select(id,
         heterozygosity_obs_ind,
         pop) %>%
  ggplot(aes(x=id,
             y=heterozygosity_obs_ind,
             color=pop)) +
  geom_point()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# heterozygosity by site (There's a problem here, think I fixed it)
# genotypes %>%
#   mutate(chrom_pos = str_c(chrom,
#                            pos,
#                            sep="_")) %>%
#   separate(id,
#            into=(c(NA,
#                    "era",
#                    NA,
#                    NA,
#                    NA,
#                    NA)),
#            remove = FALSE) %>%
#   group_by(chrom,
#            pos,
#            chrom_pos,
#            genotype,
#            era) %>%
#   summarize(num_pos = n()) %>%
#   pivot_wider(names_from = "genotype",
#               values_from = "num_pos") %>% 
#   mutate(across(hetero:homo_alt,
#                 ~replace_na(.,
#                             0)),
#          heterozygosity_obs_locus = hetero/(hetero + homo_ref + homo_alt)) %>%
#   select(chrom,
#          pos,
#          chrom_pos,
#          era,
#          heterozygosity_obs_locus) %>%
#   ggplot(aes(x=chrom_pos,
#              y=heterozygosity_obs_locus,
#              color=era)) +
#   geom_point()

# read in data rather than run the next line of code if rds exists  
# all_data <-
#   readRDS(file = str_c(outDIR,
#                     "/",
#                     outFilePREFIX,
#                     "alldata.rds",
#                     sep=""))

all_data <-
  allele_depths %>%
  # left_join(geno_likes) %>%
  left_join(genotypes) %>%
  left_join(genotypes %>%
              group_by(id,
                       genotype) %>%
              summarize(num_pos = n()) %>%
              pivot_wider(names_from = "genotype",
                          values_from = "num_pos") %>%
              mutate(across(hetero:homo_alt,
                            ~replace_na(.,
                                        0)),
                     heterozygosity_obs_ind = hetero/(num_loci - `NA`),
                     pop = case_when(heterozygosity_obs_ind >= het_cutoff_pop3 ~ 3,
                                     heterozygosity_obs_ind >= het_cutoff_pop1 ~ 1,
                                     TRUE ~ 2)) %>%
              select(id,
                     heterozygosity_obs_ind,
                     pop)) %>%
  separate(id,
           into=(c(NA,
                   "era",
                   NA,
                   NA,
                   NA,
                   NA)),
           remove = FALSE) %>%
  left_join(genotypes %>%
              mutate(chrom_pos = str_c(chrom,
                                       pos,
                                       sep="_")) %>%
              separate(id,
                       into=(c(NA,
                               "era",
                               NA,
                               NA,
                               NA,
                               NA)),
                       remove = FALSE) %>%
              group_by(chrom,
                       pos,
                       chrom_pos,
                       genotype,
                       era) %>%
              summarize(num_pos = n()) %>%
              pivot_wider(names_from = "genotype",
                          values_from = "num_pos") %>% 
              mutate(across(hetero:homo_alt,
                            ~replace_na(.,
                                        0)),
                     heterozygosity_obs_locus = hetero/(hetero + homo_ref + homo_alt)) %>% 
              select(chrom,
                     pos,
                     chrom_pos,
                     era,
                     heterozygosity_obs_locus)) %>%
  mutate(qual = as.double(qual),
         across(contains(c("num_reads",
                           "like",
                           "pos")),
                ~as.numeric(.)),
         num_reads = num_reads_alt + num_reads_ref,
         # genotype = case_when(state_1 != state_2 ~ "hetero",
         #                          is.na(state_1) & is.na(state_2) ~ NA_character_,
         #                          state_1 == ref ~ "homo_ref",
         #                          state_1 == alt ~ "homo_alt",
         #                          TRUE ~ "error"),
         allele_balance = case_when(genotype == "hetero" ~ num_reads_alt / num_reads,
                                    TRUE ~ NA_real_),
         id = str_remove(id,
                         modIdPATTERN),
         chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) 

# saveRDS(all_data, 
#         file=str_c(outDIR,
#                    "/",
#                    outFilePREFIX,
#                    "alldata.rds",
#                    sep=""))

#### SUMMARIZE DATA BY CHROM POS ERA POP ####

hetero_data_chrom_pos_era <-
  all_data %>%
  mutate(chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) %>%
  filter(genotype == "hetero") %>% 
  group_by(chrom, 
           pos,
           chrom_pos,
           # pop,
           era) %>%
  summarize(mode_allele_balance = getMode(allele_balance),
            median_allele_balance = median(allele_balance),
            mean_allele_balance = mean(allele_balance),
            mean_heterozygosity_obs_locus = mean(heterozygosity_obs_locus),
            n = n())
  

#### AB By Era ####

# all_data %>%
#   ggplot(aes(x=like_hetero)) +
#   geom_histogram() +
#   facet_grid(genotype ~ era)

all_data %>%
  filter(genotype == "hetero") %>%
  ggplot(aes(x=allele_balance,
             fill = era)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept = c(1/8,
                            1/6,
                            2/8,
                            2/6,
                            3/8,
                            4/8,
                            5/8,
                            4/6,
                            6/8,
                            5/6,
                            7/8),
             color="grey",
             linetype="dashed") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "Histograms of Allele Balance",
       subtitle = "Individuals x Position, Most Granular") +
  facet_grid(era ~ .,
             scales = "free")

ggsave(paste(outDIR, 
             outFilePREFIX,
             'HIST-AB-INDxPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

hetero_data_chrom_pos_era %>% 
  ggplot(aes(x=median_allele_balance,
             fill = era)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept = c(1/8,
                            1/6,
                            2/8,
                            2/6,
                            3/8,
                            4/8,
                            5/8,
                            4/6,
                            6/8,
                            5/6,
                            7/8),
             color="grey",
             linetype="dashed") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "Histograms of Allele Balance",
       subtitle = "Medians by Position") +
  facet_grid(era ~ .,
             scales = "free")

ggsave(paste(outDIR, 
             "/",
             outFilePREFIX,
             'HIST-AB-medPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

hetero_data_chrom_pos_era %>%
  ggplot(aes(x=mode_allele_balance,
             fill = era)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept = c(1/8,
                            1/6,
                            2/8,
                            2/6,
                            3/8,
                            4/8,
                            5/8,
                            4/6,
                            6/8,
                            5/6,
                            7/8),
             color="grey",
             linetype="dashed") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "Histograms of Allele Balance",
       subtitle = "Modes by Position") +
  facet_grid(era ~ .,
             scales = "free")

ggsave(paste(outDIR, 
             "/",
             outFilePREFIX,
             'HIST-AB-modePOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

#### AB By Num Hets ####

hetero_data_chrom_pos_era %>%
  # filter(pop != 3) %>%
  ggplot(aes(x=mean_allele_balance,
             y = n,
             color = era)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "n heterozygotes") +
  facet_grid(era ~ .,
             scales="free_y")

hetero_data_chrom_pos_era %>%
  ggplot(aes(x=mode_allele_balance,
             y = n,
             color = era)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "n heterozygotes") +
  facet_grid(era ~ .,
             scales="free_y")

hetero_data_chrom_pos_era %>%
  # filter(pop != 3) %>%
  ggplot(aes(x=mean_allele_balance,
             y = mean_heterozygosity_obs_locus,
             color = era)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "prop heterozygotes") +
  facet_grid(era ~ .,
             scales="free_y")

hetero_data_chrom_pos_era %>%
  filter(mean_heterozygosity_obs_locus < 0.6,
         mean_allele_balance > 0.375 & mean_allele_balance < 0.625) %>%
  ggplot(aes(x=mean_allele_balance,
             y = mean_heterozygosity_obs_locus,
             color = era)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "prop heterozygotes") +
  facet_grid(era ~ .,
             scales="free_y")


hetero_data_chrom_pos_era %>%
  filter(era != "a",
         mean_heterozygosity_obs_locus < 0.65,
         mean_allele_balance > 0.375 & mean_allele_balance < 0.625) %>%
  ungroup() %>%
  dplyr::select(chrom,
         pos) %>%
  distinct() %>%
  write_tsv(str_c(outDIR, 
                  "/",
                  outFilePREFIX,
                  'chrom_pos_include.tsv', 
                  sep = ""),
            col_names = FALSE)

hetero_data_chrom_pos_era %>%
  filter(era != "a") %>%
  select(chrom,
         pos) %>%
  distinct()

all_data %>%
  mutate(chrom_pos = str_c(chrom,
                           pos,
                           sep="_")) %>%
  filter(genotype != "hetero") %>% 
  # group_by(chrom, 
  #          pos,
  #          chrom_pos,
  #          # pop,
  #          era) %>%
  # summarize(mode_allele_balance = getMode(allele_balance),
  #           median_allele_balance = median(allele_balance),
  #           mean_allele_balance = mean(allele_balance),
  #           mean_heterozygosity_obs_locus = mean(heterozygosity_obs_locus),
  #           n = n())
  filter(era != "a") %>%
  select(chrom,
         pos) %>%
  distinct() %>%
  view()

#### INDIV AB ####

all_data %>%
  filter(genotype == "hetero",
         era == "a") %>%
  mutate(id = str_remove(id,
                         modIdPATTERN)) %>%
  ggplot(aes(x=allele_balance,
             fill = era)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept = c(1/8),
             color="black",
             linetype="solid") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "Histograms of Allele Balance",
       subtitle = "Individuals x Position, Most Granular, Vert Line = 1/8") +
  facet_wrap(era ~ id,
             scales = "free")

ggsave(paste(outDIR, 
             outFilePREFIX,
             'HIST-AB-INDxPOS-ALB.png', 
             sep = ""), 
       height = 9, 
       width = 6.5)


all_data %>%
  filter(genotype == "hetero",
         era != "a") %>%
  mutate(id = str_remove(id,
                         modIdPATTERN)) %>%
  ggplot(aes(x=allele_balance,
             fill = era)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept = c(1/8),
             color="black",
             linetype="solid") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "Histograms of Allele Balance",
       subtitle = "Individuals x Position, Most Granular, Vert Line = 1/8") +
  facet_wrap(era ~ id,
             scales = "free")

ggsave(paste(outDIR, 
             outFilePREFIX,
             'HIST-AB-INDxPOS-CONT.png', 
             sep = ""), 
       height = 9, 
       width = 6.5)

# all_data %>%
#   group_by(id, era) %>%
#   summarize(median_allele_balance = median(allele_balance,
#                                            na.rm = TRUE)) %>%
#   ggplot(aes(x=median_allele_balance,
#              fill = era)) +
#   geom_histogram() +
#   scale_x_continuous(limits = c(0, 1)) +
#   facet_grid(era ~ .,
#              scales = "free")

all_data %>%
  ggplot(aes(x=id,
             y=allele_balance,
             fill=era)) +
  geom_boxplot() +
  geom_hline(yintercept = c(1/8,
                            1/6,
                            1/5,
                            1/4),
             color="black",
             linetype="solid") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  labs(title = "Boxplots of Allele Balance",
       subtitle = "Individuals x Position, Most Granular, HLines at 1/8, 1/6, 1/5, 1/4",
       x = "Indiviudal ID") 

ggsave(paste(outDIR, 
             outFilePREFIX,
             'BOXPL-AB-INDxPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

all_data %>%
  ggplot(aes(x=id,
             y=log10(num_reads),
             fill=era)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  labs(title = "Boxplots of Read Depth",
       subtitle = "Individuals x Position, Most Granular",
       x = "Indiviudal ID") 

ggsave(paste(outDIR, 
             outFilePREFIX,
             'BOXP-DP-INDxPOS.png', 
             sep = ""), 
       height = 9, 
       width = 6.5)

# all_data %>%
#   filter(genotype == "hetero") %>%
#   group_by(id, era) %>%
#   summarize(n_hetero_pos = n()) %>%
#   ggplot(aes(x=id,
#              y=n_hetero_pos,
#              fill = era)) +
#   theme(axis.text.x = element_text(angle = 90, hjust=1)) +
#   geom_col()
# 
# all_data %>%
#   filter(is.na(genotype)) %>%
#   group_by(id, era) %>%
#   summarize(n_missdat_pos = n()) %>%
#   ggplot(aes(x=id,
#              y=n_missdat_pos,
#              fill = era)) +
#   theme(axis.text.x = element_text(angle = 90, hjust=1)) +
#   geom_col()

all_data %>%
  ggplot(aes(x=id,
             fill = genotype)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) 
  

ggsave(paste(outDIR, 
             outFilePREFIX,
             'BARPL-GT-INDxPOS.png', 
             sep = ""), 
       height = 6.5, 
       width = 9)

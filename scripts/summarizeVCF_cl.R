#### Initialize ####
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

list.of.packages <- c("vcfR", 
                      "tidyverse",
                      "janitor",
                      "ggforce")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(vcfR)
library(tidyverse)
library(janitor)
library(ggforce)


#### USER DEFINED VARIABLES ####
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
  vcfPATH = "Pfalcifer.Hlepto.ind.3.7.Fltr20.1.randSNPperLoc.neutral.recode.vcf"
  # popMapPATH = ""
} else {
  vcfPATH = args[1]
  # popMapPATH = args[2]
}

#### OTHER VARIABLES ####
outPATH = str_remove(vcfPATH,
                     "\\.vcf$")

#### PDF ####
pdf(str_c(outPATH,
          ".pdf"),
           width = 9,
           height = 6.5)

#### READ IN VCF ####

vcf <- read.vcfR(vcfPATH)

#### WRANGLE GENOTYPES ####
# gt <- 
#   extract.gt(vcf) %>%
#   as_tibble(rownames = "chrom_pos") %>%
#   separate(chrom_pos,
#            into = c("ddocent",
#                     "contig",
#                     "contig_id",
#                     "pos"),
#            sep = "_",
#            remove=FALSE) %>%
#   mutate(chrom= str_c(ddocent,contig,contig_id,sep="_")) %>%
#   rename_with(.,
#               ~ str_remove(.,
#                          "\\.3\\.7"),
#               contains(".3.7")) %>%
#   select(chrom_pos,
#          chrom,
#          ddocent:pos,
#          At_Pfa034:last_col()) %>%
#   pivot_longer(cols = At_Pfa034:last_col(),
#                names_to = "individual",
#                values_to = "gt") %>%
#   # colors from Fig 2 of Bird et al, STructure & PCA
#   mutate(genetic_cluster = case_when(individual == "St_Pfa017" ~ "purple",
#                                      individual == "Kr_Pfa002" ~ "purple",
#                                      individual == "Kr_Pfa006" ~ "purple",
#                                      individual == "Kr_Pfa009" ~ "purple",
#                                      individual == "Kr_Pfa017" ~ "purple",
#                                      individual == "Kr_Pfa019" ~ "purple",
#                                      individual == "Kr_Pfa023" ~ "purple",
#                                      individual == "Kr_Pfa026" ~ "purple",
#                                      individual == "Kr_Pfa028" ~ "purple",
#                                      individual == "Kr_Pfa034" ~ "purple",
#                                      individual == "Kr_Pfa036" ~ "purple",
#                                      individual == "Kr_Pfa038" ~ "purple",
#                                      individual == "At_Pfa052" ~ "dark_blue",
#                                      individual == "At_Pfa055" ~ "dark_blue",
#                                      individual == "At_Pfa077" ~ "dark_blue",
#                                      individual == "At_Pfa086" ~ "dark_blue",
#                                      individual == "At_Pfa089" ~ "dark_blue",
#                                      individual == "At_Pfa095" ~ "dark_blue",
#                                      str_detect(individual, "At_") ~ "light_blue",
#                                      TRUE ~ "dark_blue"),
#          site = str_remove(individual,
#                            "_Pfa.*"))

#### WRANGLE DP ####

vcf_tidy_0 <-
  vcfR2tidy(vcf)

remove(vcf)

vcf_tidy <-
  vcf_tidy_0$gt %>%
  left_join(vcf_tidy_0$fix) %>%
  clean_names() %>%
  mutate(indiv = str_remove(indiv,
                            "\\.3\\.7")) %>%
  rename(individual = indiv) 
  # left_join(gt %>%
  #             select(-ddocent:-contig) %>%
  #             mutate(pos = as.numeric(pos)),
  #           by = c("pos",
  #                  "chrom",
  #                  "gt_gt" = "gt",
  #                  "individual")) %>%
  # select(chrom_key,
  #        chrom_pos,
  #        chrom,
  #        contig_id,
  #        pos,
  #        genetic_cluster:site,
  #        individual,
  #        gt_gt:gt_gt_alleles,
  #        id:an)

remove(vcf_tidy_0)

#### WRANGLE MISSINGNESS ####
vcf_tidy %>%
  # filter(!is.na(gt_gt)) %>%
  # filter(individual != "Kr_Pfa015",
  #        individual != "Kr_Pfa024",
  #        individual != "Kr_Pfa043",
  #        individual != "Kr_Pfa046",
  #        individual != "St_Pfa004",
  #        individual != "St_Pfa006",
  #        individual != "St_Pfa007",
  #        individual != "St_Pfa018",
  #        individual != "St_Pfa024",
  #        individual != "St_Pfa041",
#        individual != "St_Pfa047",
#        individual != "St_Pfa048",
#        individual != "St_Pfa010",
#        individual != "St_Pfa013",
#        individual != "St_Pfa019",) %>%
mutate(gt_gt = replace_na(gt_gt, 
                          "missing")) %>%
  group_by(chrom_key,
           chrom, 
           pos,
           gt_gt) %>%
  summarize(n = n()) %>%
  # mutate(n = replace_na(n, 0)) %>% 
  pivot_wider(names_from = gt_gt,
              values_from = n) %>%
  replace_na(list(`0/0` = 0,
                  `0/1` = 0,
                  `1/1` = 0,
                  missing = 0)) %>%
  mutate(prop_gt_missing = missing / (`0/0` + `0/1` + `1/1` + missing)) ->
  prop_missing_by_site

vcf_tidy %>%
  # filter(!is.na(gt_gt)) %>%
  group_by( # genetic_cluster,
    # site,
    individual,
    gt_gt) %>%
  summarize(n = n(),
            dpmean = mean(gt_dp),
            dpsd = sd(gt_dp)) %>%
  pivot_wider(names_from = gt_gt,
              values_from = n:dpsd) %>%
  mutate(prop_gt_missing = n_NA / (`n_0/0` + `n_0/1` + `n_1/1` + n_NA)) ->
  prop_missing_by_individual


#### VISUALIZE GENOTYPES ####

vcf_tidy %>%
  # filter(is.na(gt)) %>%
  ggplot(aes(x=individual,
             fill = gt_gt)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  labs(title = "All Individuals",
       x = "Individual ID",
       y = "Number of SNP Loci") 

# ggsave(str_c(outPATH,
#              '_counts_gt-snp-ind.png'),
#        width = 9,
#        height = 6.5,
#        units = "in")

vcf_tidy %>%
  # filter(is.na(gt)) %>%
  left_join(prop_missing_by_individual %>%
              select(individual,
                     prop_gt_missing) %>%
              rename(prop_gt_missing_indiv = prop_gt_missing)) %>%
  filter(prop_gt_missing_indiv > 0.05) %>%
  ggplot(aes(x=individual,
             fill = gt_gt)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  labs(title = "Individuals w/ More than 5% Missing Genotypes",
       x = "Individual ID",
       y = "Number of SNP Loci") 

vcf_tidy %>%
  # filter(is.na(gt)) %>%
  left_join(prop_missing_by_individual %>%
              select(individual,
                     prop_gt_missing) %>%
              rename(prop_gt_missing_indiv = prop_gt_missing)) %>%
  filter(prop_gt_missing_indiv > 0.10) %>%
  ggplot(aes(x=individual,
             fill = gt_gt)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  labs(title = "Individuals w/ More than 10% Missing Genotypes",
       x = "Individual ID",
       y = "Number of SNP Loci") 

vcf_tidy %>%
  # filter(is.na(gt)) %>%
  left_join(prop_missing_by_individual %>%
              select(individual,
                     prop_gt_missing) %>%
              rename(prop_gt_missing_indiv = prop_gt_missing)) %>%
  filter(prop_gt_missing_indiv > 0.25) %>%
  ggplot(aes(x=individual,
             fill = gt_gt)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  labs(title = "Individuals w/ More than 25% Missing Genotypes",
       x = "Individual ID",
       y = "Number of SNP Loci") 

vcf_tidy %>%
  # filter(is.na(gt)) %>%
  left_join(prop_missing_by_individual %>%
              select(individual,
                     prop_gt_missing) %>%
              rename(prop_gt_missing_indiv = prop_gt_missing)) %>%
  filter(prop_gt_missing_indiv > 0.5) %>%
  ggplot(aes(x=individual,
             fill = gt_gt)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  labs(title = "Individuals w/ More than 50% Missing Genotypes",
       x = "Individual ID",
       y = "Number of SNP Loci") 


vcf_tidy %>%
  group_by(individual,
           gt_gt) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from = "gt_gt",
              values_from = "n") %>%
  mutate(missing_gt_prop = `NA` / sum(`0/0`,
                                      `0/1`,
                                      `1/1`,
                                      `NA`)) %>%
  ggplot(aes(x=missing_gt_prop)) +
  geom_histogram() 

# ggsave(str_c(outPATH,
#              '_hist_propmissgt-ind.png'),
#        width = 9,
#        height = 6.5,
#        units = "in")




#### VISUALIZE MISSING DATA BY SNP LOCUS ####
prop_missing_by_site %>%
  ggplot(aes(x = prop_gt_missing)) +
  geom_histogram() +
  labs(y = "Number of SNP Loci",
       x = "Proportion of Genotypes Filtered")

prop_missing_by_site %>%
  mutate(chrom_pos = str_c(chrom,
                           pos,
                           sep = "_")) %>%
  ggplot(aes(x=fct_reorder(chrom_pos,
                           prop_gt_missing),
             y=prop_gt_missing)) +
  geom_col() +
  theme(axis.text.x=element_blank()) +
  labs(y = "Proportion of Genotypes Filtered or Missing Per SNP Locus",
       x = "SNP Locus Sorted by Proportion of Genotypes Filtered or Missing")

#### HISTOGRAMS OF COVERAGE ####

vcf_tidy %>%
  # filter(!is.na(gt_gt)) %>%
  ggplot(aes(x=gt_dp,
             fill = gt_gt)) +
  geom_histogram(bins=50,
                 alpha = 0.7,
                 position = "identity") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  scale_x_continuous(trans='log10') +
  labs(title = "Histogram of Coverage Depth Across All Libraries",
       x = "Depth of Coverage",
       y = "Number of Genotypes",
       fill = "Genotype") 
  # facet_wrap(site ~ genetic_cluster,
  #            scales = "free")
# ggsave('puntioplites_dp_site-cluster.png',
#        width = 9,
#        height = 6.5,
#        units = "in")

vcf_tidy %>%
  # filter(!is.na(gt_gt)) %>%
  ggplot(aes(x=gt_dp,
             fill = individual)) +
  geom_histogram(bins=50,
                 alpha = 0.7,
                 position = "identity") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  scale_x_continuous(trans='log10') +
  labs(title = "Histogram of Coverage Depth Across All Libraries",
       x = "Depth of Coverage",
       y = "Number of Genotypes",
       fill = "Library") 
# facet_wrap(site ~ genetic_cluster,
#            scales = "free")

vcf_tidy %>%
  pull(individual) %>%
  unique() %>%
  length() ->
  num_individuals

num_pages = ceiling(num_individuals / (8*4))

vcf_tidy %>%
  # filter(!is.na(gt_gt),
  #        genetic_cluster == "dark_blue") %>%
  ggplot(aes(x=gt_dp,
             fill = gt_gt)) +
  geom_histogram(bins=50,
                 alpha = 0.8,
                 position = "identity") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  scale_x_continuous(trans='log10') +
  labs(title = "",
       fill = "Genotype",
       x = "Depth of Coverage",
       y = "Number of Genotypes") ->
  p
  
  for(i in 1:num_pages){
    print(
      p + 
        facet_wrap_paginate(individual ~ .,
                            nrow = 4,
                            ncol = 8,
                            page = i))
      Sys.sleep(2)
      # ggsave('puntioplites_dp_indiv_1.png',
      #        width = 9,
      #        height = 6.5,
      #        units = "in")
  }

#### BOXPLOTS OF COVERAGE BY GENOTYPE ####
 
vcf_tidy %>%
  filter(!is.na(gt_gt)) %>%
  
  ggplot(aes(x=gt_gt,
             y=gt_dp,
             fill = gt_gt)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5)) +
  scale_y_continuous(trans='log10')
# facet_wrap(site ~ genetic_cluster,
#            scales = "free")

vcf_tidy %>%
  filter(!is.na(gt_gt)) %>%
  ggplot(aes(x=individual,
             y=gt_dp,
             fill = gt_gt)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5))+
  scale_y_continuous(trans='log10')
  # facet_wrap(site ~ genetic_cluster,
  #            scales = "free")

#### MEAN COVERAGE VS PROP GENOTYPE MISSING####
vcf_tidy %>%
  # filter(!is.na(gt_gt)) %>%
  group_by( # genetic_cluster,
           # site,
           individual,
           gt_gt) %>%
  summarize(n = n(),
            dpmean = mean(gt_dp),
            dpsd = sd(gt_dp)) %>%
  pivot_wider(names_from = gt_gt,
              values_from = n:dpsd) %>%
  mutate(prop_gt_missing = n_NA / (`n_0/0` + `n_0/1` + `n_1/1` + n_NA)) %>%
  select(-starts_with("n_"),
         -contains("_NA")) %>%
  pivot_longer(cols = `dpmean_0/0`:`dpsd_1/1`,
               names_to = c("metric",
                            "gt_gt"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_wider(names_from = metric,
              values_from = value) %>%
  ggplot(aes(x=prop_gt_missing,
             y=dpmean,
             color = gt_gt)) +
  geom_point() +
  geom_errorbar(aes(ymin = dpmean-dpsd,
                    ymax = dpmean+dpsd)) +
  geom_smooth(se=FALSE) +
  scale_y_continuous(trans='log10') +
  # ylim(0,1000) +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust=0.5))
  # facet_wrap(site ~ genetic_cluster)

#### missing data ####
# vcf_tidy %>%
#   group_by( # genetic_cluster,
#            # site,
#            individual,
#            gt_gt) %>%
#   summarize(n = n(),
#             dpmean = mean(gt_dp),
#             dpsd = sd(gt_dp)) %>%
#   pivot_wider(names_from = gt_gt,
#               values_from = n:dpsd) %>%
#   mutate(pct_gt_missing = n_NA / (`n_0/0` + `n_0/1` + `n_1/1` + n_NA)) %>%
#   select(-starts_with("n_"),
#          -contains("_NA")) %>%
#   pivot_longer(cols = `dpmean_0/0`:`dpsd_1/1`,
#                names_to = c("metric",
#                             "gt_gt"),
#                names_pattern = "(.*)_(.*)") %>%
#   pivot_wider(names_from = metric,
#               values_from = value)

#### END PDF ####
dev.off()
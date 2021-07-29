#process output from vcftools --het --missing-indv --missing-site
#assumes ddocent naming format, pop_ind

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### libraries ####
library(tidyverse)
library(janitor)
library(stringr)

#### user-defined variables ####
hetFILE = "out.het"
imissFILE = "out.imiss"
lmissFILE = "out.lmiss"

#### read in data ####
het <-
  read_tsv(hetFILE) %>%
  clean_names() %>%
  mutate(o_het = n_sites - o_hom,
         heterozygosity_obs = o_het/n_sites) %>%
  separate(indv,
           c("site",
             NA),
           sep = "_",
           remove = "FALSE")

imiss <-
  read_tsv(imissFILE) %>%
  clean_names() %>%
  rename(prop_missing_data = f_miss) %>%
  separate(indv,
           c("site",
             NA),
           sep = "_",
           remove = "FALSE")

lmiss <-
  read_tsv(lmissFILE) %>%
  clean_names() %>%
  rename(prop_missing_data = f_miss) %>%
  separate(indv,
           c("site",
             NA),
           sep = "_",
           remove = "FALSE")



#### visualize het data ####
pdf("het.pdf")
het %>%
  ggplot(aes(x = indv,
             y = heterozygosity_obs,
             fill = heterozygosity_obs)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title="Heterozygosity")

het %>%
  arrange(heterozygosity_obs) %>%
  ggplot(aes(x = indv,
             y = heterozygosity_obs,
             color = heterozygosity_obs)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title="Heterozygosity")

het %>%
  ggplot(aes(y = heterozygosity_obs,
             x = site)) +
  geom_boxplot() +
  labs(title="Heterozygosity")
dev.off()

#### visualize imiss data ####
pdf("imiss.pdf")
imiss %>%
  ggplot(aes(x = indv,
             y = prop_missing_data,
             fill = prop_missing_data)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title="Missing Data By Individual")

imiss %>%
  arrange(prop_missing_data) %>%
  ggplot(aes(x = indv,
             y = prop_missing_data,
             color = prop_missing_data)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title="Missing Data By Individual")

imiss %>%
  ggplot(aes(y = prop_missing_data,
             x = site)) +
  geom_boxplot() +
  labs(title="Missing Data By Individual")
dev.off()

#### visualize lmiss data ####
pdf("lmiss.pdf")
lmiss %>%
  arrange(prop_missing_data) %>%
  ggplot(aes(x = chr,
             y = prop_missing_data,
             color = prop_missing_data)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Missing Data By Locus")

lmiss %>%
  ggplot(aes(y = prop_missing_data)) +
  geom_boxplot() +
  labs(title="Missing Data By Locus")
dev.off()

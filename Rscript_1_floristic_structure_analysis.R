### R scripts for analyses in:
### Aticle: Grade of Membership models reveal geographical and environmental correlates of floristic structure in a temperate biodiversity hotspot
### Authors: Qin Li, Hang Sun, David E. Boufford, Bruce Bartholomew, Peter W. Fritsch, Jiahui Chen, Tao Deng, Richard H. Ree
### DOI:
### compiled by Qin Li 20210329

#------------------------------------------------------------------------------------------------
# analyses include:
# 1. fit grade of membership model via ecostructure package
# 2. calculate motif overlap over K (including stability & overlap)
# 3. conduct PCA with environmental variables for K=9
# 4. calculate species/genus contribution and similarity network for K=9
# 5. estimate variable importance of motif correlates by random forest for K=c(4,6,9,HD)
#------------------------------------------------------------------------------------------------


#setwd("your/path")
source("Rscript_0_functions.R") # custom functions

#------------------------------------------------------------------------------------------------
# 1. fit Grade of Membership model (GoM)
# refer to R package ecostructure application at https://kkdey.github.io/ecostructure/
#------------------------------------------------------------------------------------------------
library(ecostructure)
library(parallel)
library(tidyverse)

# load species-by-site matrix
df = read.csv("data/species_by_site_matrix.csv") # and turn the 1st column of site to row names
df = df %>% column_to_rownames(var="site_id")
k_max = 20
Kvec = c(2:k_max) # K values: prior settings of motif number

# function to run GoM for parallel computing
ecos_fit_K = function(n){
    # K = n, the number of the motif to model
    pres_ab_fit = ecos_fit(df, K = n, tol = 10, num_trials = 100)
    saveRDS(pres_ab_fit, file=paste0("data/gom_fit_K",n,".rds"))
}

print("start model fitting:")
# initiate cluster
cl = makeCluster(length(Kvec), type="FORK") # or set a reasonable value for cluster
parLapply(cl, Kvec, ecos_fit_K)
stopCluster(cl) # shutdown cluster
print("Finish! YEAH!!!")


#------------------------------------------------------------------------------------------------
# 2. calculate motif overlap over K (including stability & overlap)
#------------------------------------------------------------------------------------------------

K = 2:12
fitK = lapply(paste0("data/gom_fit_K",K,".rds"), readRDS)
# extract motif contribution matrix (omega)
motif_list = lapply(fitK, function(x) x$omega %>% as.data.frame() %>% 
                         rename_all(list(~paste0("m_", .))))

# overlap between motifs across K by Jaccard index (see function cal.overlap())
motif_bi_list = list()
for(i in 1:length(motif_list)){
    print(i)
    
    df = motif_list[[i]]
    df_bi = df
    df_bi[] = 0
    
    # assigning a site with the dominant motif
    df$group = sapply(1:nrow(df), function(x) names(df)[which.max(df %>% slice(x))])
    
    for(g_id in 1:nrow(df)){
        row2gp = df$group[g_id]
        df_bi[g_id, which(names(df_bi) == row2gp)] = 1
    }
    motif_bi_list[[i]] = df_bi
}

K_comb = data.frame(K1=2:11, K2=3:12)
overlap_list = list()
for(i in 1:nrow(K_comb)){
    print(i)
    overlap_list[[i]] = cal.overlap(motif_bi_list=motif_bi_list, K_no=as.numeric(K_comb[i,]))
}

#saveRDS(overlap_list, file="data/comb_motif_overlap.rds")

#------------------------------------------------------------------------------------------------
# 3. environmental PCA for K = 9
#------------------------------------------------------------------------------------------------
library(ade4)

# load model fitting results
K = 9
fitK = readRDS(file=paste0("data/gom_fit_K",K,".rds"))
motif_df = as.data.frame(fitK$omega) %>% rename_all(list(~paste0("motif_", .)))

# extract the dominant motif for each site
max_motif = data.frame(ctEB = row.names(motif_df), group=NA, stringsAsFactors = F)
max_motif$group = sapply(1:nrow(motif_df), function(x) names(motif_df)[which.max(motif_df %>% slice(x))])
max_motif$m_contri = sapply(1:nrow(max_motif), function(x) motif_df[x, which(names(motif_df) == max_motif$group[x])])

# load environmental variables
# P_annual & P_cqtr & Roughness: log-transformed
var_ctEB_avg = read.csv("data/site_enVar_avg.csv",header=T)
var_sel <- c("P_annual","T_mean","SR","T_season","P_season","T_cold","P_cqtr","Roughness")
motif_var_m29 = var_ctEB_mu %>% left_join(max_motif) %>% filter(group != "motif_1")

# PCA with motif contribution as weight
motif_var_m29$m_weight = motif_var_m29$m_contri / sum(motif_var_m29$m_contri)
m29_pca = dudi.pca(motif_var_m29[,var_sel], nf=length(var_sel), center=T, scale=T, scannf=F, 
                   row.w = motif_var_m29$m_weight)
#saveRDS(m29_pca, file="data/m29_var_pca.rds")


#------------------------------------------------------------------------------------------------
# 4. species composition and similarity network: for K = 9
#------------------------------------------------------------------------------------------------
# load species occ. prob. & exclusivity
K = 9
fit_K9 = readRDS(file=paste0("data/gom_fit_K",K,".rds"))
sp_prob = as.data.frame(fit_K9$theta) %>% select(-1) %>% 
    rename_all(list(~paste0("motif_", .)))
sp_prob = sp_prob[-which(rowSums(sp_prob) == 0),]
sp_exclusivity = sp_prob / rowSums(sp_prob) # based on equation-1

# taxonomy rank: species-genus
taxa_df = sp_prob %>% rownames_to_column(var="species") %>% 
    mutate(genus=word(species,1,1,sep=fixed(" "))) %>% 
    select(species, genus)

# load motif id and name
motif_name_dt = read.csv(file="data/motif_name_K_dtset.csv")
motif_name = motif_name_dt %>% filter(dt_set == "K9") %>% 
    mutate(motif_id = as.numeric(sub("motif_","",motif)))

# calculate the contribution of a species to a motif based on equation-2
exc_thd = 0.1
prob_thd = 0.1
sp_contrib = vector("list", length = 9) # the contribution

# extract top five genera based on summed magnitude/contribution and save into lists
genus_order = vector("list", length = 9) # top five genera with their summed contribution
motif_taxa_g5 = vector("list", length = 9) # species in top five genera with their summed contribution

for(n in 2:K){
    #print(n)
    m_id = n
    print(motif_name[which(motif_name$motif_id == n),])
    
    m_col = paste0("motif_",n)
    
    temp_df = sp_exclusivity %>% select(one_of(m_col)) %>% 
        rename(exc=1) %>% tibble::rownames_to_column(var="species") %>% 
        left_join(sp_prob %>% select(one_of(m_col)) %>% 
                      rename(prob=1) %>% tibble::rownames_to_column(var="species")) %>% 
        filter(exc > exc_thd & prob > prob_thd)
    temp_df$sp_contrib = sqrt(temp_df$exc^2 + temp_df$prob^2) # based on equation-2
    
    temp_df$motif_name = motif_name$motif_name[which(motif_name$motif_id == n)]
    sp_contrib[[m_id]] = temp_df
    
    # extract species
    motif_taxa = sp_contrib[[m_id]] %>% select(species, sp_contrib) %>% left_join(taxa_df)
    
    genus_order[[m_id]] = motif_taxa %>% group_by(genus) %>% 
        summarise(genus_contrib = sum(sp_contrib)) %>% 
        left_join(motif_taxa %>% group_by(genus) %>% summarise(no_sp=n())) %>% 
        arrange(desc(genus_contrib)) 
    
    genus_order[[m_id]] = genus_order[[m_id]] %>% 
        filter(genus_contrib >= genus_order[[m_id]]$genus_contrib[5])
    genus_order[[m_id]]$g_order = paste0("g_",1:nrow(genus_order[[m_id]]))
    genus_order[[m_id]]$motif_id = m_id
    genus_order[[m_id]]$motif_name = motif_name$motif_name[which(motif_name$motif_id_new == m_id)]
    
    motif_taxa_g5[[m_id]] = motif_taxa %>% filter(genus %in% genus_order[[m_id]]$genus) %>% 
        select(species, genus, sp_contrib) %>% mutate(motif_id = m_id)
}
# top five genera with their contrbutions in each motif and number of species in each genus
genus_order_list = do.call(rbind,genus_order)
# extract genera with overlap line among motifs
g2over = genus_order_list %>% group_by(genus) %>% summarise(m_no=n()) %>% 
    as.data.frame() %>% filter(m_no > 1) %>% arrange(m_no)

# calculate overlap between motifs for each genus based on equation-3
g_overlap = list()
for(j in 1:nrow(g2over)){
    g_test = g2over$genus[j] # genus to calculate overlap between motifs
    
    # motif id with the genus in (with updated motif id)
    m2over = genus_order_list %>% filter(genus == g_test) %>% pull(motif_id)
    
    # combination matrix for motifs
    g_comb = as.data.frame(t(combn(m2over, 2)))
    g_comb$genus = g_test
    g_comb$sp_overlap = 0
    
    for(i in 1:nrow(g_comb)){
        sp_v1 = motif_taxa_g5[[g_comb$V1[i]]] %>% filter(genus == g_test)
        sp_v2 = motif_taxa_g5[[g_comb$V2[i]]] %>% filter(genus == g_test)
        sp_contrib = sp_v1 %>% select(-motif_id) %>% rename(m1_contrib = sp_contrib) %>% 
            full_join(sp_v2 %>% select(-motif_id) %>% rename(m2_contrib = sp_contrib))
        sp_contrib[is.na(sp_contrib)] = 0
        sp_contrib$min_contrib = ifelse(sp_contrib$m1_contrib >= sp_contrib$m2_contrib,
                                     sp_contrib$m2_contrib, sp_contrib$m1_contrib)
        max_contrib = max(sp_contrib$min_contrib)
        g_comb$sp_overlap[i] = sum(sp_contrib$min_contrib/max_contrib)/nrow(sp_contrib) # based on equation-3
    }
    g_overlap[[j]] = g_comb
}
g_overlap_df = do.call(rbind,g_overlap) %>% as.data.frame() %>% 
    arrange(genus,desc(sp_overlap)) %>% filter(sp_overlap > 0)

#save(genus_order_list, g_overlap_df, file="data/K9_genus_overlap.RData")


#------------------------------------------------------------------------------------------------
# 5. random forest: variable importance of motif correlates
#------------------------------------------------------------------------------------------------
# run the same code for datasets with different K values

library(randomForest)
library(party)

# load site-motif-predictor data for random forest
# each row: a site with corresponding motif designation
# predictors: latitude, longitude, elevation, environmental PC1 & PC2
dt_set = c("K4","K6","K9","HD")
for(k in 1:length(dt_set)){
    motif_pred_df = read.csv(file=paste0("data/motif_pred_",dt_set[k],".csv"))
    # scale variables
    motif_pred_df = motif_pred_df %>% mutate(envPC1=scale(envPC1),envPC2=scale(envPC2),
                                             lat=scale(lat),lon=scale(lon),alt=scale(alt))
    
    # run random forest and calculate variable importance via conditional permutation
    rf_formula = as.formula("motif ~ envPC1 + envPC2 + lat + lon + alt")
    crf_res = list()
    cvi_res = list()
    crf_try = function(df){cforest(rf_formula, data=df, control=cforest_unbiased(mtry=2, ntree=1000))}
    for(i in 1:100){
        set.seed(123+i)
        crf_res[[i]]=crf_try(motif_pred_df)
        cvi_res[[i]]=varimp(crf_res[[i]],conditional=T)
    }
    
    # extract variable importance from 100 runs by cvi_extract()
    # and calculat the mean and sd
    cvi_df = cvi_extract(cvi_res)
    cvi_df = cvi_df %>% 
        group_by(var) %>% 
        summarise(cvi_mu = mean(cvi),
                  cvi_sd = sd(cvi)) %>% ungroup() %>% 
        mutate(dt_set = dt_set[k])
    
    #write.csv(cvi_df, file=paste0("data/cRF_varImp_",dt_set[k],".csv), row.names = F, quote = F)
}




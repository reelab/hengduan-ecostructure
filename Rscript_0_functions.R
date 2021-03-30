### R scripts of functions used in analyses in
### Aticle: Grade of Membership models reveal geographical and environmental correlates of floristic structure in a temperate biodiversity hotspot
### Authors: Qin Li, Hang Sun, David E. Boufford, Bruce Bartholomew, Peter W. Fritsch, Jiahui Chen, Tao Deng, Richard H. Ree
### DOI:
### compiled by Qin Li 20210329


# --- FUNCTIONS -----------------------------------------------

# function to calculate overlap between motifs
cal.overlap <- function(motif_bi_list, K_no=c(2:3)){
    
    df1_K <- K_no[1]
    df2_K <- K_no[2]
    
    df1 <- motif_bi_list[[c(df1_K-1)]]
    df2 <- motif_bi_list[[c(df2_K-1)]]
    
    m_comb <- expand.grid(df1_col=1:df1_K,df2_col=1:df2_K)
    m_comb$df1_motif <- paste0("K",df1_K,"_m",m_comb$df1_col)
    m_comb$df2_motif <- paste0("K",df2_K,"_m",m_comb$df2_col)
    
    for(j in 1:nrow(m_comb)){
        df_temp <- apply(cbind(df1[,m_comb$df1_col[j]], df2[,m_comb$df2_col[j]]), 1, sum)
        df_count <- as.data.frame(table(df_temp))
        no_AB <- ifelse("1" %in% df_count$df_temp,df_count$Freq[df_count$df_temp =="1"],0)
        no_C <- ifelse("2" %in% df_count$df_temp,df_count$Freq[df_count$df_temp =="2"],0)
        no_D <- ifelse("0" %in% df_count$df_temp,df_count$Freq[df_count$df_temp =="0"],0)
        
        m_comb$overlap_1[j] <- (no_C) / (no_AB + no_C) # Jaccard index (used in main analysis)
        m_comb$overlap_2[j] <- (no_D + no_C) / (no_D + no_AB + no_C)
        m_comb$overlap_3[j] <- (sqrt(no_C * no_D) + no_C) / (sqrt(no_C * no_D) + no_AB + no_C)
        
    }
    return(m_comb)
}

# function to switch motif id within each K for better visualization
m.switch <- function(df, m2s=c(1,2), K.no=2){
    m_1 <- paste0("K",K.no,"_m",m2s[1])
    m_2 <- paste0("K",K.no,"_m",m2s[2])
    df[df == m_1] <- "temp"
    df[df == m_2] <- m_1
    df[df == "temp"] <- m_2
    return(df)
}

# function to extract variable importance from multiple runs of random forest
cvi_extract <- function(vi_list){
    vi_matrix = matrix(NA, nrow = 5, ncol = length(vi_list))
    for(i in 1:length(vi_list)){
        vi_matrix[,i] = vi_list[[i]]
    }
    rownames(vi_matrix) = names(vi_list[[i]])
    vi_matrix = as.data.frame(vi_matrix) %>% rownames_to_column(var="var")
    return(vi_matrix)
}


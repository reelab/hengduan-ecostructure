### R scripts of visualization
### Aticle: Grade of Membership models reveal geographical and environmental correlates of floristic structure in a temperate biodiversity hotspot
### Authors: Qin Li, Hang Sun, David E. Boufford, Bruce Bartholomew, Peter W. Fritsch, Jiahui Chen, Tao Deng, Richard H. Ree
### DOI:
### compiled by Qin Li 20210329

# contents:
# 1. plot similarity path/nested motif
# 2. plot motif contribution to site map
# 3. plot PCA result
# 4. plot motif taxonomic similarity network
# 5. plot variable importance via random forest

#setwd("your/path")
source("Rscript_0_functions.R")

#------------------------------------------------------------------------------------------------
# 1. plot similarity path/nested motif
#------------------------------------------------------------------------------------------------
library(tidyverse)
library(ggnetwork)
### load motif overlap data
# overlap_list = readRDS("data/comb_motif_overlap.rds")
# overlap_list = do.call(rbind, overlap_list)
# overlap_list = overlap_list %>% filter(df1_col!= 1 & df2_col!= 1) %>% 
#     filter(!grepl("K2",df1_motif))
# 
### list motif id for later pairwise position adjust/switch for better visualization
# K_sw = c(6,7,7,8,10,10,11,11,11,11,12,12,12,12)
# m_sw_pair = list(c(3,4),c(3,4),c(4,5),c(6,7),
#                   c(6,7),c(8,9),c(5,7),c(9,11),
#                   c(4,6),c(5,6),c(4,5),c(5,8),
#                   c(9,11),c(10,12))
# names(m_sw_pair) = K_sw
# 
# # switch motif order if necessary
# for(i in 1:length(K_sw)){
#     print(paste("m2s =",paste(m_sw_pair[[i]],collapse = " & "), "- within K.no =",K_sw[i]))
#     overlap_list = m.switch(overlap_list, m2s=m_sw_pair[[i]], K.no=K_sw[i])
# }

overlap_list$overlap = overlap_list$overlap_1 # used Jaccard overlap index
overlap_thr = 0.1

pos_dat = rbind(overlap_list %>% select(df1_motif) %>% rename(node=df1_motif),
                 overlap_list %>% select(df2_motif) %>% rename(node=df2_motif)) %>% unique()
pos_dat$x = as.numeric(sub("m","",word(pos_dat$node,2,2,sep=fixed("_"))))
pos_dat$y = as.numeric(sub("K","",word(pos_dat$node,1,1,sep=fixed("_")))) * (-1)

# create another set of x coordinate to make the plot symmetric
pos_dat$x_2 = sapply(1:nrow(pos_dat), function(k) pos_dat$x[k]+(12+pos_dat$y[k])/2)

seg_dat = overlap_list %>% select(overlap, df1_motif, df2_motif) %>% 
    rename(from = df1_motif, to = df2_motif)
seg_dat = seg_dat %>% left_join(pos_dat %>% rename(from = node)) %>% # obtain starting xy
    left_join(pos_dat %>% rename(to = node, xend = x, yend = y, xend_2 = x_2)) # obtain ending xy

m_net = pos_dat %>% rename(from=node) %>% 
    mutate(to=from, xend=x, yend=y, xend_2=x_2, overlap=NA) %>% 
    rbind(seg_dat %>% arrange(overlap))
m_net$node_label = sub("m","",word(m_net$from,2,2,sep=fixed("_")))

# add node labels for K=4/6/9, otherwise use motif id
motif_name_dt = readRDS(file="data/motif_name_K_dtset.rds")
motif_node_label = motif_name_dt %>% mutate(from = paste0(dt_set,"_",motif)) %>% 
    mutate(from = sub("motif_","m",from)) %>% select(from, motif_name)

m_net = m_net %>% left_join(motif_node_label)
m_net$node_label = ifelse(is.na(m_net$motif_name), m_net$node_label, "")
m_net$motif_name = sub(" ","\n",m_net$motif_name)

m_net_plot = ggplot(m_net, aes(x = x_2, y = y, xend = xend_2, yend = yend)) +
    geom_edges(aes(size=overlap, colour=overlap), alpha=0.8, 
               data = function(x){x[which(x$overlap >= overlap_thr),]}) +
    scale_colour_gradient(low="gray90", high="gray20",
                          name="overlap\nindex",
                          limits=c(overlap_thr,1.0),
                          breaks=c(0.1,0.4,0.7,1.0)) +
    geom_nodes(size=rel(10), shape = 21, color = "gray10", fill = "white") +
    geom_nodes(size=rel(13), shape = 21, color = "gray10", fill = "white") +
    geom_nodetext(aes(label=node_label), fontface="bold") + 
    geom_nodetext(aes(label=motif_name), size = rel(4), fontface="bold") + 
    annotate("text",x=1,y=c(-3:-12), label=paste("K =",3:12),fontface = "bold") + 
    scale_size_continuous(guide = FALSE) +
    theme_blank() +
    guides(colour = guide_colourbar(title.position="top",ticks = FALSE, title.vjust=1))

#ggsave(filename="figures/motif_overlap_over_K.pdf",m_net_plot,width=10,height=9)


#------------------------------------------------------------------------------------------------
# 2. plot motif contribution to site map
#------------------------------------------------------------------------------------------------

library(raster)
library(rgdal)
# load the polygon shapefile of selected counties
ct_sel = readOGR(dsn=path.expand("data"),layer="county_shp.shp")
ct_sel_list = as.character(ct_sel@data$COUNTY_ID) # extract county id

# load DEM data downloaded (1-km resolution of NASA SRTM digital elevation model)
# with latitude: 17-50; longitude 72-112
#DEM_raster = raster(path/to/DEM/file)
raster_mask = DEM_raster
raster_mask[] = 0 # for later motif contribution plot

# define raster cells with corresponding county-by-elevation-band sites
elev_ct = list() # list to save each county's raster cell values
for(i in 1:length(ct_sel_list)){
    print(i)
    ct_one = ct_sel[row.names(ct_sel@data[which(ct_sel@data$COUNTY_ID==ct_sel_list[i]),]), ]
    elev_ct[[i]] = raster::extract(x=DEM_raster,y=ct_one,df=TRUE,cellnumbers=T)
    elev_ct[[i]]$county = ct_sel_list[i]
}

site_raster = bind_rows(elev_ct) %>% arrange(cell) 
coordDEM = coordinates(DEM_raster)[site_raster[,2],] # coordinates of cells
site_raster = cbind(site_raster, coordDEM)
site_raster$elevBand = sapply(site_raster$DEM_raster, function(x) paste0("eb_",ceiling(x/500)), simplify=T)
site_raster$ctEB = paste(site_raster$county, site_raster$elevBand, sep="_")

# assign motif contribution to raster cells
# cells belonging to the same site have the same values
for(K in c(4,6,9)){
    fitK = readRDS(paste0("data/gom_fit_K",K,".rds"))
    omega_df = fitK$omega %>% as.data.frame() %>% rename_all(list(~paste0("m_", .)))
    
    for(i in 1:K){
        print(i)
        motif_df = site_raster %>% 
            left_join(omega_df[i] %>% tibble::rownames_to_column(var="ctEB") %>% 
                          rename(m_contrib=2))
        motif_r = DEM_mask
        motif_r[motif_df$cell] = motif_df$m_contrib
        writeRaster(motif_r,filename=paste0("figures/raster_K",K,"_m",i,"_map.tif"),format="GTiff",overwrite=TRUE)
    }
}

# then raster files could be read in QGIS to produce combined motif maps for each K.


#------------------------------------------------------------------------------------------------
# 3. plot PCA result
#------------------------------------------------------------------------------------------------
library(ade4)
library(factoextra)
library(colorspace)

#m29_pca = readRDS(file="data/m29_var_pca.rds")
eig_v = round(m29_pca$eig/sum(m29_pca$eig) * 100, 1)

# variable loadings and contribution for first 2 PCs
#fviz_pca_var(m29_pca,col.var="contrib",gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)

motif_name_dt = readRDS(file="data/motif_name_K_dtset.rds")
motif_name = motif_name_dt %>% filter(dt_set="K9")

motif_var_m29 = motif_var_m29 %>% 
    left_join(motif_name %>% select(motif_name, motif) %>% rename(group=motif))
motif_var_m29$motif_name = factor(motif_var_m29$motif_name,
                                  levels=motif_name$motif_name[motif_name$motif_order])
motif_color = lighten(hue_pal()(8), 0.15) # or manually define colors

# plot confidence/concentration ellipses in normal probability
m29_pca_ellipse_plot = fviz_pca_biplot(m29_pca, axes = c(1,2),
                                       geom.ind = "point", alpha.ind = 1, pointsize = 1,
                                       habillage = motif_var_m29$motif_name, palette = motif_color,
                                       addEllipses=TRUE, ellipse.level = 0.9, ellipse.alpha = 0.01,
                                       label = "var", col.var = "grey20", repel = TRUE,
                                       xlab = paste0("PC 1 (", eig_v[1] ,"%)"), 
                                       ylab = paste0("PC 2 (", eig_v[2] ,"%)"),
                                       axes.linetype=NA,
                                       title = NULL) +
    coord_fixed(ratio = 1 ) +
    labs(fill=" motifs", shape=" motifs", colour=" motifs") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid = element_blank())

#ggsave(filename="figures/K9_niche_8var.pdf",m29_pca_ellipse_plot,width=8,height=6)


#------------------------------------------------------------------------------------------------
# 4. plot taxonomic similarity network
#------------------------------------------------------------------------------------------------
library(scales)

load("data/K9_genus_overlap.RData") # genus_order_list, g_overlap_df

# order genus by species richness for each motif, and assign x-axis coordinate
genus_order_list = genus_order_list %>% 
    arrange(motif_id, g_order) %>% group_by(motif_id) %>% 
    mutate(xright = cumsum(no_sp)) %>% as.data.frame() %>% 
    group_by(motif_id) %>% 
    mutate(xleft = lag(xright)) %>% as.data.frame()
genus_order_list[is.na(genus_order_list)] = 0

# extract genus list and assign colors
motif_g5 = genus_order_list %>% select(motif_id,xleft,xright,g_order,genus,motif_name)
genus_col = data.frame(genus=unique(motif_g5$genus), genus_col=NA, stringsAsFactors=F)
genus_col$genus_col = hue_pal()(nrow(genus_col)) # or mannually define colors
motif_g5 = motif_g5 %>% left_join(genus_col)

# calculate link line position/x-axis coordinates
dist_vec = c(0:4,-4:-1)
dist_m = matrix(NA, ncol=9, nrow=9, dimnames=list(c(2:10),c(2:10)))
dist_m[1,] = dist_vec
for(i in 2:nrow(dist_m)){
    dist_m[i,] = c(dist_vec[(length(dist_vec)-i+2):length(dist_vec)], dist_vec[1:(length(dist_vec)-i+1)])
}

dist_m_level = as.character(c(-1:-4, 4:1)) # for ordering link lines

g_overlap_df = g_overlap_df %>% arrange(sp_overlap)
g_overlap_df2 = g_overlap_df %>% mutate(V1_x = NA, V2_x = NA)

# generate link matrix and its x-axis coordinates
g_v_freq = g_overlap_df %>% select(1,3) %>% rename(v =V1) %>% 
    rbind(g_overlap_df %>% select(2,3) %>% rename(v =V2)) %>% 
    group_by(genus,v) %>% summarise(v_freq = n()) %>% as.data.frame()

for(i in 1:nrow(g_v_freq)){
    g_temp = g_v_freq$genus[i]
    m_temp = g_v_freq$v[i]
    freq_link = g_v_freq$v_freq[i]
    
    g_x = motif_g5 %>% filter(genus == g_temp & motif_id == m_temp) %>% select(xleft,xright)
    g_x_breaks = seq(g_x[1,1],g_x[1,2], length.out = 2+freq_link)
    
    row_id =  which(g_overlap_df2$genus == g_temp & (g_overlap_df2$V1 == m_temp | g_overlap_df2$V2 == m_temp))
    
    motif_dist = sapply(row_id, function(x) ifelse(g_overlap_df2$V1[x] == m_temp,
                                                    dist_m[which(row.names(dist_m)==m_temp),
                                                           which(names(dist_m)==g_overlap_df2$V2[x])],
                                                    dist_m[which(row.names(dist_m)==m_temp),
                                                           which(names(dist_m)==g_overlap_df2$V1[x])]))
    x_order = rank(as.numeric(factor(as.character(motif_dist), levels = dist_m_level)))
    
    # check rows for V1_x or V2_x, should be the same length as freq_link
    # xleft = breaks, xright = xleft + 1
    for(j in 1:length(row_id)){
        
        if(g_overlap_df2$V1[row_id[j]] == m_temp){
            g_overlap_df2$V1_x[row_id[j]] = g_x_breaks[-c(1,length(g_x_breaks))][x_order[j]]
        }
        
        if(g_overlap_df2$V2[row_id[j]] == m_temp){
            g_overlap_df2$V2_x[row_id[j]] = g_x_breaks[-c(1,length(g_x_breaks))][x_order[j]]
        }
    }
}

# create link lines with coordinates
link_df1 = g_overlap_df2 %>% 
    rename(xleft = V1_x) %>% 
    mutate(motif_id = as.character(V1)) %>% 
    mutate(xright = xleft + 100 * sp_overlap) %>% 
    select(motif_id,xleft,xright,sp_overlap,genus)

link_df2 = g_overlap_df2 %>% 
    rename(xleft = V2_x) %>% 
    mutate(motif_id = as.character(V2)) %>% 
    mutate(xright = xleft + 100 * sp_overlap) %>% 
    select(motif_id,xleft,xright,sp_overlap,genus)

# use subset with overlap > 0.01 for plotting
link_df1_sub = link_df1 %>% filter(sp_overlap > 0.01)
link_df2_sub = link_df2 %>% filter(sp_overlap > 0.01)

# plot network
library(circlize)

motif_g5$value = as.numeric(factor(motif_g5$genus)) # for plot genomicTrack

# adjust line color: degree of gray
ind_overlap = round(link_df1_sub$sp_overlap * 100)
col_overlap = sapply(ind_overlap/max(ind_overlap),function(x) adjustcolor("gray20", alpha.f = x))

pdf(file = "figures/K9_motif_genus_similarity_network.pdf",width = 8, height = 8)
# general settings for motif gap
circos.par("gap.degree" = rep(5, length(unique(motif_g5$motif_id))))

# generate plot area
circos.initializeWithIdeogram(motif_g5, plotType = NULL)

# add genus label
circos.genomicLabels(motif_g5, labels.column = 5, padding = 0.0, cex = 0.8,
                     connection_height = convert_height(5, "mm"),
                     side = "outside", col = motif_g5$genus_col, 
                     line_col = motif_g5$genus_col, line_lwd = 3)

# add motif name
circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA, 
             panel.fun = function(x, y) {
                 chr = unique(motif_g5$motif_name[which(motif_g5$motif_id == CELL_META$sector.index)])
                 xlim = CELL_META$xlim
                 ylim = CELL_META$ylim
                 circos.text(mean(xlim), mean(ylim), 
                             chr, 
                             cex = 1, col="black",
                             facing = "inside", niceFacing = TRUE)
             })

# add genus rectangular
circos.genomicTrack(motif_g5, track.height = 0.05, bg.border = NA,
                    panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, 
                                           col = value$genus_col,
                                           border = "gray")
                    })

# add link lines (similarity)
circos.genomicLink(link_df1_sub, 
                   link_df2_sub, 
                   h.ratio = 0.7, border = NA, lwd = 3, col = col_overlap)
dev.off()
circos.clear()


#------------------------------------------------------------------------------------------------
# 5. plot variable importance via random forest
#------------------------------------------------------------------------------------------------
# combine variable importance data of different K values (column: dt_set)
# and re-order by the rank of variables

dt_set = c("K4","K6","K9","HD")
cvi_df_comb = lapply(paste0("data/cRF_varImp_",dt_set,".csv"), read.csv) %>% 
    do.call(rbind, .)

cvi_df_comb = cvi_df_comb %>% 
    mutate(var2 = reorder_within(var, desc(cvi_mu), dt_set))
p_cvi_rank = ggplot(cvi_ldf, aes(x=var2, y=cvi_mu)) +
    facet_wrap(~ dt_set, scales = "free", nrow=1) +
    geom_point(aes(col=var)) + 
    geom_errorbar(aes(ymin=cvi_mu-cvi_sd,ymax=cvi_mu+cvi_sd,col=var),width=0.2) +
    scale_x_reordered() + 
    xlab("") +
    ylab("Mean decrease in accuracy") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "figures/variable_importance.pdf", p_cvi_rank, width = 10, height = 3)


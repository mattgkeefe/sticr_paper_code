### Fig 1b -- main UMAP with subcluster identities
ls_synthesis = SetIdent(ls_synthesis, value = "subcluster_identity")
table(ls_synthesis@meta.data$subcluster_identity)
ls_synthesis@active.ident = factor(ls_synthesis@active.ident, levels=c(
  "RG",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "EX_IPCs",
  "ENs",
  "IN_IPCs",
  "IN_local",
  "IN_CGE",
  "IN_OB",
  "IN_MGE"
))
DimPlot(ls_synthesis, pt.size=0.5) +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=20))
image = DimPlot(ls_synthesis, pt.size=0.5) +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=20))
ggsave(file="full_umap_subclust_idents.png", plot=image, width=12, height=8)
ggsave(file="full_umap_subclust_idents.svg", plot=image, width=12, height=8)


### Fig 1c -- full UMAP by sample
DimPlot(ls_synthesis, pt.size=0.1, group.by="sample_id", raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=20), plot.title = element_blank())
image = DimPlot(ls_synthesis, pt.size=0.5, group.by="sample_id", raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=20), plot.title = element_blank())
ggsave(file="full_umap_sample_idents.png", plot=image, width=12, height=8)
ggsave(file="full_umap_sample_idents.svg", plot=image, width=12, height=8)

### Fig 1d -- full UMAP by pseudoage
young_col = '#ffc012'
old_col = '#00a318'
DimPlot(ls_synthesis, pt.size=3, group.by="age_binned", cols = c(young_col, old_col), raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
image = DimPlot(ls_synthesis, pt.size=3, group.by="age_binned", cols = c(young_col, old_col), raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
ggsave(file="full_umap_pseudoage_idents_revised.png", plot=image, width=12, height=8)
ggsave(file="full_umap_pseudoage_idents_revised.svg", plot=image, width=12, height=8)

### Fig 1e -- sample by cluster contribution
ls_synthesis@meta.data$subcluster_identity = factor(ls_synthesis@meta.data$subcluster_identity, levels=c(
  "RG",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "EX_IPCs",
  "ENs",
  "IN_IPCs",
  "IN_local",
  "IN_CGE",
  "IN_OB",
  "IN_MGE"
))
ggplot(ls_synthesis@meta.data, aes(x=sample_id, fill=subcluster_identity)) +
  geom_bar(position="fill") +
  labs(fill="Cell type") +
  ylab("Proportion of cluster") +
  xlab("Sample") +
  ggtitle("Cluster composition by sample") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(ls_synthesis@meta.data, aes(x=sample_id, fill=subcluster_identity)) +
  geom_bar(position="fill") +
  labs(fill="Cell type") +
  ylab("Proportion of cluster") +
  xlab("Sample") +
  ggtitle("Cluster composition by sample") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave('stacked_barchart_cluster_compositon_by_sample_all_samples.svg', image, width=10, height=8)

### Fig 1f -- cluster composition by sample age
ggplot(ls_synthesis@meta.data, aes(x=subcluster_identity, fill=factor(age_binned, levels=c(">GW20", "<=GW20")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c(old_col, young_col)) +
  labs(fill="Sample age") +
  ylab('Ratio of cluster') +
  xlab('Cell type') +
  ggtitle('Cluster composition by sample age') +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(ls_synthesis@meta.data, aes(x=subcluster_identity, fill=factor(age_binned, levels=c(">GW20", "<=GW20")))) +
  geom_bar(position="fill") +
  scale_fill_manual(values=c(old_col, young_col)) +
  labs(fill="Sample age") +
  ylab('Ratio of cluster') +
  xlab('Cell type') +
  ggtitle('Cluster composition by sample age') +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave('stacked_barchart_cluster_compositon_by_age_all_samples.svg', image, width=10, height=8)

ggplot(ls_synthesis@meta.data, aes(x=1, fill=factor(age_binned, levels=c(">GW20", "<=GW20")))) +
  geom_bar(position = "fill", width = 1) +
  xlim(0,11) +
  scale_fill_manual(values=c(old_col, young_col)) +
  theme_bw()
image = ggplot(ls_synthesis@meta.data, aes(x=1, fill=factor(age_binned, levels=c(">GW20", "<=GW20")))) +
  geom_bar(position = "fill", width = 1) +
  xlim(0,11) +
  scale_fill_manual(values=c(old_col, young_col)) +
  theme_bw()
ggsave('stacked_barchart_cluster_compositon_by_age_single_bar.svg', image, width=10, height=8)

### Fig 1g -- Progenitors and glia DimPlot with subcluster identities
ls_synthesis_prog@meta.data$subcluster_identity = factor(ls_synthesis_prog@meta.data$subcluster_identity, 
                                                         levels=c(
                                                          "bRG",
                                                          "tRG",
                                                          "Early OPCs",
                                                          "Astrocytes (dense smooth)",
                                                          "Astrocytes (dense bulbous)", 
                                                          "Dividing", 
                                                          "EX IPCs"
                                                          ))
ls_synthesis_prog = SetIdent(ls_synthesis_prog, value = "subcluster_identity")
DimPlot(ls_synthesis_prog, pt.size=1.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30, family="Helvetica"))
image = DimPlot(ls_synthesis_prog, pt.size=1.5) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30, family="Helvetica"))
ggsave("ls_synthesis_prog_umap_seurat_clusters_revision.png", image, width=12, height=8)
ggsave("ls_synthesis_prog_umap_seurat_clusters_revision.svg", image, width=12, height=8)

### Fig 1h -- Barchart of progenitors and glia subcluster identity by VZ/OSVZ
ls_synthesis_prog@meta.data$GZ_origin = factor(ls_synthesis_prog@meta.data$GZ_origin, 
                                                         levels=c("OSVZ","VZ"))

vz_col = "#05c4c7"
osvz_col = "#ff6666"
ls_synthesis_prog@meta.data$subcluster_identity = factor(ls_synthesis_prog@meta.data$subcluster_identity, 
                                                         levels=c("tRG", "EX IPCs", "Astrocytes (dense bulbous)", "bRG", "Dividing", "Early OPCs", "Astrocytes (dense smooth)"))
ggplot(subset(ls_synthesis_prog@meta.data, GZ_origin %in% c("VZ", "OSVZ")), aes(x=subcluster_identity, fill=GZ_origin)) +
  geom_bar(position="fill") +
  ggtitle('GZ contribution to glial subclusters') +
  scale_fill_manual(values=c(osvz_col, vz_col))+
  ylab("Proportion of cells in cluster") +
  xlab("Glial subcluster") +
  theme_bw() +
  labs(fill = "GZ") +
  theme(plot.title = element_text(hjust = 0.5, size=25), legend.text = element_text(size=15), legend.title = element_text(size=20),
        axis.text = element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
image = ggplot(subset(ls_synthesis_prog@meta.data, GZ_origin %in% c("VZ", "OSVZ")), aes(x=subcluster_identity, fill=GZ_origin)) +
  geom_bar(position="fill") +
  ggtitle('GZ contribution to glial subclusters') +
  scale_fill_manual(values=c(osvz_col, vz_col))+
  ylab("Proportion of cells in cluster") +
  xlab("Glial subcluster") +
  theme_bw() +
  labs(fill = "GZ") +
  theme(plot.title = element_text(hjust = 0.5, size=25), legend.text = element_text(size=15), legend.title = element_text(size=20),
        axis.text = element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
ggsave('ls_synthesis_prog_stacked_barchart_subcluster_comp_by_GZ_astrocyte_subtypes.svg', image, width=10, height=8)

### Fig 1i -- Progenitors and glia DimPlot with pseudoage
DimPlot(ls_synthesis_prog, pt.size=3, group.by="age_binned", cols = c(young_col, old_col), raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
image = DimPlot(ls_synthesis_prog, pt.size=3, group.by="age_binned", cols = c(young_col, old_col), raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
ggsave(file="prog_glia_umap_pseudoage_idents_revised.png", plot=image, width=12, height=8)
ggsave(file="prog_glia_umap_pseudoage_idents_revised.svg", plot=image, width=12, height=8, device='svg')

### Fig 1j -- DEGs between <=GW20 and >GW20 radial glia (REVISION)
de.markers.glia_no_clonal_restriction = FindMarkers(subset(ls_prog, subclust_idents_definitive_final %in% c("RG", "bRG", "aRG")), ident.1 = ">GW20", ident.2 = "<=GW20")
head(de.markers.glia_no_clonal_restriction, n=20)

de.markers.glia_no_clonal_restriction_filt = de.markers.glia_no_clonal_restriction[-grep("MT", rownames(de.markers.glia_no_clonal_restriction)),]
de.markers.glia_no_clonal_restriction_filt = de.markers.glia_no_clonal_restriction_filt[-grep("RPL", rownames(de.markers.glia_no_clonal_restriction_filt)),]
de.markers.glia_no_clonal_restriction_filt = de.markers.glia_no_clonal_restriction_filt[-grep("XIST", rownames(de.markers.glia_no_clonal_restriction_filt)),]
de.markers.glia_no_clonal_restriction_filt = de.markers.glia_no_clonal_restriction_filt[-grep("RPS", rownames(de.markers.glia_no_clonal_restriction_filt)),]
dim(de.markers.glia_no_clonal_restriction)
#[1] 14899     5
dim(de.markers.glia_no_clonal_restriction_filt)
#[1] 14566     5
head(de.markers.glia_no_clonal_restriction_filt)
ggplot(data=de.markers.glia_no_clonal_restriction_filt, aes(x=avg_log2FC, y=-log10(p_val_adj), label=rownames(de.markers.glia_no_clonal_restriction_filt)))+
  geom_point()+
  theme_minimal()+
  geom_text_repel(
    box.padding = 0.4
  )+
  ggtitle("Differential expression between <=GW20 glia (left) and >GW20 glia (right) -- tRG, oRG, and RG clone size >=3")+
  theme(plot.title = element_text(hjust = 0.5))

de.markers.glia_no_clonal_restriction_filt$plot_genenames = 'do_not_plot'
de.markers.glia_no_clonal_restriction_filt$genenames_to_plot = rownames(de.markers.glia_no_clonal_restriction_filt)
de.markers.glia_no_clonal_restriction_filt$plot_genenames[rownames(de.markers.glia_no_clonal_restriction_filt) %in% c(young_enriched_genes_for_volcano, old_enriched_genes_for_volcano)] = 'do_plot'
de.markers.glia_no_clonal_restriction_filt$genenames_to_plot[de.markers.glia_no_clonal_restriction_filt$plot_genenames %in% c('do_not_plot')] = NA
head(de.markers.glia_no_clonal_restriction_filt)

ggplot(data=de.markers.glia_no_clonal_restriction_filt, aes(x=avg_log2FC, y=-log10(p_val_adj), label=genenames_to_plot))+
  geom_point()+
  theme_minimal()+
  geom_label_repel(
    box.padding = 0.8, alpha=0.9
  )+
  ggtitle("Differential expression between <=GW20 glia (left) and >GW20 glia (right) -- tRG, oRG, and RG")+
  theme(plot.title = element_text(hjust = 0.5))
write.csv(de.markers.glia_no_clonal_restriction_filt, 'de_markers_no_clonal_restrict_for_volcano_ls_prog_rg_clusts_young_vs_old_filt_box.csv')


young_enriched_genes_volcano_filtered = c('PAX6', 'ROBO1', 'FEZF2', 'EMX1', 'DCX', 'EOMES', 'NEUROG2', 'NEUROD2', 'NEUROD6', 'PPP1R17')
old_enriched_genes_volcano_filtered = c('OLIG2', 'HES5', 'HES1', 'S100B', 'SOX9', 'NR2F2', 'NR4A1', 'SPARCL1', 'GAD1', 'GAD2', 'METTL7B', 'DLX1', 'GJA1', 'DLX2', 'AQP4')
de.markers.glia_no_clonal_restriction_filt$plot_genenames = 'do_not_plot'
de.markers.glia_no_clonal_restriction_filt$genenames_to_plot = rownames(de.markers.glia_no_clonal_restriction_filt)
de.markers.glia_no_clonal_restriction_filt$plot_genenames[rownames(de.markers.glia_no_clonal_restriction_filt) %in% c(young_enriched_genes_volcano_filtered, old_enriched_genes_volcano_filtered)] = 'do_plot'
de.markers.glia_no_clonal_restriction_filt$genenames_to_plot[de.markers.glia_no_clonal_restriction_filt$plot_genenames %in% c('do_not_plot')] = NA
ggplot(data=de.markers.glia_no_clonal_restriction_filt, aes(x=avg_log2FC, y=-log10(p_val_adj), label=genenames_to_plot))+
  geom_point()+
  theme_minimal()+
  geom_label_repel(
    box.padding = 0.8, alpha=0.9
  )+
  ggtitle("Differential expression between <=GW20 glia (left) and >GW20 glia (right) -- tRG, oRG, and RG")+
  theme(plot.title = element_text(hjust = 0.5))
image = ggplot(data=de.markers.glia_no_clonal_restriction_filt, aes(x=avg_log2FC, y=-log10(p_val_adj), label=genenames_to_plot))+
  geom_point()+
  theme_minimal()+
  geom_label_repel(
    box.padding = 0.8, alpha=0.9
  )+
  ggtitle("Differential expression between <=GW20 glia (left) and >GW20 glia (right) -- tRG, oRG, and RG")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('volcano_ls_prog_rg_clusts_young_vs_old_filt_manual_genes_v3.svg', image, width=10, height=8)


#### Figure 2 ####
### Fig 2a-b -- upset plots
clone_df <- ls_synthesis@meta.data %>%
  filter(!is.na(Clone_IDs)) %>%
  group_by(Clone_IDs) %>%
  summarise(Clone_size = n(), .groups = "drop")

table(clone_df$Clone_size)
#     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    19    22    23    25 
# 39851  3700  1357   610   300   179    85    59    27    26    23    13     6     4     2     6     1     1     1     1     1
clone_df_v2 = data.frame(clone_df)

## v1 -- using the final lock subcluster identities, tRG and bRG are both RG
compiled_clust_df_v1 = data.frame(Clone_IDs=multicell_clones)
compiled_clust_df_v1$ENs = rep(0,length(multicell_clones))
compiled_clust_df_v1$EX_IPCs = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_local = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_CGE = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_MGE = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_OB = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_IPCs = rep(0,length(multicell_clones))
compiled_clust_df_v1$Astrocytes = rep(0,length(multicell_clones))
compiled_clust_df_v1$OPCs = rep(0,length(multicell_clones))
compiled_clust_df_v1$Oligodendrocytes = rep(0,length(multicell_clones))
compiled_clust_df_v1$RG = rep(0,length(multicell_clones))
compiled_clust_df_v1$tRG = rep(0,length(multicell_clones))
compiled_clust_df_v1$bRG = rep(0,length(multicell_clones))
compiled_clust_df_v1$ENs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$EX_IPCs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_local_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_CGE_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_MGE_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_OB_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$IN_IPCs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$Astrocytes_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$Oligodendrocytes_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$OPCS_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$RG_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$tRG_pct = rep(0,length(multicell_clones))
compiled_clust_df_v1$bRG_pct = rep(0,length(multicell_clones))

for(clone in multicell_clones) {
  clone_size = clone_df$Clone_size[clone_df$Clone_IDs %in% clone][1]
  clone_clusts = table(ls_synthesis@meta.data$subcluster_identity[ls_synthesis@meta.data$Clone_IDs %in% clone])
  clone_row = which(compiled_clust_df_v1$Clone_IDs %in% clone)
  compiled_clust_df_v1[clone_row,"ENs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "ENs"])
  compiled_clust_df_v1[clone_row,"EX_IPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX_IPCs"])
  compiled_clust_df_v1[clone_row,"IN_local"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_local"])
  compiled_clust_df_v1[clone_row,"IN_CGE"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_CGE"])
  compiled_clust_df_v1[clone_row,"IN_MGE"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_MGE"])
  compiled_clust_df_v1[clone_row,"IN_OB"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_OB"])
  compiled_clust_df_v1[clone_row,"IN_IPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_IPCs"])
  compiled_clust_df_v1[clone_row,"Astrocytes"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Astrocytes"])
  compiled_clust_df_v1[clone_row,"Oligodendrocytes"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Oligodendrocytes"])
  compiled_clust_df_v1[clone_row,"OPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])
  compiled_clust_df_v1[clone_row,"RG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "RG"])
  compiled_clust_df_v1[clone_row,"tRG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "tRG"])
  compiled_clust_df_v1[clone_row,"bRG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "bRG"])
  compiled_clust_df_v1[clone_row,"ENs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "ENs"])/clone_size
  compiled_clust_df_v1[clone_row,"EX_IPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX_IPCs"])/clone_size
  compiled_clust_df_v1[clone_row,"IN_local_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_local"])/clone_size
  compiled_clust_df_v1[clone_row,"IN_CGE_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_CGE"])/clone_size
  compiled_clust_df_v1[clone_row,"IN_MGE_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_MGE"])/clone_size
  compiled_clust_df_v1[clone_row,"IN_OB_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_OB"])/clone_size
  compiled_clust_df_v1[clone_row,"IN_IPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_IPCs"])/clone_size
  compiled_clust_df_v1[clone_row,"Astrocytes_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Astrocytes"])/clone_size
  compiled_clust_df_v1[clone_row,"Oligodendrocytes_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Oligodendrocytes"])/clone_size
  compiled_clust_df_v1[clone_row,"OPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])/clone_size
  compiled_clust_df_v1[clone_row,"RG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "RG"])/clone_size
  compiled_clust_df_v1[clone_row,"tRG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "tRG"])/clone_size
  compiled_clust_df_v1[clone_row,"bRG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "bRG"])/clone_size
}
compiled_clust_df_v1
head(compiled_clust_df_v1)
dim(compiled_clust_df_v1)
# [1] 6402   28
dim(compiled_clust_df_v1[complete.cases(compiled_clust_df_v1),])
# [1] 6402   28
upset_notations = c()
for(r in 1:nrow(compiled_clust_df_v1)) {
  upset_notations=c(upset_notations, paste(colnames(compiled_clust_df_v1[r,2:14] %>% select(where(~ any(. != 0)))), collapse="&"))
}
upset_notations
compiled_clust_df_v1$upset = upset_notations
table(compiled_clust_df_v1$upset)
dim(left_join(compiled_clust_df_v1, ls_synthesis@meta.data[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first"))
# [1] 6402   33
compiled_clust_df_v1 = left_join(compiled_clust_df_v1, ls_synthesis@meta.data[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first")
head(compiled_clust_df_v1)

# identify clone types only in one sample
all_clone_types_v1 = unique(compiled_clust_df_v1$upset)
all_clone_types_sample_counts_v1 = data.frame("clone_types" = all_clone_types_v1, "sample_counts" = 0)
all_clone_types_sample_counts_v1
for (ct in all_clone_types_v1){
  all_clone_types_sample_counts_v1$sample_counts[all_clone_types_sample_counts_v1$clone_types %in% ct] = length(table(compiled_clust_df_v1$sample_id[compiled_clust_df_v1$upset %in% ct]))
}
table(all_clone_types_sample_counts_v1$sample_counts)
#  1  2  3  4  5  6  7  8 
# 29 15 11 12  6  9 11  6
single_sample_clone_types_v1 = all_clone_types_sample_counts_v1$clone_types[all_clone_types_sample_counts_v1$sample_counts <= 1]
multi_sample_clone_types_v1 = all_clone_types_sample_counts_v1$clone_types[all_clone_types_sample_counts_v1$sample_counts > 1]
single_sample_clone_type_clone_ids_v1 = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$upset %in% single_sample_clone_types_v1]
multi_sample_clones_v1 = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1]

upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))

dim(left_join(compiled_clust_df_v1, ls_synthesis@meta.data[,c('age_binned' ,'Clone_IDs')], multiple="first"))
compiled_clust_df_v1 = left_join(compiled_clust_df_v1, ls_synthesis@meta.data[,c('age_binned' ,'Clone_IDs')], multiple="first")

pdf('C:/Users/mattg/Desktop/250429_FINAL_FIGURES/Fig2/upset_plot_all_cells_young.pdf', width=12, height=8)
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$age_binned == "<=GW20"])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()
pdf('C:/Users/mattg/Desktop/250429_FINAL_FIGURES/Fig2/upset_plot_all_cells_old.pdf', width=12, height=8)
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$age_binned == ">GW20"])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

### Fig 2a and Fig 2b -- chord diagrams
# install.packages("circlize")
library(circlize)
subclust_names = rownames(table(ls_synthesis@meta.data$subcluster_identity))
subclust_names
#[1] "RG"               "Astrocytes"       "OPCs"             "Oligodendrocytes" "EX_IPCs"         
#[6] "ENs"              "IN_IPCs"          "IN_local"         "IN_CGE"           "IN_OB"           
#[11] "IN_MGE"
length(subclust_names)
#[1] 11
chord_matrix = matrix(rep(0, length(subclust_names)*length(subclust_names)), ncol=length(subclust_names))
dim(chord_matrix)
#[1] 11 11
rownames(chord_matrix) = c("RG","Astrocytes","OPCs","Oligodendrocytes","EX_IPCs","ENs","IN_IPCs","IN_local","IN_CGE","IN_OB","IN_MGE")
colnames(chord_matrix) = c("RG","Astrocytes","OPCs","Oligodendrocytes","EX_IPCs","ENs","IN_IPCs","IN_local","IN_CGE","IN_OB","IN_MGE")
multicell_clones = clone_df$Clone_IDs[clone_df$Clone_size>=2]
sum(!is.na(multicell_clones))
# [1] 6402
multicell_clones = multicell_clones[!is.na(multicell_clones)]
ls_synthesis_meta_multicell = ls_synthesis@meta.data[ls_synthesis@meta.data$Clone_IDs %in% multicell_clones,]
dim(ls_synthesis_meta_multicell)
# [1] 18830    14
colnames(compiled_clust_df_v1)[2:12]
# [1] "ENs"              "EX_IPCs"          "IN_local"         "IN_CGE"           "IN_MGE"           "IN_OB"            "IN_IPCs"         
# [8] "Astrocytes"       "OPCs"             "Oligodendrocytes" "RG

## generate the chord diagram for <=GW20 samples
for (clust_id in c(subclust_names)) {
  clone_list_temp = ls_synthesis_meta_multicell$Clone_IDs[ls_synthesis_meta_multicell$subcluster_identity %in% clust_id & ls_synthesis_meta_multicell$age_binned == "<=GW20"]
  clone_list_temp = clone_list_temp[!is.na(clone_list_temp)]
  for (clust_name in subclust_names){
    if (clust_id != clust_name){
      chord_matrix[rownames(chord_matrix) %in% clust_id, colnames(chord_matrix) %in% clust_name] = sum(compiled_clust_df_v1[compiled_clust_df_v1$Clone_IDs %in% clone_list_temp, colnames(compiled_clust_df_v1) %in% clust_name] >=1)
    } else{
      chord_matrix[rownames(chord_matrix) %in% clust_id, colnames(chord_matrix) %in% clust_name] = sum(compiled_clust_df_v1[compiled_clust_df_v1$Clone_IDs %in% clone_list_temp, colnames(compiled_clust_df_v1) %in% clust_name] >1) #for self-interactions, only count it if there are multiple of the given cell type
    }
  }
}
chord_matrix_upper = upper.tri(chord_matrix, diag=T) * chord_matrix
chord_matrix_upper
celltype_cols_chord = c("#F8766D", "#DB8E00", "#AEA200", "#64B200", "#00BD5C", "#00C1A7", "#00BADE", "#00A6FF", "#B385FF", "#EF67EB", "#FF63B6")
chordDiagram(chord_matrix_upper, grid.col = celltype_cols_chord)
pdf('C:/Users/mattg/Desktop/250429_FINAL_FIGURES/Fig2/chord_diagram_young.pdf', width=8, height=8)
chordDiagram(chord_matrix_upper, grid.col = celltype_cols_chord)
dev.off()

## generate the chord diagram for >GW20 samples
for (clust_id in c(subclust_names)) {
  clone_list_temp = ls_synthesis_meta_multicell$Clone_IDs[ls_synthesis_meta_multicell$subcluster_identity %in% clust_id & ls_synthesis_meta_multicell$age_binned == ">GW20"]
  clone_list_temp = clone_list_temp[!is.na(clone_list_temp)]
  for (clust_name in subclust_names){
    if (clust_id != clust_name){
      chord_matrix[rownames(chord_matrix) %in% clust_id, colnames(chord_matrix) %in% clust_name] = sum(compiled_clust_df_v1[compiled_clust_df_v1$Clone_IDs %in% clone_list_temp, colnames(compiled_clust_df_v1) %in% clust_name] >=1)
    } else{
      chord_matrix[rownames(chord_matrix) %in% clust_id, colnames(chord_matrix) %in% clust_name] = sum(compiled_clust_df_v1[compiled_clust_df_v1$Clone_IDs %in% clone_list_temp, colnames(compiled_clust_df_v1) %in% clust_name] >1) #for self-interactions, only count it if there are multiple of the given cell type
    }
  }
}
chord_matrix_upper = upper.tri(chord_matrix, diag=T) * chord_matrix
chord_matrix_upper
celltype_cols_chord = c("#F8766D", "#DB8E00", "#AEA200", "#64B200", "#00BD5C", "#00C1A7", "#00BADE", "#00A6FF", "#B385FF", "#EF67EB", "#FF63B6")
chordDiagram(chord_matrix_upper, grid.col = celltype_cols_chord)
dev.off()

### Fig 2c -- UMAP of ls_in with subcluster identities
ls_synthesis_in@meta.data$subcluster_identity = factor(ls_synthesis_in@meta.data$subcluster_identity, levels=c(
  "IN_local",
  "IN_IPCs",
  "IN_MGE",
  "IN_CGE",
  "IN_OB"
))
ls_synthesis_in = SetIdent(ls_synthesis_in, value = "subcluster_identity")

DimPlot(ls_synthesis_in, pt.size=0.5, group.by="subcluster_identity", raster=F)
DimPlot(ls_synthesis_in, pt.size=2, group.by="subcluster_identity ", raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
image = DimPlot(ls_synthesis_in, pt.size=2, group.by="subcluster_identity", raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
ggsave(file="in_umap_subcluster_idents.png", plot=image, width=12, height=8)
ggsave(file="in_umap_subcluster_idents.svg", plot=image, width=12, height=8)

### Fig 2d -- UMAP of ls_in with pseudoage
DimPlot(ls_synthesis_in, pt.size=3, group.by="age_binned", cols = c(young_col, old_col), raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
image = DimPlot(ls_synthesis_in, pt.size=3, group.by="age_binned", cols = c(young_col, old_col), raster=F, shuffle=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=40), plot.title = element_blank())
ggsave(file="in_umap_pseudoage_idents_revised.png", plot=image, width=12, height=8)
ggsave(file="in_umap_pseudoage_idents_revised.svg", plot=image, width=12, height=8)

### Fig 2e -- UMAP with multicellular clones
DimPlot(ls_synthesis_in, sizes.highlight=2, cells.highlight = rownames(ls_synthesis_in@meta.data)[ls_synthesis_in@meta.data$Clone_size >=2], cols.highlight = "#c44e33") +
  ggtitle('EX and IN shared clones (Clone size >= 2)') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.title = element_text(size=25, hjust=0.5))
image = DimPlot(ls_synthesis_in, sizes.highlight=2, cells.highlight = rownames(ls_synthesis_in@meta.data)[ls_synthesis_in@meta.data$Clone_size >=2], cols.highlight = "#c44e33") +
  ggtitle('EX and IN shared clones (Clone size >= 2)') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.title = element_text(size=25, hjust=0.5))
ggsave(file="in_umap_multicell_clones.png", plot=image, width=12, height=8)
ggsave(file="in_umap_multicell_clones.svg", plot=image, width=12, height=8)

### Fig 2f -- UMAP with shared glutamatergic and inhibitory clones
ex_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$ENs>=1 | compiled_clust_df_v1$EX_IPCs>=1]
length(ex_clones)
# [1] 3483
in_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$IN_local>=1 | compiled_clust_df_v1$IN_CGE>=1 | compiled_clust_df_v1$IN_MGE>=1 | compiled_clust_df_v1$IN_OB>=1 | compiled_clust_df_v1$IN_IPCs>=1]
length(in_clones)
# [1] 2128
ex_in_shared_clones = ex_clones[ex_clones %in% in_clones]
length(ex_in_shared_clones)

DimPlot(ls_synthesis_in, sizes.highlight=2, cells.highlight = rownames(ls_synthesis_in@meta.data)[ls_synthesis_in@meta.data$Clone_IDs %in% ex_in_shared_clones & ls_synthesis_in@meta.data$Clone_IDs %in% multi_sample_clones_v1], cols.highlight = "#f27527") +
  ggtitle('EX and IN shared clones (Clone size >= 2)') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.title = element_text(size=25, hjust=0.5))
image = DimPlot(ls_synthesis_in, sizes.highlight=2, cells.highlight = rownames(ls_synthesis_in@meta.data)[ls_synthesis_in@meta.data$Clone_IDs %in% ex_in_shared_clones & ls_synthesis_in@meta.data$Clone_IDs %in% multi_sample_clones_v1], cols.highlight = "#f27527") +
  ggtitle('EX and IN shared clones (Clone size >= 2)') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.title = element_text(size=25, hjust=0.5))
ggsave(file="in_umap_glutamatergic_shared_clones.png", plot=image, width=12, height=8)
ggsave(file="in_umap_glutamatergic_shared_clones.svg", plot=image, width=12, height=8)


### Fig 2g -- dotplot for ls_in
ls_synthesis_in = SetIdent(ls_synthesis_in, value = "subcluster_identity")
ls_synthesis_in@active.ident = factor(ls_synthesis_in@active.ident, levels=c(
  "IN_IPCs",
  "IN_local",
  "IN_CGE",
  "IN_OB",
  "IN_MGE"
))
DotPlot(ls_synthesis_in, features=c("GAD2", "DLX1", "DLX2", "DLX5", "FOXG1", "ERBB4",
                                    "NR2F2", "ADARB2", "PROX1", "CALB2", "VIP", "RELN", 
                                    "LHX6", "MAF", "SST",
                                    "MEIS2", "SP8",
                                    "PBX3", "TSHZ1",
                                    "SOX6", "NR2F1", "SCGN", "PAX6", "ST18", "NPY", "MKI67")) +
  ggtitle('Interneuron marker gene expression')+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=18, angle=30, hjust=1, vjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1)
  )
image = DotPlot(ls_synthesis_in, features=c("GAD2", "DLX1", "DLX2", "DLX5", "FOXG1", "ERBB4",
                                            "NR2F2", "ADARB2", "PROX1", "CALB2", "VIP", "RELN", 
                                            "LHX6", "MAF", "SST",
                                            "MEIS2", "SP8",
                                            "PBX3", "TSHZ1",
                                            "SOX6", "NR2F1", "SCGN", "PAX6", "ST18", "NPY", "MKI67")) +
  ggtitle('Interneuron marker gene expression')+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=18, angle=30, hjust=1, vjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1)
  )
ggsave(file="in_marker_dotplot_with_st18.svg", plot=image, width=14, height=8)


### Fig 2j -- scatterplot with trendline from stains, triple positives/DLX2+SCGN+ barchart (REVISION)
scgn_avg_counts_revision = read.csv('/media/chang/HDD-4/mgkeefe/R_work_in_progress_cajal/local_STICR_synthesis_figs/240727_revision_figs/240728_local_STICR_scgn_average_counts.csv')
scgn_raw_counts_revision = read.csv('/media/chang/HDD-4/mgkeefe/R_work_in_progress_cajal/local_STICR_synthesis_figs/240727_revision_figs/240728_local_STICR_scgn_raw_counts.csv')
colnames(scgn_raw_counts_revision)[1] = "image"
colnames(scgn_avg_counts_revision)[1] = "image"

scgn_raw_counts_revision$age_pseudo = NA
scgn_raw_counts_revision$age_pseudo[scgn_raw_counts_revision$age <= 20] = "<=GW20"
scgn_raw_counts_revision$age_pseudo[scgn_raw_counts_revision$age > 20] = ">GW20"
scgn_avg_counts_revision$age_pseudo = NA
scgn_avg_counts_revision$age_pseudo[scgn_avg_counts_revision$age <= 20] = "<=GW20"
scgn_avg_counts_revision$age_pseudo[scgn_avg_counts_revision$age > 20] = ">GW20"


ggplot() +
  #geom_point(data = scgn_avg_counts_revision, aes(x=age, y=DLX2_SCGN_of_DLX2), color="darkorange") +
  geom_jitter(data = scgn_raw_counts_revision, aes(x=age, y=DLX2_SCGN_of_DLX2, size=2), color="#ffc561", width=0.1) +
  geom_smooth(data = scgn_avg_counts_revision, aes(x=age, y=DLX2_SCGN_of_DLX2), color="darkorange", method = "lm") +
  #geom_point(data = scgn_avg_counts_revision, aes(x=age, y=triple_of_DLX2_SCGN), color="#1a6d78") +
  #geom_jitter(data = scgn_raw_counts_revision, aes(x=age, y=triple_of_DLX2, size=2), color="cyan", width=0.1) +
  #geom_smooth(data = scgn_avg_counts_revision, aes(x=age, y=triple_of_DLX2), color="#1a6d78", method = "lm") +
  ylab('Ratio of cells') +
  xlab('Sample age') +
  ggtitle('Emergence of DLX2+/SCGN+ cells') +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5))

ggplot(data = scgn_avg_counts_revision, aes(x=age_pseudo, y=triple_of_DLX2_SCGN, fill=age_pseudo, color=age_pseudo)) +
  geom_bar(stat="summary", alpha=0.6) +
  geom_jitter(width=0.1, height=0, size=3, stroke=1, colour='black', fill='gray') +
  scale_fill_manual(values=c("#c9c9c9", "#999999")) +
  scale_color_manual(values=c("#c9c9c9", "#999999")) +
  theme_bw() +
  ggtitle("Emergence of DLX2+/SCGN+/PAX6+ cells") +
  labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Sample ages (GW)") +
  ylab("Ratio of DLX2+SCGN+ cells") +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
image = ggplot(data = scgn_avg_counts_revision, aes(x=age_pseudo, y=triple_of_DLX2_SCGN, fill=age_pseudo, color=age_pseudo)) +
  geom_bar(stat="summary", alpha=0.6) +
  geom_jitter(width=0.1, height=0, stroke=1) +
  scale_fill_manual(values=c("#c9c9c9", "#999999")) +
  scale_color_manual(values=c("#c9c9c9", "#999999")) +
  theme_bw() +
  ggtitle("Emergence of DLX2+/SCGN+/PAX6+ cells") +
  labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Sample ages (GW)") +
  ylab("Ratio of DLX2+SCGN+ cells") +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5))
ggsave('revision_scgn_counts_emergence_of_scgn_pax6_of_dlx2_scgn.svg', image, width=10, height=8)

summary(scgn_avg_counts_revision$triple_of_DLX2_SCGN[scgn_avg_counts_revision$age <=20])
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.01325 0.02649 0.02149 0.03223 0.03797
sd(scgn_avg_counts_revision$triple_of_DLX2_SCGN[scgn_avg_counts_revision$age <=20])/sqrt(length(scgn_avg_counts_revision$triple_of_DLX2_SCGN[scgn_avg_counts_revision$age <=20]))
#[1] 0.011244

summary(scgn_avg_counts_revision$triple_of_DLX2_SCGN[scgn_avg_counts_revision$age >20])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06145 0.07431 0.09274 0.14358 0.20796 0.28142
sd(scgn_avg_counts_revision$triple_of_DLX2_SCGN[scgn_avg_counts_revision$age >20])/sqrt(length(scgn_avg_counts_revision$triple_of_DLX2_SCGN[scgn_avg_counts_revision$age <=20]))
#[1] 0.05573253

means_with_sds_barcharts = data.frame("ages" = c("young", "old"), "means" = c(0.02149, 0.14358), "sds" = c(0.011244, 0.05573253))
means_with_sds_barcharts
ggplot(means_with_sds_barcharts, aes(x=ages, y=means)) +
  geom_bar(stat="identity", alpha=0.6) +
  geom_jitter(width=0.1, height=0, stroke=1) +
  geom_errorbar(aes(x=ages, ymin=means-sds, ymax=means+sds), width=0.4, alpha=0.9, size=1.3)
image = ggplot(means_with_sds_barcharts, aes(x=ages, y=means)) +
  geom_bar(stat="identity", alpha=0.6) +
  geom_jitter(width=0.1, height=0, stroke=1) +
  geom_errorbar(aes(x=ages, ymin=means-sds, ymax=means+sds), width=0.4, alpha=0.9, size=1.3)
ggsave('revision_scgn_pax6_barplot_with_sem_error_bars.svg' , image, width=10, height=8)

install.packages("ggpubr", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(ggpubr)
ggscatter(scgn_raw_counts, x = "age", y = "triple_of_DLX2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson") +
  ylab('Ratio of DLX2+ cells') +
  xlab('Sample age') +
  ggtitle('Emergence of DLX2+/SCGN+/PAX6+ cells') +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=20), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5, ))
image = ggscatter(scgn_raw_counts_revision, x = "age", y = "triple_of_DLX2", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson") +
  ylab('Ratio of DLX2+ cells') +
  xlab('Sample age') +
  ggtitle('Emergence of DLX2+/SCGN+/PAX6+ cells') +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=20), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5, ))
ggsave('revision_scgn_counts_emergence_of_dlx2_scgn_pax6_only_lineplot_with_corr.svg', image, width=10, height=8)


#### Figure 3 ####
### Fig 3b -- multicellular clones from >GW20 samples
osvz_col = "#FF6666"
vz_col = "#00C3C6"
DimPlot(ls_synthesis, cells.highlight = 
          rownames(ls_synthesis@meta.data)[
            ls_synthesis@meta.data$Clone_IDs %in% multicell_clones &
              ls_synthesis@meta.data$age_binned %in% ">GW20" &
              ls_synthesis@meta.data$GZ_origin %in% "VZ"])

DimPlot(ls_synthesis, cells.highlight = list(
  rownames(ls_synthesis@meta.data)[
  ls_synthesis@meta.data$Clone_IDs %in% multicell_clones &
    ls_synthesis@meta.data$age_binned %in% ">GW20" &
    ls_synthesis@meta.data$GZ_origin %in% "VZ"],
  rownames(ls_synthesis@meta.data)[
    ls_synthesis@meta.data$Clone_IDs %in% multicell_clones &
      ls_synthesis@meta.data$age_binned %in% ">GW20" &
      ls_synthesis@meta.data$GZ_origin %in% "OSVZ"]),
        cols.highlight = c(osvz_col, vz_col), na.value = "lightgray", sizes.highlight = 3, ) +
  ggtitle("") +
  labs(color="GZ labeled")+
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20), legend.title = element_text(size=25))
image = 
  DimPlot(ls_synthesis, cells.highlight = list(
    rownames(ls_synthesis@meta.data)[
      ls_synthesis@meta.data$Clone_IDs %in% multicell_clones &
        ls_synthesis@meta.data$age_binned %in% ">GW20" &
        ls_synthesis@meta.data$GZ_origin %in% "VZ"],
    rownames(ls_synthesis@meta.data)[
      ls_synthesis@meta.data$Clone_IDs %in% multicell_clones &
        ls_synthesis@meta.data$age_binned %in% ">GW20" &
        ls_synthesis@meta.data$GZ_origin %in% "OSVZ"]),
    cols.highlight = c(osvz_col, vz_col), na.value = "lightgray", sizes.highlight = 3, ) +
  ggtitle("") +
  labs(color="GZ labeled")+
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20), legend.title = element_text(size=25))
ggsave("full_umap_>GW20_clones_by_gz.png", plot = image, width=10, height=8)
ggsave("full_umap_>GW20_clones_by_gz.svg", plot = image, width=10, height=8)


### Fig 3c -- EX gene expression in full UMAP
FeaturePlot(ls_synthesis, features=c("SATB2"), order=T, pt.size=0.8, min.cutoff="q40", cols=c("#e6e8e7","#08519C")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("SATB2"), order=T, pt.size=0.4, min.cutoff="q40", cols=c("#e6e8e7","#08519C")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("revision_satb2_full_umap.png", plot=image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("TBR1"), order=T, pt.size=0.4, min.cutoff="q40", cols=c("#e6e8e7","#08519C")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("TBR1"), order=T, pt.size=0.4, min.cutoff="q40", cols=c("#e6e8e7","#08519C")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("revision_tbr1_full_umap.png", plot=image, width=4, height=3.2)
image = FeaturePlot(ls_synthesis, features=c("NEUROD2"), order=T, pt.size=0.4, min.cutoff="q40", cols=c("#e6e8e7","#08519C")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("revision_neurod2_full_umap.png", plot=image, width=4, height=3.2)
image = FeaturePlot(ls_synthesis, features=c("SLC17A7"), order=T, pt.size=0.4, min.cutoff="q40", cols=c("#e6e8e7", "#08519C")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("revision_slc17a7_full_umap.png", plot=image, width=4, height=3.2)

### Fig 3e -- tRG clones highlighted on full UMAP
brg_col = '#ff8587'
trg_col = '#00c4c7'
brg_clones_revision = unique(ls_synthesis_prog_glia@meta.data$Clone_IDs[ls_synthesis_prog_glia@meta.data$subcluster_identity %in% 'bRG'])
trg_clones_revision = unique(ls_synthesis_prog_glia@meta.data$Clone_IDs[ls_synthesis_prog_glia@meta.data$subcluster_identity %in% 'tRG'])
length(brg_clones_revision)
#[1] 1716
length(trg_clones_revision)
#[1] 335
trg_clones_revision = trg_clones_revision[complete.cases(trg_clones_revision)]
brg_clones_revision = brg_clones_revision[complete.cases(brg_clones_revision)]
DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$Clone_IDs %in% trg_clones_revision],
        cols.highlight=c(trg_col), sizes.highlight=3, cols = "lightgrey") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
image = DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$Clone_IDs %in% trg_clones_revision],
                cols.highlight=c(trg_col), sizes.highlight=3, cols = "lightgrey") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
ggsave('full_umap_trg_clones_highlighted.svg', image, width=10, height=8)
ggsave('full_umap_trg_clones_highlighted8.png', image, width=10, height=8)


### Fig 3e -- tRG clones highlighted on full UMAP
brg_col = '#ff8587'
trg_col = '#00c4c7'
brg_clones_revision = unique(ls_synthesis_prog_glia@meta.data$Clone_IDs[ls_synthesis_prog_glia@meta.data$subcluster_identity %in% 'bRG'])
trg_clones_revision = unique(ls_synthesis_prog_glia@meta.data$Clone_IDs[ls_synthesis_prog_glia@meta.data$subcluster_identity %in% 'tRG'])
length(brg_clones_revision)
#[1] 1716
length(trg_clones_revision)
#[1] 335
trg_clones_revision = trg_clones_revision[complete.cases(trg_clones_revision)]
brg_clones_revision = brg_clones_revision[complete.cases(brg_clones_revision)]
DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$Clone_IDs %in% trg_clones_revision],
        cols.highlight=c(trg_col), sizes.highlight=3, cols = "lightgrey") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
image = DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$Clone_IDs %in% trg_clones_revision],
                cols.highlight=c(trg_col), sizes.highlight=3, cols = "lightgrey") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
ggsave('full_umap_trg_clones_highlighted.svg', image, width=10, height=8)
ggsave('full_umap_trg_clones_highlighted8.png', image, width=10, height=8)

### Fig 3f -- tRG upset plot
## v2 -- bRG as RG
ls_synthesis_meta = ls_synthesis@meta.data
ls_synthesis_meta$subcluster_identity_rg = ls_synthesis_meta$subcluster_identity
ls_synthesis_meta$subcluster_identity_rg[rownames(ls_synthesis_meta) %in% rownames(ls_synthesis_prog_glia@meta.data)[ls_synthesis_prog_glia@meta.data$subcluster_identity  %in% "tRG"]] = "tRG"
ls_synthesis_meta$subcluster_identity_rg[rownames(ls_synthesis_meta) %in% rownames(ls_synthesis_prog_glia@meta.data)[ls_synthesis_prog_glia@meta.data$subcluster_identity  %in% "bRG"]] = "bRG"
table(ls_synthesis_meta$subcluster_identity_rg)

compiled_clust_df_v2 = data.frame(Clone_IDs=multicell_clones)
compiled_clust_df_v2$ENs = rep(0,length(multicell_clones))
compiled_clust_df_v2$EX_IPCs = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_local = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_CGE = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_MGE = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_OB = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_IPCs = rep(0,length(multicell_clones))
compiled_clust_df_v2$Astrocytes = rep(0,length(multicell_clones))
compiled_clust_df_v2$OPCs = rep(0,length(multicell_clones))
compiled_clust_df_v2$Oligodendrocytes = rep(0,length(multicell_clones))
compiled_clust_df_v2$RG = rep(0,length(multicell_clones))
compiled_clust_df_v2$tRG = rep(0,length(multicell_clones))
compiled_clust_df_v2$bRG = rep(0,length(multicell_clones))
compiled_clust_df_v2$ENs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$EX_IPCs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_local_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_CGE_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_MGE_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_OB_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$IN_IPCs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$Astrocytes_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$Oligodendrocytes_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$OPCS_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$RG_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$tRG_pct = rep(0,length(multicell_clones))
compiled_clust_df_v2$bRG_pct = rep(0,length(multicell_clones))

for(clone in multicell_clones) {
  clone_size = clone_df$Clone_size[clone_df$Clone_IDs %in% clone][1]
  clone_clusts = table(ls_synthesis_meta$subcluster_identity_rg[ls_synthesis_meta$Clone_IDs %in% clone])
  clone_row = which(compiled_clust_df_v2$Clone_IDs %in% clone)
  compiled_clust_df_v2[clone_row,"ENs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "ENs"])
  compiled_clust_df_v2[clone_row,"EX_IPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX_IPCs"])
  compiled_clust_df_v2[clone_row,"IN_local"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_local"])
  compiled_clust_df_v2[clone_row,"IN_CGE"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_CGE"])
  compiled_clust_df_v2[clone_row,"IN_MGE"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_MGE"])
  compiled_clust_df_v2[clone_row,"IN_OB"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_OB"])
  compiled_clust_df_v2[clone_row,"IN_IPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_IPCs"])
  compiled_clust_df_v2[clone_row,"Astrocytes"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Astrocytes"])
  compiled_clust_df_v2[clone_row,"Oligodendrocytes"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Oligodendrocytes"])
  compiled_clust_df_v2[clone_row,"OPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])
  compiled_clust_df_v2[clone_row,"RG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% c("bRG", "RG")])
  compiled_clust_df_v2[clone_row,"tRG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "tRG"])
  #compiled_clust_df_v2[clone_row,"bRG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "bRG"])
  compiled_clust_df_v2[clone_row,"ENs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "ENs"])/clone_size
  compiled_clust_df_v2[clone_row,"EX_IPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX_IPCs"])/clone_size
  compiled_clust_df_v2[clone_row,"IN_local_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_local"])/clone_size
  compiled_clust_df_v2[clone_row,"IN_CGE_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_CGE"])/clone_size
  compiled_clust_df_v2[clone_row,"IN_MGE_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_MGE"])/clone_size
  compiled_clust_df_v2[clone_row,"IN_OB_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_OB"])/clone_size
  compiled_clust_df_v2[clone_row,"IN_IPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_IPCs"])/clone_size
  compiled_clust_df_v2[clone_row,"Astrocytes_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Astrocytes"])/clone_size
  compiled_clust_df_v2[clone_row,"Oligodendrocytes_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Oligodendrocytes"])/clone_size
  compiled_clust_df_v2[clone_row,"OPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])/clone_size
  compiled_clust_df_v2[clone_row,"RG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% c("bRG", "RG")])/clone_size
  compiled_clust_df_v2[clone_row,"tRG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "tRG"])/clone_size
  #compiled_clust_df_v2[clone_row,"bRG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "bRG"])/clone_size
}
compiled_clust_df_v2
head(compiled_clust_df_v2)
dim(compiled_clust_df_v2)
#[1] 6402   28
compiled_clust_df_v2 = compiled_clust_df_v2[complete.cases(compiled_clust_df_v2),] #drop NA row

upset_notations = c()
for(r in 1:nrow(compiled_clust_df_v2)) {
  upset_notations=c(upset_notations, paste(colnames(compiled_clust_df_v2[r,2:14] %>% select(where(~ any(. != 0)))), collapse="&"))
}
upset_notations
compiled_clust_df_v2$upset = upset_notations
table(compiled_clust_df_v2$upset)
dim(left_join(compiled_clust_df_v2, ls_synthesis_meta[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first"))
#[1] 6402   33
compiled_clust_df_v2 = left_join(compiled_clust_df_v2, ls_synthesis_meta[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first")
head(compiled_clust_df_v2)

# drop clone types only in one sample
all_clone_types_v2 = unique(compiled_clust_df_v2$upset)
all_clone_types_sample_counts_v2 = data.frame("clone_types" = all_clone_types_v2, "sample_counts" = 0)
all_clone_types_sample_counts_v2
for (ct in all_clone_types_v2){
  all_clone_types_sample_counts_v2$sample_counts[all_clone_types_sample_counts_v2$clone_types %in% ct] = length(table(compiled_clust_df_v2$sample_id[compiled_clust_df_v2$upset %in% ct]))
}
table(all_clone_types_sample_counts_v2$sample_counts)
# 1  2  3  4  5  6  7  8 
#45 21 14 13  7 11 11  6
length(all_clone_types_sample_counts_v2$clone_types[all_clone_types_sample_counts_v2$sample_counts <= 1])
#[1] 45
length(all_clone_types_sample_counts_v2$clone_types[all_clone_types_sample_counts_v2$sample_counts > 1])
#[1] 83

## instead of using those, going to use the multi sample correction from when RG are grouped for consistency
multi_sample_clone_types_v2 = unique(compiled_clust_df_v2$upset[compiled_clust_df_v2$Clone_IDs %in% multi_sample_clones_v1])

#install.packages("UpSetR")
#library(UpSetR)

upset(fromExpression(table(compiled_clust_df_v2$upset[compiled_clust_df_v2$upset %in% multi_sample_clone_types_v2 & compiled_clust_df_v2$Clone_IDs %in% trg_clones_revision])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dim(left_join(compiled_clust_df_v2, ls_synthesis_meta[,c('age_binned' ,'Clone_IDs')], multiple="first"))
compiled_clust_df_v2 = left_join(compiled_clust_df_v2, ls_synthesis_meta[,c('age_binned' ,'Clone_IDs')], multiple="first")

pdf('C:/Users/mattg/Desktop/250510_FINAL_FIGURES/Fig3_upset_tRG_clones_bRG_as_RG.pdf', width=12, height=8)
upset(fromExpression(table(compiled_clust_df_v2$upset[compiled_clust_df_v2$upset %in% multi_sample_clone_types_v2 & compiled_clust_df_v2$Clone_IDs %in% trg_clones_revision])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

#### Figure 4
### Fig 4d,f -- RNAscope quantifications (REVISION)
#rnascope_ratios = read.csv('/media/chang/HDD-4/mgkeefe/R_work_in_progress_cajal/local_STICR_synthesis_figs/240727_revision_figs/240728_rnascope_pivot_ratios.csv')
table(rnascope_ratios$stain)
rnascope_ratios$stain = factor(rnascope_ratios$stain, levels=c("TBR1_of_DAPI", "CPLX3_TBR1_of_TBR1", "NURR1_TBR1_of_TBR1", "CTGF_TBR1_of_TBR1"))
rnascope_ratios

rnascope_ratios_tbr1_dapi = rnascope_ratios[rnascope_ratios$stain %in% 'TBR1_of_DAPI',]
rnascope_ratios_tbr1_dapi

rnascope_ratios_not_tbr1_dapi = rnascope_ratios[rnascope_ratios$stain %in% c("CPLX3_TBR1_of_TBR1", "NURR1_TBR1_of_TBR1", "CTGF_TBR1_of_TBR1"),]
rnascope_ratios_not_tbr1_dapi

ggplot(data = rnascope_ratios_tbr1_dapi, aes(x=stain, y=Ratio, fill=region)) +
  geom_bar(stat="summary", alpha=0.6, position='dodge') +
  geom_point(size=3, stroke=1, colour='black', position=position_jitterdodge(dodge.width=0.75, jitter.width = 0.1)) +
  scale_fill_manual(values=c("#c9c9c9", "#999999")) +
  scale_color_manual(values=c("#c9c9c9", "#999999")) +
  theme_bw() +
  ggtitle("Co-expressed markers in TBR1+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
image = ggplot(data = rnascope_ratios_tbr1_dapi, aes(x=stain, y=Ratio, fill=region)) +
  geom_bar(stat="summary", alpha=0.6, position='dodge') +
  geom_point(size=3, stroke=1, colour='black', position=position_jitterdodge(dodge.width=0.75, jitter.width = 0.1)) +
  scale_fill_manual(values=c("#c9c9c9", "#999999")) +
  scale_color_manual(values=c("#c9c9c9", "#999999")) +
  theme_bw() +
  ggtitle("Co-expressed markers in DAPI+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
ggsave('rnascope_tbr1_of_dapi_barchart.svg', image, width=10, height=8, device = 'svg', dpi = 300)

ggplot(data = rnascope_ratios_not_tbr1_dapi, aes(x=stain, y=Ratio, fill=region)) +
  geom_bar(stat="summary", alpha=0.6, position='dodge') +
  geom_point(size=3, stroke=1, colour='black', position=position_jitterdodge(dodge.width=0.75, jitter.width = 0.1)) +
  scale_fill_manual(values=c("#c9c9c9", "#999999")) +
  scale_color_manual(values=c("#c9c9c9", "#999999")) +
  theme_bw() +
  ggtitle("Co-expressed markers in TBR1+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
image = ggplot(data = rnascope_ratios_not_tbr1_dapi, aes(x=stain, y=Ratio, fill=region)) +
  geom_bar(stat="summary", alpha=0.6, position='dodge') +
  geom_point(size=3, stroke=1, colour='black', position=position_jitterdodge(dodge.width=0.75, jitter.width = 0.1)) +
  scale_fill_manual(values=c("#c9c9c9", "#999999")) +
  scale_color_manual(values=c("#c9c9c9", "#999999")) +
  theme_bw() +
  ggtitle("Co-expressed markers in TBR1+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
ggsave('rnascope_sp_markers_of_tbr1_barchart.svg', image, width=10, height=8, device='svg', dpi=300)

## SP RNAscope plot: make separate graphs for SP and VZ
ggplot(data = subset(rnascope_ratios_not_tbr1_dapi, rnascope_ratios_not_tbr1_dapi$region %in% 'SP'), aes(x=stain, y=Ratio)) +
  geom_bar(stat="summary", alpha=0.6, fill="#c9c9c9", colour="#999999") +
  geom_point(size=3, stroke=1, colour='black', position=position_jitter(width = 0.1)) +
  scale_fill_manual(values=c("#c9c9c9")) +
  scale_color_manual(values=c("#c9c9c9")) +
  theme_bw() +
  ggtitle("Co-expressed markers in SP TBR1+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
image = ggplot(data = subset(rnascope_ratios_not_tbr1_dapi, rnascope_ratios_not_tbr1_dapi$region %in% 'SP'), aes(x=stain, y=Ratio)) +
  geom_bar(stat="summary", alpha=0.6, fill="#c9c9c9", colour="#999999") +
  geom_point(size=3, stroke=1, colour='black', position=position_jitter(width = 0.1)) +
  scale_fill_manual(values=c("#c9c9c9")) +
  scale_color_manual(values=c("#c9c9c9")) +
  theme_bw() +
  ggtitle("Co-expressed markers in SP TBR1+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
ggsave('rnascope_sp_markers_of_tbr1_barchart_SP_only.svg', image, width=10, height=8, device='svg', dpi=300)

ggplot(data = subset(rnascope_ratios_not_tbr1_dapi, rnascope_ratios_not_tbr1_dapi$region %in% 'VZ'), aes(x=stain, y=Ratio)) +
  geom_bar(stat="summary", alpha=0.6, fill="#999999", colour="#666666") +
  geom_point(size=3, stroke=1, colour='black', position=position_jitter(width = 0.1)) +
  scale_fill_manual(values=c("#999999")) +
  scale_color_manual(values=c("#999999")) +
  theme_bw() +
  ggtitle("Co-expressed markers in VZ TBR1+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
image = ggplot(data = subset(rnascope_ratios_not_tbr1_dapi, rnascope_ratios_not_tbr1_dapi$region %in% 'VZ'), aes(x=stain, y=Ratio)) +
  geom_bar(stat="summary", alpha=0.6, fill="#999999", colour="#666666") +
  geom_point(size=3, stroke=1, colour='black', position=position_jitter(width = 0.1)) +
  scale_fill_manual(values=c("#999999")) +
  scale_color_manual(values=c("#999999")) +
  theme_bw() +
  ggtitle("Co-expressed markers in VZ TBR1+ cells") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker coexpression") +
  ylab("Ratio of cells") +
  ylim(c(-0.01,0.61))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
ggsave('rnascope_sp_markers_of_tbr1_barchart_VZ_only.svg', image, width=10, height=8, device='svg', dpi=300)

### Figure4i -- BrdU quantifications (REVISION)
## BrdU plot: change to show all 9 points and color by individual, rather than showing the averages for individuals
setwd('/media/chang/HDD-4/mgkeefe/R_work_in_progress_cajal/local_STICR_synthesis_figs/240806_revision_figs')

#brdu_avg_counts_revision = read.csv('/media/chang/HDD-4/mgkeefe/R_work_in_progress_cajal/local_STICR_synthesis_figs/240727_revision_figs/240728_brdu_quant_avg.csv')
brdu_raw_counts_revision = read.csv('/media/chang/HDD-4/mgkeefe/R_work_in_progress_cajal/local_STICR_synthesis_figs/240806_revision_figs/240806_brdu_quant_raw.csv')

brdu_raw_counts_revision

ggplot(data = brdu_avg_counts_revision, aes(x=stain, y=double_positive_ratio_of_BrdU, fill=stain, color=stain)) +
  geom_bar(stat="summary", alpha=0.6) +
  #geom_jitter(width=0.1, height=0, size=3, stroke=1, colour='black', fill='gray') +
  geom_jitter(width=0.2, height=0, size=3, stroke=1, colour='black', fill='gray') +
  scale_fill_manual(values=c("#00c703", "#cc00cc")) +
  scale_color_manual(values=c("#00c703", "#cc00cc")) +
  theme_bw() +
  ggtitle("Generation of excitatory neurons in cultured slices") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker") +
  ylab("Ratio of BrdU+ cells") +
  ylim(c(-0.002,0.085))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
image = ggplot(data = brdu_avg_counts_revision, aes(x=stain, y=double_positive_ratio_of_BrdU, fill=stain, color=stain)) +
  geom_bar(stat="summary", alpha=0.6) +
  #geom_jitter(width=0.1, height=0, size=3, stroke=1, colour='black', fill='gray') +
  geom_jitter(width=0.2, height=0, size=3, stroke=1, colour='black', fill='gray') +
  scale_fill_manual(values=c("#00c703", "#cc00cc")) +
  scale_color_manual(values=c("#00c703", "#cc00cc")) +
  theme_bw() +
  ggtitle("Generation of excitatory neurons in cultured slices") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker") +
  ylab("Ratio of BrdU+ cells") +
  ylim(c(-0.002,0.085))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
ggsave('brdu_tbr1_satb2_quant_avg_barplot_revision.svg', image, width=10, height=8, device = 'svg', dpi = 300)


ggplot(data = brdu_raw_counts_revision, aes(x=Stain, y=double_positive_ratio_of_BrdU, fill=Sample, color=Sample)) +
  #geom_bar(stat="summary", alpha=0.6) +
  #geom_jitter(width=0.1, height=0, size=3, stroke=1, colour='black', fill='gray') +
  geom_jitter(width=0.2, height=0, size=3, stroke=1) +
  #scale_fill_manual(values=c("#00c703", "#cc00cc")) +
  #scale_color_manual(values=c("#00c703", "#cc00cc")) +
  theme_bw() +
  ggtitle("Generation of excitatory neurons in cultured slices") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker") +
  ylab("Ratio of BrdU+ cells") +
  ylim(c(-0.002,0.085))+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
image = ggplot(data = brdu_raw_counts_revision, aes(x=Stain, y=double_positive_ratio_of_BrdU, fill=Sample, color=Sample)) +
  #geom_bar(stat="summary", alpha=0.6) +
  #geom_jitter(width=0.1, height=0, size=3, stroke=1, colour='black', fill='gray') +
  geom_jitter(width=0.2, height=0, size=3, stroke=1) +
  #scale_fill_manual(values=c("#00c703", "#cc00cc")) +
  #scale_color_manual(values=c("#00c703", "#cc00cc")) +
  theme_bw() +
  ggtitle("Generation of excitatory neurons in cultured slices") +
  #labs(color = "Sample ages", fill ="Sample ages") +
  xlab("Marker") +
  ylab("Ratio of BrdU+ cells") +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), #note: doing this with legend title to make space to put real legend title in
        plot.title = element_text(size=30, hjust=0.5),
        text=element_text(family="Helvetica"))
ggsave('brdu_tbr1_satb2_quant_raw_dotplot_colored_by_sample_revision.svg', image, width=10, height=8, device = 'svg', dpi = 300)



#### Figure ED1
### Fig ED1a -- single cell qc by sample
ls_synthesis = SetIdent(ls_synthesis, value = "sample_id")
ls_synthesis@active.ident = factor(ls_synthesis@active.ident, levels=c(
  "GW16",
  "GW18",
  "GW19",
  "GW20_PFC",
  "GW20_V1",
  "GW21",
  "GW22",
  "GW23",
  "GW24"
))
VlnPlot(ls_synthesis, features=c("percent.mt"), pt.size=0) + 
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
                                                                                 legend.title = element_text(size=25), legend.text = element_text(size=25), 
                                                                                 plot.title = element_text(size=30, hjust=0.5),
                                                                                 axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = VlnPlot(ls_synthesis, features=c("percent.mt"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave(file='qc_violin_mt_pct_by_sample.svg', plot=image, width=12, height=8)

VlnPlot(ls_synthesis, features=c("nCount_RNA"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = VlnPlot(ls_synthesis, features=c("nCount_RNA"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave(file='qc_violin_ncount_rna_by_sample.svg', plot=image, width=12, height=8)

VlnPlot(ls_synthesis, features=c("nFeature_RNA"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = VlnPlot(ls_synthesis, features=c("nFeature_RNA"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave(file='qc_violin_nfeature_rna_by_sample.svg', plot=image, width=12, height=8)

### Fig ED1b -- single cell qc by cell type
ls_synthesis = SetIdent(ls_synthesis, value = "subcluster_identity")
table(ls_synthesis@meta.data$subcluster_identity)
ls_synthesis@active.ident = factor(ls_synthesis@active.ident, levels=c(
  "RG",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "EX_IPCs",
  "ENs",
  "IN_IPCs",
  "IN_local",
  "IN_CGE",
  "IN_OB",
  "IN_MGE"
))
VlnPlot(ls_synthesis, features=c("percent.mt"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = VlnPlot(ls_synthesis, features=c("percent.mt"), pt.size=0)+
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave(file='qc_violin_mt_pct.svg', plot=image, width=12, height=8)

VlnPlot(ls_synthesis, features=c("nCount_RNA"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = VlnPlot(ls_synthesis, features=c("nCount_RNA"), pt.size=0) + 
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave(file='qc_violin_ncount_rna.svg', plot=image, width=12, height=8)

VlnPlot(ls_synthesis, features=c("nFeature_RNA"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = VlnPlot(ls_synthesis, features=c("nFeature_RNA"), pt.size=0) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=25, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=25), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave(file='qc_violin_nfeature_rna.svg', plot=image, width=12, height=8)


## Fig ED1c -- heatmap of gene expression by subcluster
ls_synthesis = SetIdent(ls_synthesis, value = "subcsluter_identity")
ls_synthesis@active.ident = factor(ls_synthesis@active.ident, levels=c(
  "RG",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "EX_IPCs",
  "ENs",
  "IN_IPCs",
  "IN_local",
  "IN_CGE",
  "IN_OB",
  "IN_MGE"
))

ls_synthesis.markers <- FindAllMarkers(ls_synthesis, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(ls_synthesis.markers)
ls_synthesis_clustermarkers_unfilt = as.data.frame(ls_synthesis.markers %>% group_by(cluster))
write.csv(ls_synthesis_clustermarkers_unfilt,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_synthesis_clustermarkers_v1_unfilt.csv',row.names=F)
ls_synthesis.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
ls_synthesis_clustermarkers = as.data.frame(ls_synthesis.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
write.csv(ls_synthesis_clustermarkers,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_synthesis_clustermarkers_v1.csv',row.names=F)

ls_synthesis_top10_clustermarkers = as.data.frame(ls_synthesis.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))

setwd('/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/local_STICR_synthesis_figs/240118_figs')
random_cell_list = c()
for(sc in unique(ls_synthesis@meta.data$subcsluter_identity)){
  cells_in_sc = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$subcsluter_identity %in% sc]
  random_cell_list = c(random_cell_list, sample(cells_in_sc, size = min(length(cells_in_sc),1000), replace = F))
}
DoHeatmap(ls_synthesis, features=ls_synthesis_top10_clustermarkers$gene, cells = random_cell_list)
image = DoHeatmap(ls_synthesis, features=ls_synthesis_top10_clustermarkers$gene, cells = random_cell_list)
ggsave('all_clusters_heatmap_top10_genes_1000cells.svg', image, width=8, height=12)
ggsave('all_clusters_heatmap_top10_genes_1000cells.png', image, width=8, height=12)


## Fig ED1d -- UMAPs of marker genes

FeaturePlot(ls_synthesis, features=c("HES1"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image =FeaturePlot(ls_synthesis, features=c("HES1"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_HES1_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("MKI67"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("MKI67"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_MKI67_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("EOMES"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("EOMES"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_EOMES_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("SLC17A7"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("SLC17A7"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_SLC17A7_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("GAD2"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("GAD2"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_GAD2_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("SPARCL1"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("SPARCL1"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_SPARCL1_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("OLIG2"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("OLIG2"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_OLIG2_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("MBP"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("MBP"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_MBP_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis, features=c("LHX6"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("LHX6"), order=T, pt.size=1, min.cutoff='q50') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_umap_LHX6_large_dots_q50_ordered.png", plot=image, width=4, height=3.2)



## Fig ED1e-f -- pie charts of sample composition by cell types
ls_synthesis@meta.data$subcluster_identity = factor(ls_synthesis@meta.data$subcluster_identity, levels=c(
  "RG",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "EX_IPCs",
  "ENs",
  "IN_IPCs",
  "IN_local",
  "IN_CGE",
  "IN_OB",
  "IN_MGE"
))

pie_chart_age = "<=GW20"
pie_chart_df = data.frame((table(subset(ls_synthesis@meta.data, age_binned %in% pie_chart_age)$subcluster_identity))/sum(table(subset(ls_synthesis@meta.data, age_binned %in% pie_chart_age)$subcluster_identity)))
pie_chart_df$Freq = round((pie_chart_df$Freq*100), 1)

ggplot(pie_chart_df, aes(x=0, y=Freq, fill=Var1, color="white", label=paste(Freq, "%", sep=""))) +
  geom_bar(stat="identity", position="fill") +
  geom_text_repel(size=6, position = position_fill(vjust=0.5), color="black") +
  scale_color_manual(values = "white") +
  coord_polar("y", start=0) +
  labs(fill="Cell type") +
  ggtitle(paste("Cluster composition", pie_chart_age)) +
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size=30, hjust=0.5))
image = ggplot(pie_chart_df, aes(x=0, y=Freq, fill=Var1, color="white", label=paste(Freq, "%", sep=""))) +
  geom_bar(stat="identity", position="fill") +
  geom_text_repel(size=6, position = position_fill(vjust=0.5), color="black") +
  scale_color_manual(values = "white") +
  coord_polar("y", start=0) +
  labs(fill="Cell type") +
  ggtitle(paste("Cluster composition", pie_chart_age)) +
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size=30, hjust=0.5))
ggsave('pie_chart_young.svg', image, width=10, height=8)

pie_chart_age = ">GW20"
pie_chart_df = data.frame((table(subset(ls_synthesis@meta.data, age_binned %in% pie_chart_age)$subcluster_identity))/sum(table(subset(ls_synthesis@meta.data, age_binned %in% pie_chart_age)$subcluster_identity)))
pie_chart_df$Freq = round((pie_chart_df$Freq*100), 1)
ggplot(pie_chart_df, aes(x=0, y=Freq, fill=Var1, color="white", label=paste(Freq, "%", sep=""))) +
  geom_bar(stat="identity", position="fill") +
  geom_text_repel(size=6, position = position_fill(vjust=0.5), color="black") +
  scale_color_manual(values = "white") +
  coord_polar("y", start=0) +
  labs(fill="Cell type") +
  ggtitle(paste("Cluster composition", pie_chart_age)) +
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size=30, hjust=0.5))
image = ggplot(pie_chart_df, aes(x=0, y=Freq, fill=Var1, color="white", label=paste(Freq, "%", sep=""))) +
  geom_bar(stat="identity", position="fill") +
  geom_text_repel(size=6, position = position_fill(vjust=0.5), color="black") +
  scale_color_manual(values = "white") +
  coord_polar("y", start=0) +
  labs(fill="Cell type") +
  ggtitle(paste("Cluster composition", pie_chart_age)) +
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size=30, hjust=0.5))
ggsave('pie_chart_old.svg', image, width=10, height=8)

pie_chart_age = c("<=GW20", ">GW20")
pie_chart_df = data.frame((table(subset(ls_synthesis@meta.data, age_binned %in% pie_chart_age)$subcluster_identity))/sum(table(subset(ls_synthesis@meta.data, age_binned %in% pie_chart_age)$subcluster_identity)))
pie_chart_df$Freq = round((pie_chart_df$Freq*100), 1)
ggplot(pie_chart_df, aes(x=0, y=Freq, fill=Var1, color="white", label=paste(Freq, "%", sep=""))) +
  geom_bar(stat="identity", position="fill") +
  geom_text_repel(size=6, position = position_fill(vjust=0.5), color="black") +
  scale_color_manual(values = "white") +
  coord_polar("y", start=0) +
  labs(fill="Cell type") +
  ggtitle(paste("Cluster composition", "all")) +
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size=30, hjust=0.5))
image = ggplot(pie_chart_df, aes(x=0, y=Freq, fill=Var1, color="white", label=paste(Freq, "%", sep=""))) +
  geom_bar(stat="identity", position="fill") +
  geom_text_repel(size=6, position = position_fill(vjust=0.5), color="black") +
  scale_color_manual(values = "white") +
  coord_polar("y", start=0) +
  labs(fill="Cell type") +
  ggtitle(paste("Cluster composition", "all")) +
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size=30, hjust=0.5))
ggsave('pie_chart_all_ages.svg', image, width=10, height=8)

## ED2a -- barcode collisions across samples
ls_synthesis_meta = ls_synthesis@meta.data
ls_synthesis_meta$Clone_barcode_compare = ls_synthesis_meta$Clone_barcode
ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW24" & !is.na(ls_synthesis_meta$Clone_barcode_compare)] = paste(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW24" & !is.na(ls_synthesis_meta$Clone_barcode_compare)], "t", sep="")

gw16_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW16" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw18_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW18" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw19_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW19" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw20_pfc_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW20_PFC" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw20_v1_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW20_V1" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw21_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW21" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw22_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW22" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw23_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW23" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw24_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW24" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
length(c(gw16_bcs,gw19_bcs,gw20_pfc_bcs,gw20_v1_bcs,gw21_bcs,gw22_bcs,gw23_bcs,gw24_bcs))
# [1] 46243
length(unique(c(gw16_bcs,gw19_bcs,gw20_pfc_bcs,gw20_v1_bcs,gw21_bcs,gw22_bcs,gw23_bcs,gw24_bcs)))
# [1] 46226

bc_collision_df = data.frame("Clone_barcode_compare" = c(gw16_bcs, gw19_bcs, gw20_pfc_bcs, gw20_v1_bcs, gw21_bcs, gw22_bcs, gw23_bcs, gw24_bcs),
                             "sample_id" = c(rep("GW16", length(gw16_bcs)),
                                             rep("GW19", length(gw19_bcs)),
                                             rep("GW20_PFC", length(gw20_pfc_bcs)),
                                             rep("GW20_V1", length(gw20_v1_bcs)),
                                             rep("GW21", length(gw21_bcs)),
                                             rep("GW22", length(gw22_bcs)),
                                             rep("GW23", length(gw23_bcs)),
                                             rep("GW24", length(gw24_bcs))))

sample_idents = c("GW16","GW19", "GW20_PFC", "GW20_V1", "GW21", "GW22", "GW23", "GW24")
length(sample_idents)
bc_matrix_counts = matrix(rep(0, length(sample_idents)*length(sample_idents)), ncol=length(sample_idents))
dim(bc_matrix_counts)
rownames(bc_matrix_counts) = sample_idents
colnames(bc_matrix_counts) = sample_idents

for (sample_id in sample_idents) {
  temp_collision_count = c()
  for (sample_id2 in sample_idents){
    bc_collision_df_filt = bc_collision_df[bc_collision_df$sample_id %in% c(sample_id, sample_id2),]
    n_collisions = length(bc_collision_df_filt$Clone_barcode_compare) - length(unique(bc_collision_df_filt$Clone_barcode_compare))
    temp_collision_count = c(temp_collision_count, n_collisions)
  }
  bc_matrix_counts[rownames(bc_matrix_counts) %in% sample_id,] = as.numeric(temp_collision_count)
}
bc_matrix_counts

heatmap(bc_matrix_counts, Rowv = NA, Colv = NA)


## ED2a -- barcode collisions across samples
ls_synthesis_meta = ls_synthesis@meta.data
ls_synthesis_meta$Clone_barcode_compare = ls_synthesis_meta$Clone_barcode
ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW24" & !is.na(ls_synthesis_meta$Clone_barcode_compare)] = paste(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW24" & !is.na(ls_synthesis_meta$Clone_barcode_compare)], "t", sep="")

gw16_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW16" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw18_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW18" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw19_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW19" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw20_pfc_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW20_PFC" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw20_v1_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW20_V1" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw21_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW21" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw22_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW22" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw23_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW23" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
gw24_bcs = unique(ls_synthesis_meta$Clone_barcode_compare[ls_synthesis_meta$sample_id %in% "GW24" & !is.na(ls_synthesis_meta$Clone_barcode_compare)])
length(c(gw16_bcs,gw19_bcs,gw20_pfc_bcs,gw20_v1_bcs,gw21_bcs,gw22_bcs,gw23_bcs,gw24_bcs))
# [1] 46243
length(unique(c(gw16_bcs,gw19_bcs,gw20_pfc_bcs,gw20_v1_bcs,gw21_bcs,gw22_bcs,gw23_bcs,gw24_bcs)))
# [1] 46226

bc_collision_df = data.frame("Clone_barcode_compare" = c(gw16_bcs, gw19_bcs, gw20_pfc_bcs, gw20_v1_bcs, gw21_bcs, gw22_bcs, gw23_bcs, gw24_bcs),
                             "sample_id" = c(rep("GW16", length(gw16_bcs)),
                                             rep("GW19", length(gw19_bcs)),
                                             rep("GW20_PFC", length(gw20_pfc_bcs)),
                                             rep("GW20_V1", length(gw20_v1_bcs)),
                                             rep("GW21", length(gw21_bcs)),
                                             rep("GW22", length(gw22_bcs)),
                                             rep("GW23", length(gw23_bcs)),
                                             rep("GW24", length(gw24_bcs))))

sample_idents = c("GW16","GW19", "GW20_PFC", "GW20_V1", "GW21", "GW22", "GW23", "GW24")
length(sample_idents)
bc_matrix_counts = matrix(rep(0, length(sample_idents)*length(sample_idents)), ncol=length(sample_idents))
dim(bc_matrix_counts)
rownames(bc_matrix_counts) = sample_idents
colnames(bc_matrix_counts) = sample_idents

for (sample_id in sample_idents) {
  temp_collision_count = c()
  for (sample_id2 in sample_idents){
    bc_collision_df_filt = bc_collision_df[bc_collision_df$sample_id %in% c(sample_id, sample_id2),]
    n_collisions = length(bc_collision_df_filt$Clone_barcode_compare) - length(unique(bc_collision_df_filt$Clone_barcode_compare))
    temp_collision_count = c(temp_collision_count, n_collisions)
  }
  bc_matrix_counts[rownames(bc_matrix_counts) %in% sample_id,] = as.numeric(temp_collision_count)
}
bc_matrix_counts

heatmap(bc_matrix_counts, Rowv = NA, Colv = NA)

## Fig2b-c -- Sankey plot for barcode recovery
library(ggsankey)
sankey_barcode_recovery_metadata = ls_synthesis@meta.data
colnames(ls_synthesis@meta.data)
sum(!is.na(sankey_barcode_recovery_metadata$Clone_IDs))
# [1] 58681
sankey_barcode_recovery_metadata$cell = "Cell"

table(clone_df$Clone_size)
# 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    19    22    23    25 
# 39851  3700  1357   610   300   179    85    59    27    26    23    13     6     4     2     6     1     1     1     1     1
table(clone_df_v2$Clone_size)
# 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    19    22    23    25 
# 39851  3700  1357   610   300   179    85    59    27    26    23    13     6     4     2     6     1     1     1     1     1

length(multicell_clones)
#[1] 6402

sankey_barcode_recovery_metadata$barcode_recovered = NA
sankey_barcode_recovery_metadata$barcode_recovered[sankey_barcode_recovery_metadata$Clone_IDs %in% clone_df$Clone_IDs] = "Barcode recovered"
sankey_barcode_recovery_metadata$barcode_recovered[!(sankey_barcode_recovery_metadata$Clone_IDs %in% clone_df$Clone_IDs)] = "Barcode lost"
table(sankey_barcode_recovery_metadata$barcode_recovered)
# Barcode lost Barcode recovered
# 38859             58681
colnames(sankey_barcode_recovery_metadata)
table(sankey_barcode_recovery_metadata$multi_clone)
# multi_cell_clone single_cell_clone 
# 18830  39851

sankey_barcode_recovery = sankey_barcode_recovery_metadata %>% make_long(cell, barcode_recovered, age_binned, multi_clone)
ggplot(sankey_barcode_recovery, aes(x = x, next_x = next_x, node = node, next_node = next_node,fill = factor(node), label=node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = TRUE) +
  geom_sankey_label(Size = 3, 
                    color = "black", 
                    fill = "white", ) +
  theme_bw() +
  labs(fill='Nodes') +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), 
        legend.position = "none", plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(sankey_barcode_recovery, aes(x = x, next_x = next_x, node = node, next_node = next_node,fill = factor(node), label=node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = TRUE) +
  geom_sankey_label(Size = 3, 
                    color = "black", 
                    fill = "white", ) +
  theme_bw() +
  labs(fill='Nodes') +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), 
        legend.position = "none", plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("sankey_barcode_recovery_by_cells.svg", image, width=10, height=8)

table(sankey_barcode_recovery_metadata$barcode_recovered)
#     Barcode lost Barcode recovered 
#     38859             58681
table(sankey_barcode_recovery_metadata$age_binned[sankey_barcode_recovery_metadata$barcode_recovered %in% "Barcode recovered"])
#<=GW20  >GW20 
#38128  20553
table(sankey_barcode_recovery_metadata$multi_clone[sankey_barcode_recovery_metadata$age_binned %in% "<=GW20" & sankey_barcode_recovery_metadata$barcode_recovered %in% "Barcode recovered"])
# multi_cell_clone single_cell_clone 
# 12566  25562
table(sankey_barcode_recovery_metadata$multi_clone[sankey_barcode_recovery_metadata$age_binned %in% ">GW20" & sankey_barcode_recovery_metadata$barcode_recovered %in% "Barcode recovered"])
# multi_cell_clone single_cell_clone 
# 6264  14289

table(ls_synthesis@meta.data$Clone_size[ls_synthesis@meta.data$Clone_IDs %in% sankey_barcode_recovery_metadata$Clone_IDs[sankey_barcode_recovery_metadata$multi_clone %in% "multi_cell_clone"]])
#    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   19   22   23   25 
# 7400 4071 2440 1500 1074  595  472  243  260  253  156   78   56   30   96   17   19   22   23   25

## repeat for clones instead of cells
## starting with the method that was used in the first submission
# sankey_barcode_recovery_metadata = cortex_midgestation@meta.data
# head(sankey_barcode_recovery_metadata)
dim(sankey_barcode_recovery_metadata[!(is.na(sankey_barcode_recovery_metadata$Clone_IDs)),])
#[1] 58681    17
sankey_barcode_recovery_metadata = sankey_barcode_recovery_metadata[!is.na(sankey_barcode_recovery_metadata$Clone_IDs),]
dim(sankey_barcode_recovery_metadata)
#[1] 58681    17
dim(sankey_barcode_recovery_metadata[!duplicated(sankey_barcode_recovery_metadata$Clone_IDs),])
#[1] 46253    17
sankey_barcode_recovery_metadata = sankey_barcode_recovery_metadata[!duplicated(sankey_barcode_recovery_metadata$Clone_IDs),]
table(sankey_barcode_recovery_metadata$multi_clone)
#multi_cell_clone single_cell_clone 
# 6402  39851 
sankey_barcode_recovery_metadata$clones = "Clone recovered"
sankey_barcode_recovery = sankey_barcode_recovery_metadata %>% make_long(clones, age_binned, multi_clone)
ggplot(sankey_barcode_recovery, aes(x = x, next_x = next_x, node = node, next_node = next_node,fill = factor(node), label=node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = TRUE) +
  geom_sankey_label(Size = 3, 
                    color = "black", 
                    fill = "white", ) +
  theme_bw() +
  labs(fill='Nodes') +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), 
        legend.position = "none", plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(sankey_barcode_recovery, aes(x = x, next_x = next_x, node = node, next_node = next_node,fill = factor(node), label=node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = TRUE) +
  geom_sankey_label(Size = 3, 
                    color = "black", 
                    fill = "white", ) +
  theme_bw() +
  labs(fill='Nodes') +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), 
        legend.position = "none", plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("sankey_barcode_recovery_by_clones.svg", image, width=10, height=8)

table(sankey_barcode_recovery_metadata$age_binned)
# <=GW20  >GW20 
# 29771  16482
table(sankey_barcode_recovery_metadata$multi_clone[sankey_barcode_recovery_metadata$age_binned%in% ">GW20"])
# multi single 
# 2193  14289 
table(sankey_barcode_recovery_metadata$multi_clone[sankey_barcode_recovery_metadata$age_binned %in% "<=GW20"])
# multi single 
# 4209  25562


## Fig ED2d-e -- UMAP of cells with recovered barcodes or in multicellular clones
DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[!is.na(ls_synthesis_meta$Clone_IDs)], sizes.highlight = 0.5) +
  scale_color_manual(values=c("lightgray","limegreen")) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[!is.na(ls_synthesis_meta$Clone_IDs)], sizes.highlight = 0.5) +
  scale_color_manual(values=c("lightgray","limegreen")) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave('full_umap_cells_with_barcodes.png', image, width=10, height=8)
ggsave('full_umap_cells_with_barcodes.svg', image, width=10, height=8)

DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis_meta$Clone_IDs %in% multicell_clones], sizes.highlight = 0.5) +
  scale_color_manual(labels = c("Clone size <= 1", "Clone size >= 2"), values=c("lightgray","darkgreen")) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis_meta$Clone_IDs %in% multicell_clones], sizes.highlight = 0.5) +
  scale_color_manual(labels = c("Clone size <= 1", "Clone size >= 2"), values=c("lightgray","darkgreen")) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave('full_umap_multicellular_clones.png', image, width=10, height=8)
ggsave('full_umap_multicellular_clones.svg', image, width=10, height=8)

## ED2f-h -- histograms of clone size by sample
ggplot(clone_df_v2[clone_df_v2$Clone_size>=2,], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  xlim(c(1,27)) +
  scale_y_log10(limits=c(NA, 4000)) +
  ggtitle("Clone Size all samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
image = ggplot(clone_df_v2[clone_df_v2$Clone_size>=2,], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  xlim(c(1,27)) +
  scale_y_log10(oob = scales::oob_squish) +
  ggtitle("Clone Size all samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("histogram_clone_size_all_cells_log_scale.svg", image, width=10, height=8)

ggplot(clone_df_v2[clone_df_v2$Clone_size>=2,], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  xlim(c(1,27)) +
  ylim(0,100) +
  # scale_y_log10() +
  ggtitle("Clone Size all samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
image = ggplot(clone_df_v2[clone_df_v2$Clone_size>=2,], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  xlim(c(1,27)) +
  ylim(0,100) +
  # scale_y_log10(oob = scales::oob_squish) +
  ggtitle("Clone Size all samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("histogram_clone_size_all_cells_linear_scale.svg", image, width=10, height=8)

clone_df_v2 = data.frame(clone_df)
ls_synthesis_meta = ls_synthesis@meta.data
ls_synthesis_meta_unique_clones = ls_synthesis_meta
ls_synthesis_meta_unique_clones = ls_synthesis_meta_unique_clones %>% distinct(Clone_IDs, .keep_all = TRUE)
ls_synthesis_meta_unique_clones = ls_synthesis_meta_unique_clones[!is.na(ls_synthesis_meta_unique_clones$Clone_IDs),]
dim(ls_synthesis_meta_unique_clones)
# [1] 46253    15

head(clone_df_v2)
dim(clone_df_v2)
# [1] 46253     2
sum(clone_df_v2$Clone_IDs %in% ls_synthesis_meta$Clone_IDs)
# [1] 46253
dim(left_join(clone_df_v2, ls_synthesis_meta_unique_clones[,c("Clone_IDs", "GZ_origin", "age_binned")], by=c("Clone_IDs")))
# [1] 46253     4
clone_df_v2 = left_join(clone_df_v2, ls_synthesis_meta_unique_clones[,c("Clone_IDs", "GZ_origin", "age_binned")], by=c("Clone_IDs"))
head(clone_df_v2)
table(clone_df_v2$age_binned)
# <=GW20  >GW20
# 29771  16482
clone_df_v2 <- clone_df_v2 %>%
  distinct(Clone_IDs, .keep_all = TRUE)
table(clone_df_v2$age_binned)
# <=GW20  >GW20
# 29771  16482

ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% "<=GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  scale_y_log10(oob = scales::oob_squish, limits=c(NA, 4000)) +
  ggtitle("Clone Size young samples") +
  xlim(c(1,27)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
table(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% "<=GW20",]$Clone_size)
#    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   19   22   23 
# 2435  906  369  180  114   62   44   23   21   21   13    5    4    2    6    1    1    1    1
image = ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% "<=GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  scale_y_log10(oob = scales::oob_squish, limits=c(NA, 4000)) +
  ggtitle("Clone Size young samples") +
  xlim(c(1,27)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("histogram_clone_size_young_log_scale.svg", image, width=10, height=8)

ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% "<=GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  # scale_y_log10(oob = scales::oob_squish) +
  ggtitle("Clone Size all samples") +
  xlim(c(1,27)) +
  ylim(0,100) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
image = ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% "<=GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  # scale_y_log10(oob = scales::oob_squish) +
  ggtitle("Clone Size all samples") +
  xlim(c(1,27)) +
  ylim(0,100) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("histogram_clone_size_young_linear_scale.svg", image, width=10, height=8)

ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% ">GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  scale_y_log10(oob = scales::oob_squish, limits=c(NA, 4000)) +
  ggtitle("Clone Size old samples") +
  xlim(c(1,27)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
table(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% ">GW20",]$Clone_size)
#    2    3    4    5    6    7    8    9   10   11   13   25 
# 1265  451  241  120   65   23   15    4    5    2    1    1
image = ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% ">GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  scale_y_log10(oob = scales::oob_squish, limits=c(NA, 4000)) +
  ggtitle("Clone Size old samples") +
  xlim(c(1,27)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("histogram_clone_size_old_log_scale.svg", image, width=10, height=8)

ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% ">GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  # scale_y_log10(oob = scales::oob_squish) +
  ggtitle("Clone Size all samples") +
  xlim(c(1,27)) +
  ylim(0,100) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
image = ggplot(clone_df_v2[clone_df_v2$Clone_size>=2 & clone_df_v2$age_binned %in% ">GW20",], aes(x = Clone_size)) +
  geom_histogram(position = "dodge", binwidth=1, ) +
  # scale_fill_manual(values=c(osvz_col, vz_col)) +
  xlab("Clone Size") +
  ylab("Counts") +
  # scale_y_log10(oob = scales::oob_squish) +
  ggtitle("Clone Size all samples") +
  xlim(c(1,27)) +
  ylim(0,100) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("histogram_clone_size_old_linear_scale.svg", image, width=10, height=8)


## Figure ED3
## ED3a -- subsetting big UMAP to progenitor clusters
FeaturePlot(ls_synthesis, features=c("HES1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("HES1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_HES1.png", image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("EOMES")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("EOMES")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_EOMES.png", image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("OLIG2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("OLIG2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_OLIG2.png", image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_MKI67.png", image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("NPY")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("NPY")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_NPY.png", image, width=4, height=3.2)

DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[rownames(ls_synthesis@meta.data) %in% rownames(ls_prog@meta.data)]) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
image = DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[rownames(ls_synthesis@meta.data) %in% rownames(ls_prog@meta.data)]) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
ggsave("ls_full_umap_ls_prog_highlighted.png", image, width=4, height=3.2)

## ED3b -- subsetting ls_prog to ls_synthesis_prog_glia
FeaturePlot(ls_prog, features=c("FOXG1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("FOXG1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_FOXG1.png", image, width=4, height=3.2)
FeaturePlot(ls_prog, features=c("HES1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("HES1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_HES1.png", image, width=4, height=3.2)
FeaturePlot(ls_prog, features=c("VIM")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("VIM")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_VIM.png", image, width=4, height=3.2)
FeaturePlot(ls_prog, features=c("OLIG2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("OLIG2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_OLIG2.png", image, width=4, height=3.2)
FeaturePlot(ls_prog, features=c("EMX2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("EMX2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_EMX2.png", image, width=4, height=3.2)
FeaturePlot(ls_prog, features=c("DLX2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("DLX2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_DLX2.png", image, width=4, height=3.2)
FeaturePlot(ls_prog, features=c("EOMES")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("EOMES")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_EOMES.png", image, width=4, height=3.2)
FeaturePlot(ls_prog, features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_prog, features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_prog_umap_MKI67.png", image, width=4, height=3.2)

DimPlot(ls_prog, cells.highlight = rownames(ls_prog@meta.data)[rownames(ls_prog@meta.data) %in% rownames(ls_synthesis_prog_glia@meta.data)]) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
image = DimPlot(ls_prog, cells.highlight = rownames(ls_prog@meta.data)[rownames(ls_prog@meta.data) %in% rownames(ls_synthesis_prog_glia@meta.data)]) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
ggsave("ls_prog_umap_ls_synthesis_prog_glia_highlighted.png", image, width=4, height=3.2)


## ED3c -- UMAP of ls_synthesis_prog_glia with seurat clusters
DimPlot(ls_synthesis_prog_glia, pt.size=1.5, group.by="seurat_clusters", label=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30))
image = DimPlot(ls_synthesis_prog_glia, pt.size=1.5, group.by="seurat_clusters", label=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30))
ggsave("ls_synthesis_prog_glia_umap_seurat_clusters.png", image, width=10, height=8)
ggsave("ls_synthesis_prog_glia_umap_seurat_clusters.svg", image, width=10, height=8)

## ED3d -- UMAP of ls_synthesis_prog_glia by sample
DimPlot(ls_synthesis_prog_glia, pt.size=1, group.by="sample_id") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30), plot.title = element_blank())
image = DimPlot(ls_synthesis_prog_glia, pt.size=1, group.by="sample_id") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30), plot.title = element_blank())
ggsave("ls_synthesis_prog_glia_umap_grouped_sample.png", image, width=10, height=8)
ggsave("ls_synthesis_prog_glia_umap_grouped_sample.svg", image, width=10, height=8)

## ED3e -- UMAP of ls_synthesis_prog_glia by subcluster
table(ls_synthesis_prog_glia@meta.data$subcluster_identity)
ls_synthesis_prog_glia@meta.data$subcluster_identity = factor(ls_synthesis_prog_glia@meta.data$subcluster_identity, levels=c(
  "bRG",
  "tRG",
  "Early OPCs",
  "Astrocytes (dense smooth)",
  "Astrocytes (dense bulbous)",
  "Dividing", 
  "EX IPCs"))
DimPlot(ls_synthesis_prog_glia, pt.size=1.5, group.by="subcluster_identity") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30), plot.title = element_blank())
image = DimPlot(ls_synthesis_prog_glia, pt.size=1.5, group.by="subcluster_identity") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30), plot.title = element_blank())
ggsave("ls_synthesis_prog_glia_umap_grouped_subclusts.png", image, width=10, height=8)
ggsave("ls_synthesis_prog_glia_umap_grouped_subclusts.svg", image, width=10, height=8)


## ED3f -- marker gene expression for ls_synthesis_prog_glia
FeaturePlot(ls_synthesis_prog_glia, features=c("HES1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("HES1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_prog_glia_umap_HES1.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_HES1.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_prog_glia, features=c("VIM")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("VIM")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_prog_glia_umap_VIM.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_VIM.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_prog_glia, features=c("INPP1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("INPP1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_prog_glia_umap_INPP1.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_INPP1.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_prog_glia, features=c("ANXA1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("ANXA1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_prog_glia_umap_ANXA1.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_ANXA1.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_prog_glia, features=c("CRYAB")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("CRYAB")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_prog_glia_umap_CRYAB.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_CRYAB.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_prog_glia, features=c("PDGFRA")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("PDGFRA")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_prog_glia_umap_PDGFRA.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_PDGFRA.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_prog_glia, features=c("GJA1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("GJA1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_prog_glia_umap_GJA1.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_GJA1.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_prog_glia, features=c("SPARCL1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis_prog_glia, features=c("SPARCL1")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_prog_glia_umap_SPARCL1.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_SPARCL1.png", plot=image, width=4, height=3.2)

FeaturePlot(subset(ls_synthesis_prog_glia, seurat_clusters %in% c(5), invert=T), features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(subset(ls_synthesis_prog_glia, seurat_clusters %in% c(5), invert=T), features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave(file="ls_synthesis_prog_glia_umap_MKI67.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_prog_glia_umap_MKI67.png", plot=image, width=4, height=3.2)

## ED3g -- dotplot of glial gene marker expression
DotPlot(ls_synthesis_prog_glia, features=c("PAX6", "FOXG1", "EMX2", 
                                 "OTX2", "LHX6",
                                 "VIM", "HES1", "AQP4", "GFAP", "GLI3", "NES", "HOPX",
                                 "CRYAB", "ANXA1",
                                 "INPP1", "PPM1K",
                                 "APOE", "CD44", "SPARCL1", "GJA1",
                                 "TIMP1", "S100A11", "EMP3",
                                 "ANGPTL4", "TIMP3", "NTRK2",
                                 "MKI67", "CENPE",
                                 "OLIG2", "PDGFRA", "EOMES")) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = DotPlot(ls_synthesis_prog_glia, features=c("PAX6", "FOXG1", "EMX2", 
                                         "OTX2", "LHX6",
                                         "VIM", "HES1", "AQP4", "GFAP", "GLI3", "NES", "HOPX",
                                         "CRYAB", "ANXA1",
                                         "INPP1", "PPM1K",
                                         "APOE", "CD44", "SPARCL1", "GJA1",
                                         "TIMP1", "S100A11", "EMP3",
                                         "ANGPTL4", "TIMP3", "NTRK2",
                                         "MKI67", "CENPE",
                                         "OLIG2", "PDGFRA", "EOMES")) +
  theme(axis.text.y = element_text(size=2), axis.text.x = element_text(size=20, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("ls_synthesis_prog_glia_dotplot_by_defined_subtype_revision_astrocyte_subtypes.svg", image, width=12, height=8)







## Figure ED4
## ED4b -- OPCs
pdf(paste(savepath,'/upset_plot_OPCs_all.pdf', sep=""), width=12, height=8)
focus_age = c("<=GW20", ">GW20")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$OPCs>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_OPCs_young.pdf', sep=""), width=12, height=8)
focus_age = "<=GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$OPCs>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_OPCs_old.pdf', sep=""), width=12, height=8)
focus_age = ">GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$OPCs>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

## ED4c -- Astrocytes
pdf(paste(savepath,'/upset_plot_Astrocytes_all.pdf', sep=""), width=12, height=8)
focus_age = c("<=GW20", ">GW20")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$Astrocytes>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_Astrocytes_young.pdf', sep=""), width=12, height=8)
focus_age = "<=GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$Astrocytes>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_Astrocytes_old.pdf', sep=""), width=12, height=8)
focus_age = ">GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$Astrocytes>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

## ED4d -- RG
pdf(paste(savepath,'/upset_plot_RG_all.pdf', sep=""), width=12, height=8)
focus_age = c("<=GW20", ">GW20")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$RG>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_RG_young.pdf', sep=""), width=12, height=8)
focus_age = "<=GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$RG>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_RG_old.pdf', sep=""), width=12, height=8)
focus_age = ">GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$RG>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

## ED4e -- IN_IPCs
pdf(paste(savepath,'/upset_plot_IN_IPCs_all.pdf', sep=""), width=12, height=8)
focus_age = c("<=GW20", ">GW20")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$IN_IPCs>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_IN_IPCs_young.pdf', sep=""), width=12, height=8)
focus_age = "<=GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$IN_IPCs>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_IN_IPCs_old.pdf', sep=""), width=12, height=8)
focus_age = ">GW20"
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$IN_IPCs>=1 & compiled_clust_df_v1$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$upset %in% multi_sample_clone_types_v1 & compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

## ED4f -- bRG
## v3 -- tRG as RG
ls_synthesis_meta = ls_synthesis@meta.data
ls_synthesis_meta$subcluster_identity_rg = as.character(ls_synthesis_meta$subcluster_identity)
ls_synthesis_meta$subcluster_identity_rg[rownames(ls_synthesis_meta) %in% rownames(ls_synthesis_prog_glia@meta.data)[ls_synthesis_prog_glia@meta.data$subcluster_identity  %in% "tRG"]] = "tRG"
ls_synthesis_meta$subcluster_identity_rg[rownames(ls_synthesis_meta) %in% rownames(ls_synthesis_prog_glia@meta.data)[ls_synthesis_prog_glia@meta.data$subcluster_identity  %in% "bRG"]] = "bRG"
table(ls_synthesis_meta$subcluster_identity_rg)

compiled_clust_df_v3 = data.frame(Clone_IDs=multicell_clones)
compiled_clust_df_v3$ENs = rep(0,length(multicell_clones))
compiled_clust_df_v3$EX_IPCs = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_local = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_CGE = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_MGE = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_OB = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_IPCs = rep(0,length(multicell_clones))
compiled_clust_df_v3$Astrocytes = rep(0,length(multicell_clones))
compiled_clust_df_v3$OPCs = rep(0,length(multicell_clones))
compiled_clust_df_v3$Oligodendrocytes = rep(0,length(multicell_clones))
compiled_clust_df_v3$RG = rep(0,length(multicell_clones))
compiled_clust_df_v3$tRG = rep(0,length(multicell_clones))
compiled_clust_df_v3$bRG = rep(0,length(multicell_clones))
compiled_clust_df_v3$ENs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$EX_IPCs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_local_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_CGE_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_MGE_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_OB_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$IN_IPCs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$Astrocytes_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$Oligodendrocytes_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$OPCS_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$RG_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$tRG_pct = rep(0,length(multicell_clones))
compiled_clust_df_v3$bRG_pct = rep(0,length(multicell_clones))

for(clone in multicell_clones) {
  clone_size = clone_df$Clone_size[clone_df$Clone_IDs %in% clone][1]
  clone_clusts = table(ls_synthesis_meta$subcluster_identity_rg[ls_synthesis_meta$Clone_IDs %in% clone])
  clone_row = which(compiled_clust_df_v3$Clone_IDs %in% clone)
  compiled_clust_df_v3[clone_row,"ENs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "ENs"])
  compiled_clust_df_v3[clone_row,"EX_IPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX_IPCs"])
  compiled_clust_df_v3[clone_row,"IN_local"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_local"])
  compiled_clust_df_v3[clone_row,"IN_CGE"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_CGE"])
  compiled_clust_df_v3[clone_row,"IN_MGE"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_MGE"])
  compiled_clust_df_v3[clone_row,"IN_OB"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_OB"])
  compiled_clust_df_v3[clone_row,"IN_IPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_IPCs"])
  compiled_clust_df_v3[clone_row,"Astrocytes"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Astrocytes"])
  compiled_clust_df_v3[clone_row,"Oligodendrocytes"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Oligodendrocytes"])
  compiled_clust_df_v3[clone_row,"OPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])
  compiled_clust_df_v3[clone_row,"RG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% c("tRG", "RG")])
  #compiled_clust_df_v3[clone_row,"tRG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "tRG"])
  compiled_clust_df_v3[clone_row,"bRG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "bRG"])
  compiled_clust_df_v3[clone_row,"ENs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "ENs"])/clone_size
  compiled_clust_df_v3[clone_row,"EX_IPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX_IPCs"])/clone_size
  compiled_clust_df_v3[clone_row,"IN_local_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_local"])/clone_size
  compiled_clust_df_v3[clone_row,"IN_CGE_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_CGE"])/clone_size
  compiled_clust_df_v3[clone_row,"IN_MGE_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_MGE"])/clone_size
  compiled_clust_df_v3[clone_row,"IN_OB_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_OB"])/clone_size
  compiled_clust_df_v3[clone_row,"IN_IPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN_IPCs"])/clone_size
  compiled_clust_df_v3[clone_row,"Astrocytes_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Astrocytes"])/clone_size
  compiled_clust_df_v3[clone_row,"Oligodendrocytes_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Oligodendrocytes"])/clone_size
  compiled_clust_df_v3[clone_row,"OPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])/clone_size
  compiled_clust_df_v3[clone_row,"RG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% c("tRG", "RG")])/clone_size
  #compiled_clust_df_v3[clone_row,"tRG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "tRG"])/clone_size
  compiled_clust_df_v3[clone_row,"bRG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "bRG"])/clone_size
}
compiled_clust_df_v3
head(compiled_clust_df_v3)
dim(compiled_clust_df_v3)
#[1] 6402   28

upset_notations = c()
for(r in 1:nrow(compiled_clust_df_v3)) {
  upset_notations=c(upset_notations, paste(colnames(compiled_clust_df_v3[r,2:14] %>% select(where(~ any(. != 0)))), collapse="&"))
}
upset_notations
compiled_clust_df_v3$upset = upset_notations
table(compiled_clust_df_v3$upset)
dim(left_join(compiled_clust_df_v3, ls_synthesis_meta[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first"))
#[1] 6402   33
compiled_clust_df_v3 = left_join(compiled_clust_df_v3, ls_synthesis_meta[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first")
head(compiled_clust_df_v3)

# drop clone types only in one sample
all_clone_types_v3 = unique(compiled_clust_df_v3$upset)
all_clone_types_sample_counts_v3 = data.frame("clone_types" = all_clone_types_v3, "sample_counts" = 0)
all_clone_types_sample_counts_v3
for (ct in all_clone_types_v3){
  all_clone_types_sample_counts_v3$sample_counts[all_clone_types_sample_counts_v3$clone_types %in% ct] = length(table(compiled_clust_df_v3$sample_id[compiled_clust_df_v3$upset %in% ct]))
}
table(all_clone_types_sample_counts_v3$sample_counts)
# 1  2  3  4  5  6  7  8 
#42 21 15 16  8 10 13  6
length(all_clone_types_sample_counts_v3$clone_types[all_clone_types_sample_counts_v3$sample_counts <= 1])
#[1] 42
length(all_clone_types_sample_counts_v3$clone_types[all_clone_types_sample_counts_v3$sample_counts > 1])
#[1] 89

## instead of using those, going to use the multi sample correction from when RG are grouped for consistency
multi_sample_clone_types_v3 = unique(compiled_clust_df_v3$upset[compiled_clust_df_v3$Clone_IDs %in% multi_sample_clones_v1])

length(brg_clones_revision)
#[1] 1715
length(trg_clones_revision)
#[1] 334

#install.packages("UpSetR")
#library(UpSetR)

upset(fromExpression(table(compiled_clust_df_v3$upset[compiled_clust_df_v3$upset %in% multi_sample_clone_types_v3 & compiled_clust_df_v3$Clone_IDs %in% brg_clones_revision])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dim(left_join(compiled_clust_df_v3, ls_synthesis_meta[,c('age_binned' ,'Clone_IDs')], multiple="first"))
compiled_clust_df_v3 = left_join(compiled_clust_df_v3, ls_synthesis_meta[,c('age_binned' ,'Clone_IDs')], multiple="first")

pdf('C:/Users/mattg/Desktop/250429_FINAL_FIGURES/ED4/upset_bRG_clones_tRG_as_RG.pdf', width=12, height=8)
upset(fromExpression(table(compiled_clust_df_v3$upset[compiled_clust_df_v3$upset %in% multi_sample_clone_types_v3 & compiled_clust_df_v3$Clone_IDs %in% brg_clones_revision])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()


## ED4a -- all cells
## make another compiled clust df with the zoomed out definitions of cell types...
## v4 -- broad_celltypes
compiled_clust_df_v4 = data.frame(Clone_IDs=multicell_clones)
compiled_clust_df_v4$EX = rep(0,length(multicell_clones))
compiled_clust_df_v4$IN = rep(0,length(multicell_clones))
compiled_clust_df_v4$RG = rep(0,length(multicell_clones))
compiled_clust_df_v4$OPCs = rep(0,length(multicell_clones))
compiled_clust_df_v4$Glia = rep(0,length(multicell_clones))
compiled_clust_df_v4$EX_pct = rep(0,length(multicell_clones))
compiled_clust_df_v4$IN_pct = rep(0,length(multicell_clones))
compiled_clust_df_v4$RG_pct = rep(0,length(multicell_clones))
compiled_clust_df_v4$OPCs_pct = rep(0,length(multicell_clones))
compiled_clust_df_v4$Glia_pct = rep(0,length(multicell_clones))

for(clone in multicell_clones) {
  clone_size = clone_df$Clone_size[clone_df$Clone_IDs %in% clone][1]
  clone_clusts = table(ls_synthesis@meta.data$broad_celltype[ls_synthesis@meta.data$Clone_IDs %in% clone])
  clone_row = which(compiled_clust_df_v4$Clone_IDs %in% clone)
  compiled_clust_df_v4[clone_row,"EX"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX"])
  compiled_clust_df_v4[clone_row,"IN"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN"])
  compiled_clust_df_v4[clone_row,"RG"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "RG"])
  compiled_clust_df_v4[clone_row,"OPCs"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])
  compiled_clust_df_v4[clone_row,"Glia"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Glia"])
  compiled_clust_df_v4[clone_row,"EX_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "EX"])/clone_size
  compiled_clust_df_v4[clone_row,"IN_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "IN"])/clone_size
  compiled_clust_df_v4[clone_row,"RG_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "RG"])/clone_size
  compiled_clust_df_v4[clone_row,"OPCs_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "OPCs"])/clone_size
  compiled_clust_df_v4[clone_row,"Glia_pct"] = sum(clone_clusts[dimnames(clone_clusts)[[1]] %in% "Glia"])/clone_size
}
head(compiled_clust_df_v4)
dim(compiled_clust_df_v4)
#[1] 6402   11

upset_notations = c()
for(r in 1:nrow(compiled_clust_df_v4)) {
  upset_notations=c(upset_notations, paste(colnames(compiled_clust_df_v4[r,2:6] %>% select(where(~ any(. != 0)))), collapse="&"))
}
upset_notations
compiled_clust_df_v4$upset = upset_notations
table(compiled_clust_df_v4$upset)

dim(left_join(compiled_clust_df_v4, ls_synthesis_meta[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first"))
#[1] 6402   16
compiled_clust_df_v4 = left_join(compiled_clust_df_v4, ls_synthesis_meta[,c('sample_id' ,'Clone_IDs', 'sample_age', 'GZ_origin', 'Clone_size')], multiple="first")
head(compiled_clust_df_v4)

# drop clone types only in one sample
all_clone_types_v4 = unique(compiled_clust_df_v4$upset)
all_clone_types_sample_counts_v4 = data.frame("clone_types" = all_clone_types_v4, "sample_counts" = 0)
all_clone_types_sample_counts_v4
for (ct in all_clone_types_v4){
  all_clone_types_sample_counts_v4$sample_counts[all_clone_types_sample_counts_v4$clone_types %in% ct] = length(table(compiled_clust_df_v4$sample_id[compiled_clust_df_v4$upset %in% ct]))
}
table(all_clone_types_sample_counts_v4$sample_counts)
# 1 4 5 6 7 8 
# 2 4 4 2 9 5
length(all_clone_types_sample_counts_v4$clone_types[all_clone_types_sample_counts_v4$sample_counts <= 1])
#[1] 2
length(all_clone_types_sample_counts_v4$clone_types[all_clone_types_sample_counts_v4$sample_counts > 1])
#[1] 24
single_sample_clone_types_v4 = all_clone_types_sample_counts_v4$clone_types[all_clone_types_sample_counts_v4$sample_counts <= 1]
multi_sample_clone_types_v4 = all_clone_types_sample_counts_v4$clone_types[all_clone_types_sample_counts_v4$sample_counts > 1]
single_sample_clone_type_clone_ids_v4 = compiled_clust_df_v4$Clone_IDs[compiled_clust_df_v4$upset %in% single_sample_clone_types_v4]

upset(fromExpression(table(compiled_clust_df_v4$upset[compiled_clust_df_v4$upset %in% multi_sample_clone_types_v4])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))

compiled_clust_df_v4 = left_join(compiled_clust_df_v4, ls_synthesis_meta[,c('Clone_IDs', 'age_binned')], multiple="first")

pdf(paste(savepath,'/upset_plot_all_celltypes_all.pdf', sep=""), width=12, height=8)
focus_age = c("<=GW20", ">GW20")
upset_clones = compiled_clust_df_v4$Clone_IDs[compiled_clust_df_v4$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v4$upset[compiled_clust_df_v4$upset %in% multi_sample_clone_types_v4 & compiled_clust_df_v4$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_all_celltypes_young.pdf', sep=""), width=12, height=8)
focus_age = "<=GW20"
upset_clones = compiled_clust_df_v4$Clone_IDs[compiled_clust_df_v4$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v4$upset[compiled_clust_df_v4$upset %in% multi_sample_clone_types_v4 & compiled_clust_df_v4$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_all_celltypes_old.pdf', sep=""), width=12, height=8)
focus_age = ">GW20"
upset_clones = compiled_clust_df_v4$Clone_IDs[compiled_clust_df_v4$age_binned %in% focus_age] 
upset(fromExpression(table(compiled_clust_df_v4$upset[compiled_clust_df_v4$upset %in% multi_sample_clone_types_v4 & compiled_clust_df_v4$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

## Figure ED5 -- ls_in
## Fig ED5a -- marker gene expression in full UMAP to subset down to ls_in
FeaturePlot(ls_synthesis, features=c("DLX2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("DLX2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_DLX2.png", image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("MKI67")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_MKI67.png", image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("GAD2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("GAD2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_GAD2.png", image, width=4, height=3.2)
FeaturePlot(ls_synthesis, features=c("OLIG2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
image = FeaturePlot(ls_synthesis, features=c("OLIG2")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"))
ggsave("ls_full_umap_OLIG2.png", image, width=4, height=3.2)

DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[rownames(ls_synthesis@meta.data) %in% rownames(ls_synthesis_in@meta.data)]) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
image = DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[rownames(ls_synthesis@meta.data) %in% rownames(ls_synthesis_in@meta.data)]) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30, family="Helvetica"), legend.position = "none")
ggsave("ls_full_umap_ls_synthesis_in_highlighted.png", image, width=4, height=3.2)

## Fig ED5b -- ls_in by seurat cluster
DimPlot(ls_synthesis_in, group.by="seurat_clusters", label=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=20))
image = DimPlot(ls_synthesis_in, group.by="seurat_clusters", label=T) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=20))
ggsave("ls_synthesis_in_umap_by_seurat_cluster.png", width=10, height=8)
ggsave("ls_synthesis_in_umap_by_seurat_cluster.svg", width=10, height=8)

## Fig ED5c -- ls_in by sample
DimPlot(ls_synthesis_in, pt.size=1, group.by="sample_id") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30), plot.title = element_blank())
image = DimPlot(ls_synthesis_in, pt.size=1, group.by="sample_id") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size=30), plot.title = element_blank())
ggsave("ls_synthesis_in_umap_grouped_sample.png", image, width=10, height=8)
ggsave("ls_synthesis_in_umap_grouped_sample.svg", image, width=10, height=8)


## Fig ED5d -- marker gene expression on ls_in UMAP
FeaturePlot(ls_synthesis_in, features=c("GAD2"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("GAD2"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_gad2.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_gad2.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("NPY"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("NPY"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_npy.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_npy.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("NR2F2"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("NR2F2"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_nr2f2.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_nr2f2.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("LHX6"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("LHX6"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_lhx6.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_lhx6.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("SOX6"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("SOX6"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_sox6.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_sox6.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("PAX6"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("PAX6"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_pax6.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_pax6.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("ERBB4"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("ERBB4"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_erbb4.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_erbb4.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("PBX3"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("PBX3"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_pbx3.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_pbx3.png", plot=image, width=4, height=3.2)

FeaturePlot(ls_synthesis_in, features=c("SCGN"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
image = FeaturePlot(ls_synthesis_in, features=c("SCGN"),raster=F, cols=c("#e6e8e7", "#c44e33")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size=30))
ggsave(file="ls_synthesis_in_umap_scgn.svg", plot=image, width=4, height=3.2)
ggsave(file="ls_synthesis_in_umap_scgn.png", plot=image, width=4, height=3.2)

## Fig ED5e, ED5f -- IN cells with shared dorsal barcodes ##
dorsal_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$ENs>=1 | compiled_clust_df_v1$EX_IPCs>=1 | compiled_clust_df_v1$RG>=1]
length(dorsal_clones)
# [1] 4802
dorsal_in_shared_clones = dorsal_clones[dorsal_clones %in% in_clones]
length(dorsal_in_shared_clones)
# [1] 819
DimPlot(ls_synthesis_in, group.by="seurat_clusters")

DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$Clone_IDs %in% dorsal_in_shared_clones & ls_synthesis@meta.data$Clone_IDs %in% multi_sample_clones_v1])

DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$Clone_IDs %in% dorsal_in_shared_clones & ls_synthesis@meta.data$Clone_IDs %in% multi_sample_clones_v1]) +
  ggtitle('Dorsal and IN shared clones') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.title = element_text(size=25, hjust=0.5))
image = DimPlot(ls_synthesis, cells.highlight = rownames(ls_synthesis@meta.data)[ls_synthesis@meta.data$Clone_IDs %in% dorsal_in_shared_clones & ls_synthesis@meta.data$Clone_IDs %in% multi_sample_clones_v1]) +
  ggtitle('Dorsal and IN shared clones') +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "none", plot.title = element_text(size=25, hjust=0.5))
ggsave(file="full_umap_dorsal_in_shared_clones_2cells_3umis_fixed.png", plot=image, width=12, height=8)
ggsave(file="full_umap_dorsal_in_shared_clones_2cells_3umis_fixed.svg", plot=image, width=12, height=8)

ggplot(subset(ls_synthesis_in@meta.data, Clone_IDs %in% dorsal_in_shared_clones & ls_synthesis_in@meta.data$Clone_IDs %in% multi_sample_clones_v1), aes(x=seurat_clusters)) +
  geom_bar() +
  scale_x_discrete(drop=F) +
  ylab("Number of cells") +
  xlab("Subcluster") +
  ggtitle("Cells with shared dorsal barcodes") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(subset(ls_synthesis_in@meta.data, Clone_IDs %in% dorsal_in_shared_clones & ls_synthesis_in@meta.data$Clone_IDs %in% multi_sample_clones_v1), aes(x=seurat_clusters)) +
  geom_bar() +
  scale_x_discrete(drop=F) +
  ylab("Number of cells") +
  xlab("Subcluster") +
  ggtitle("Cells with shared dorsal barcodes") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("ls_in_v2_barchart_cells_with_dorsal_barcodes.svg", image, width=10, height=8)


## Fig ED5h -- dotplot of marker genes from Delgado STICR paper for ls_in
ls_synthesis_in = SetIdent(ls_synthesis_in, value="seurat_clusters")
ls_synthesis_in@active.ident = factor(ls_synthesis_in@active.ident, levels=c(5,3,2,4,9,7))
DotPlot(subset(ls_synthesis_in, seurat_clusters %in% c(2,3,4,5,7,9)), features=c("COL1A2", "LGALS1", "NPY",
                                                                                 "MEG3", "KLHL35", "NR2F1", "NFIX", "SOX6", "PROX1",
                                                                                 "MEIS2", "RUNX1T1", "ETV1", "PBX3", "TSHZ1")) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = DotPlot(subset(ls_synthesis_in, seurat_clusters %in% c(2,3,4,5,7,9)), features=c("COL1A2", "LGALS1", "NPY",
                                                                                         "MEG3", "KLHL35", "NR2F1", "NFIX", "SOX6", "PROX1",
                                                                                         "MEIS2", "RUNX1T1", "ETV1", "PBX3", "TSHZ1")) +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=30, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("ls_synthesis_in_dotplot_local_clusters_only.svg", image, width=10, height=8)

## Fig ED5g -- heatmap by seurat cluster for ls_in genes
ls_synthesis_in = SetIdent(ls_synthesis_in, value="seurat_clusters")
ls_synthesis_in.seurat_cluster.markers <- FindAllMarkers(ls_synthesis_in, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(ls_synthesis_in.seurat_cluster.markers)
ls_synthesis_in_seurat_clustermarkers_unfilt = as.data.frame(ls_synthesis_in.seurat_cluster.markers %>% group_by(cluster))
write.csv(ls_synthesis_in_seurat_clustermarkers_unfilt,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_synthesis_in_clustermarkers_seurat_clusters_v1_unfilt.csv',row.names=F)
ls_synthesis_in.seurat_cluster.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
ls_synthesis_in_seurat_clustermarkers = as.data.frame(ls_synthesis_in.seurat_cluster.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
write.csv(ls_synthesis_in_seurat_clustermarkers,'/media/chang/HDD-11/mgkeefe/R_work_in_progress_cajal/ls_synthesis_in_clustermarkers_seurat_clusters_v1.csv',row.names=F)

ls_synthesis_in_seurat_clustermarkers_top10_clustermarkers = as.data.frame(ls_synthesis_in.seurat_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
DoHeatmap(ls_synthesis_in, features=ls_synthesis_in_seurat_clustermarkers_top10_clustermarkers$gene)
image = DoHeatmap(ls_synthesis_in, features=ls_synthesis_in_seurat_clustermarkers_top10_clustermarkers$gene)
ggsave("heatmap_ls_synthesis_in_by_seurat_cluster.png", width=8, height=12)
ggsave("heatmap_ls_synthesis_in_by_seurat_cluster.svg", width=8, height=12)


#### Fig ED9
## Fig ED9a-b -- Upset plots for GW20 PFC and GW20 V1
head(compiled_clust_df_v1)
table(compiled_clust_df_v1$sample_id)

focus_ident = c("GW20_PFC_L1")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))

focus_ident = c("GW20_V1_L1", "GW20_V1_L2")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))

pdf(paste(savepath,'/upset_plot_PFC.pdf', sep=""), width=12, height=8)
focus_ident = c("GW20_PFC_L1")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_V1.pdf', sep=""), width=12, height=8)
focus_ident = c("GW20_V1_L1", "GW20_V1_L2")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

#### Fig ED9
## Fig ED9a-b -- Upset plots for GW20 PFC and GW20 V1
head(compiled_clust_df_v1)
table(compiled_clust_df_v1$sample_id)

focus_ident = c("GW20_PFC")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))

focus_ident = c("GW20_V1")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))

pdf(paste(savepath,'/upset_plot_PFC.pdf', sep=""), width=12, height=8)
focus_ident = c("GW20_PFC")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

pdf(paste(savepath,'/upset_plot_V1.pdf', sep=""), width=12, height=8)
focus_ident = c("GW20_V1")
upset_clones = compiled_clust_df_v1$Clone_IDs[compiled_clust_df_v1$sample_id %in% focus_ident] 
upset(fromExpression(table(compiled_clust_df_v1$upset[compiled_clust_df_v1$Clone_IDs %in% upset_clones])),
      nsets=13,
      nintersects=30,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

## Fig ED10 -- clone randomization (REVISION)
## Fig ED10b -- unique cell types per clone
compiled_clust_df_v1$number_unique_celltypes <- apply(compiled_clust_df_v1[, 2:12], 1, function(x) sum(x > 0))
df_aggregated <- compiled_clust_df_v1 %>%
  group_by(Clone_size, number_unique_celltypes) %>%
  summarise(count = n(), .groups = 'drop')
df_aggregated
df_aggregated$log_count = log10(df_aggregated$count)
# Create dot plot
ggplot(df_aggregated, aes(x = Clone_size, y = number_unique_celltypes, size = log_count)) +
  geom_point(color = "blue", alpha = 0.6) +  # Add color and transparency to points
  labs(title = "Number of unique cell types per clone, by clone size",
       x = "Clone size",
       y = "Number of unique cell types in clone",
       size = "Number of clones (log10)") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size=30, family="Helvetica"),
        axis.text = element_text(size=15, family="Helvetica"),
        axis.title=element_text(size=15, family="Helvetica")) +
  scale_size_continuous(range = c(2, 8)) +  # Control dot size range
  xlim(0,26)
image = ggplot(df_aggregated, aes(x = Clone_size, y = number_unique_celltypes, size = log_count)) +
  geom_point(color = "blue", alpha = 0.6) +  # Add color and transparency to points
  labs(title = "Number of unique cell types per clone, by clone size",
       x = "Clone size",
       y = "Number of unique cell types in clone",
       size = "Number of clones (log10)") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size=30, family="Helvetica"),
        axis.text = element_text(size=15, family="Helvetica"),
        axis.title=element_text(size=15, family="Helvetica")) +
  scale_size_continuous(range = c(2, 8)) +  # Control dot size range
  xlim(0,26)
ggsave("unique_celltypes_per_clone_by_size_all_clones.svg", image, width=10, height=8)


## Fig ED10c -- Randomized unique cell types per clone
# run on full dataset:
total_cell_counts <- compiled_clust_df_v1[,1:12] %>%
  select(-Clone_IDs) %>%
  summarise(across(everything(), sum))
total_cell_counts
# ENs EX_IPCs IN_local IN_CGE IN_MGE IN_OB IN_IPCs Astrocytes OPCs Oligodendrocytes   RG
# 2855    5237     2057     12     26    11    2915        175 1919              192 3431

create_random_clone = function(cell_pool, clone_size) {
  available_pool = unlist(lapply(names(cell_pool), function(cell_type) rep(cell_type, cell_pool[[cell_type]])))
  new_clone = sample(available_pool, size = clone_size, replace = FALSE)
  clone_distribution = table(new_clone)
  full_clone_distribution = as.integer(clone_distribution[names(total_cell_counts)])
  full_clone_distribution[is.na(full_clone_distribution)] = 0
  updated_pool = cell_pool-full_clone_distribution
  return(list(distribution = full_clone_distribution, updated_pool = updated_pool))
}

cell_pool = total_cell_counts
new_data = list()
for (i in 1:nrow(compiled_clust_df_v1[,1:12])) {
  clone_size = sum(compiled_clust_df_v1[,1:12][i, -1])  # Total size of the current clone
  result = create_random_clone(cell_pool, clone_size)
  #print(result)
  cell_pool = result$updated_pool
  new_row = c(Clone_IDs = compiled_clust_df_v1$Clone_IDs[i], result$distribution)
  new_data[[i]]=new_row}
new_data
compiled_clust_df_v1_randomized <- do.call(rbind, new_data)
compiled_clust_df_v1_randomized <- as.data.frame(compiled_clust_df_v1_randomized)
colnames(compiled_clust_df_v1_randomized) <- c("Clone_IDs", names(total_cell_counts))
compiled_clust_df_v1_randomized
head(compiled_clust_df_v1_randomized)
sum(compiled_clust_df_v1$ENs)

head(compiled_clust_df_v1_randomized[,2:12])
compiled_clust_df_v1_randomized[,2:12] <- sapply(compiled_clust_df_v1_randomized[,2:12],as.numeric)
compiled_clust_df_v1_randomized <- compiled_clust_df_v1_randomized %>%
  mutate(Clone_size = rowSums(across(ENs:RG)))
table(compiled_clust_df_v1$Clone_size)
table(compiled_clust_df_v1_randomized$clone_size)

compiled_clust_df_v1_randomized = compiled_clust_df_v1_randomized %>%
  left_join(compiled_clust_df_v1 %>% select(Clone_IDs, sample_id, sample_age, GZ_origin, age_binned), 
            by = "Clone_IDs")
head(compiled_clust_df_v1_randomized)
compiled_clust_df_v1_randomized$number_unique_celltypes <- apply(compiled_clust_df_v1_randomized[, 2:12], 1, function(x) sum(x > 0))

df_aggregated <- compiled_clust_df_v1_randomized %>%
  group_by(Clone_size, number_unique_celltypes) %>%
  summarise(count = n(), .groups = 'drop')
df_aggregated
df_aggregated$log_count = log10(df_aggregated$count)
ggplot(df_aggregated, aes(x = Clone_size, y = number_unique_celltypes, size = log_count)) +
  geom_point(color = "red", alpha = 0.6) +  # Add color and transparency to points
  labs(title = "Number of unique cell types per clone -- RANDOMIZED",
       x = "Clone size",
       y = "Number of unique cell types in clone",
       size = "Number of clones (log10)") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size=30, family="Helvetica"),
        axis.text = element_text(size=15, family="Helvetica"),
        axis.title=element_text(size=15, family="Helvetica")) +
  scale_size_continuous(range = c(2, 8)) +  # Control dot size range
  xlim(0,26)
image = ggplot(df_aggregated, aes(x = Clone_size, y = number_unique_celltypes, size = log_count)) +
  geom_point(color = "red", alpha = 0.6) +  # Add color and transparency to points
  labs(title = "Number of unique cell types per clone -- RANDOMIZED",
       x = "Clone size",
       y = "Number of unique cell types in clone",
       size = "Number of clones (log10)") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size=30, family="Helvetica"),
        axis.text = element_text(size=15, family="Helvetica"),
        axis.title=element_text(size=15, family="Helvetica")) +
  scale_size_continuous(range = c(2, 8)) +  # Control dot size range
  xlim(0,26)
ggsave("unique_celltypes_per_clone_by_size_all_RANDOM_clones.svg", width=10, height=8)


## Fig ED10d -- Distribution of unique cell types per clone
get_random_clone_unique_celltypes = function(cell_pool, clone_size) {
  available_pool = unlist(lapply(names(cell_pool), function(cell_type) rep(cell_type, cell_pool[[cell_type]])))
  new_clone = sample(available_pool, size = clone_size, replace = FALSE)
  new_clone_unique_celltypes = length(unique(new_clone))
  clone_distribution = table(new_clone)
  full_clone_distribution = as.integer(clone_distribution[names(total_cell_counts)])
  full_clone_distribution[is.na(full_clone_distribution)] = 0
  updated_pool = cell_pool-full_clone_distribution
  return(list(unique_celltypes = new_clone_unique_celltypes, updated_pool = updated_pool))
}
#make new dataframe to work on
randomized_df = compiled_clust_df_v1 %>% select(Clone_IDs, Clone_size, number_unique_celltypes)
head(randomized_df)
for (i in 1:1000){
  set.seed(i*3)
  new_colname = paste("i",i,sep="")
  # Step 1: Create the character vector (the pool)
  pool = unlist(lapply(names(total_cell_counts), function(cell_type) rep(cell_type, total_cell_counts[[cell_type]])))
  # Step 2: Extract clone sizes from the existing dataframe
  num_clones <- nrow(randomized_df)
  clone_sizes <- randomized_df$Clone_size  # Use the existing column
  # Step 4: Shuffle pool to ensure random selection
  pool <- sample(pool)
  # Step 5: Compute cumulative indices for slicing
  clone_indices <- cumsum(clone_sizes)  # End indices for each clone
  start_indices <- c(1, head(clone_indices, -1) + 1)  # Start indices
  # Step 6: Generate new clones by slicing the shuffled pool
  sampled_clones <- mapply(function(start, end) pool[start:end], start_indices, clone_indices, SIMPLIFY = FALSE)
  # Step 7: Get the number of unique celltypes in each clone and add to the dataframe
  randomized_df[dim(randomized_df)[2]+1] <- sapply(sampled_clones, function(x) length(unique(x)))
  colnames(randomized_df)[dim(randomized_df)[2]] = new_colname
}
head(randomized_df)


all_clones_fdr_df = data.frame(Var1 = 1:8)
all_clones_fdr_df$Var1 = factor(all_clones_fdr_df$Var1)
for (i in c(1:1000)){
  colname=paste("i",i,sep="")
  all_clones_fdr_df = left_join(all_clones_fdr_df, data.frame(table(randomized_df[,i+3])))
  colnames(all_clones_fdr_df)[i+1] = colname
}
all_clones_fdr_df[1:5, 1:5]
all_clones_fdr_df = left_join(all_clones_fdr_df, data.frame(table(randomized_df$number_unique_celltypes)))
colnames(all_clones_fdr_df)[colnames(all_clones_fdr_df) == "Freq"] = "observed_unique_celltypes"
colnames(all_clones_fdr_df)[colnames(all_clones_fdr_df) == "Var1"] = "unique_celltypes_per_clone"
all_clones_fdr_df

## p-value is the fraction of random iterations that have counts greater than or equal to the observed count
all_clones_fdr_df$observed_unique_celltypes[7] = 0
all_clones_fdr_df$observed_unique_celltypes[8] = 0
sum(all_clones_fdr_df[1, 2:1001] >= all_clones_fdr_df$observed_unique_celltypes[1])/1000

for (r in c(1:8)){
  thresh = all_clones_fdr_df$observed_unique_celltypes[r]
  all_clones_fdr_df$p_val[r] = sum(all_clones_fdr_df[r, 2:1001] >= thresh)/1000
}
all_clones_fdr_df[,1001:1003]

all_clones_fdr_df[is.na(all_clones_fdr_df)] <- 0

# Create a new dataframe for visualization
library(tidyr)
obs_dist <- all_clones_fdr_df %>%
  select(unique_celltypes_per_clone, observed_unique_celltypes) %>%
  mutate(source = "Observed")

boot_dist <- all_clones_fdr_df %>%
  select(unique_celltypes_per_clone, starts_with("i")) %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "count") %>%
  mutate(source = "Bootstrapped")

# Combine observed and bootstrapped for side-by-side visualization
combined_dist <- bind_rows(obs_dist, boot_dist)

# Plot the distributions
ggplot(combined_dist, aes(x = unique_celltypes_per_clone, y = count, fill = source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = obs_dist, aes(x = unique_celltypes_per_clone, y = observed_unique_celltypes), color = "red", size = 2) +  # Mark observed values
  labs(title = "All clones -- Observed vs. bootstrapped # of unique cells per clone", 
       x = "Number of Unique Cell Types", y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "gray"))

# get CIs
bootstrap_results <- all_clones_fdr_df %>%
  select(unique_celltypes_per_clone, starts_with("i")) %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "count") %>%
  group_by(unique_celltypes_per_clone) %>%
  summarise(
    lower_95CI = quantile(count, 0.025),  # Lower 2.5 percentile
    upper_95CI = quantile(count, 0.975),  # Upper 97.5 percentile
    .groups = 'drop'
  )
bootstrap_results

# Calculate p-values based on whether the observed value falls outside the 95% CI
all_clones_fdr_df <- all_clones_fdr_df %>%
  left_join(bootstrap_results, by = "unique_celltypes_per_clone")
all_clones_fdr_df[,1000:1005]

for (r in c(1:8)){all_clones_fdr_df$mean[r] = mean(as.numeric(all_clones_fdr_df[r, 2:1001]))}
all_clones_fdr_df$mean

ggplot(all_clones_fdr_df, aes(x = unique_celltypes_per_clone)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "darkgray",
    size = 1
  ) +
  # Plot the observed value as a red point
  geom_point(aes(y = observed_unique_celltypes), color = "blue", size = 3) +
  # Add labels and theme
  labs(
    title = "Distribution of unique cell types per clones",
    x = "Number of unique cell types",
    y = "Number of clones",
    subtitle = "Blue: Mean of random distribution with 95% CI | Red: Observed Value"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(all_clones_fdr_df, aes(x = unique_celltypes_per_clone)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "darkgray",
    size = 1
  ) +
  # Plot the observed value as a red point
  geom_point(aes(y = observed_unique_celltypes), color = "blue", size = 3) +
  # Add labels and theme
  labs(
    title = "Distribution of unique cell types per clones",
    x = "Number of unique cell types",
    y = "Number of clones",
    subtitle = "Blue: Mean of random distribution with 95% CI | Red: Observed Value"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("all_clones_unique_celltypes_dist_with_95CI.svg", image, width=10, height=8)

## Fig ED10e -- Clone composition, tRG/bRG clones only
## definitely need to bring in some more code before this to get the rg_clone_comparison_df, etc.
## need to also figure out what the difference is between filt1 and not -- seems like the figure uses filt1, but can't tell why
compiled_clust_df_rg = compiled_clust_df_v1[compiled_clust_df_v1$Clone_IDs %in% c(brg_clones_revision, trg_clones_revision),]
dim(compiled_clust_df_rg)
compiled_clust_df_rg$rg_clone_type[compiled_clust_df_rg$Clone_IDs %in% c(brg_clones_revision)] = "bRG"
compiled_clust_df_rg$rg_clone_type[compiled_clust_df_rg$Clone_IDs %in% c(trg_clones_revision)] = "tRG"
table(compiled_clust_df_rg$rg_clone_type)

total_cell_counts <- compiled_clust_df_rg[,1:12] %>%
  select(-Clone_IDs) %>%
  summarise(across(everything(), sum))
total_cell_counts
# ENs EX_IPCs IN_local IN_CGE IN_MGE IN_OB IN_IPCs Astrocytes OPCs Oligodendrocytes  RG
#   45     247       31      0      1     0      79         51  124                5 998
cell_pool = total_cell_counts
new_data = list()
for (i in 1:nrow(compiled_clust_df_rg[,1:12])) {
  clone_size = sum(compiled_clust_df_rg[,1:12][i, -1])  # Total size of the current clone
  result = create_random_clone(cell_pool, clone_size)
  #print(result)
  cell_pool = result$updated_pool
  new_row = c(Clone_IDs = compiled_clust_df_rg$Clone_IDs[i], result$distribution)
  new_data[[i]]=new_row}
#new_data
compiled_clust_df_rg_randomized <- do.call(rbind, new_data)
compiled_clust_df_rg_randomized <- as.data.frame(compiled_clust_df_rg_randomized)
colnames(compiled_clust_df_rg_randomized) <- c("Clone_IDs", names(total_cell_counts))
compiled_clust_df_rg_randomized
head(compiled_clust_df_rg_randomized)
head(compiled_clust_df_rg[1:12])
compiled_clust_df_rg_randomized[2:12] <- lapply(compiled_clust_df_rg_randomized[2:12], as.integer)
sum(compiled_clust_df_rg$EX_IPCs)
sum(compiled_clust_df_rg_randomized$EX_IPCs)
compiled_clust_df_rg_randomized %>%
  mutate(clone_size = rowSums(across(ENs:RG)))
compiled_clust_df_rg_randomized <- compiled_clust_df_rg_randomized %>%
  mutate(clone_size = rowSums(across(ENs:RG)))

compiled_clust_df_rg_randomized = compiled_clust_df_rg_randomized %>%
  left_join(compiled_clust_df_rg %>% select(Clone_IDs, sample_id, sample_age, GZ_origin, age_binned), 
            by = "Clone_IDs")
compiled_clust_df_rg_randomized$number_unique_celltypes <- apply(compiled_clust_df_rg_randomized[, 2:12], 1, function(x) sum(x > 0))
head(compiled_clust_df_rg_randomized)

upset_notations = c()
for(r in 1:nrow(compiled_clust_df_rg_randomized)) {
  upset_notations=c(upset_notations, paste(colnames(compiled_clust_df_rg_randomized[r,2:12] %>% select(where(~ any(. != 0)))), collapse="&"))
}
head(upset_notations)
compiled_clust_df_rg_randomized$upset = upset_notations
table(compiled_clust_df_rg_randomized$upset)

upset(fromExpression(table(compiled_clust_df_rg_randomized$upset)),
      nsets=13,
      nintersects=40,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
upset(fromExpression(table(compiled_clust_df_rg$upset)),
      nsets=13,
      nintersects=40,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))

## bootstrap
randomized_df_rg = compiled_clust_df_rg %>% select(Clone_IDs, Clone_size, number_unique_celltypes, upset)
total_cell_counts <- compiled_clust_df_rg[,1:12] %>%
  select(-Clone_IDs) %>%
  summarise(across(everything(), sum))
cell_pool = total_cell_counts
new_data = list()

for (i in 1:1000){
  set.seed(i*3)
  new_colname = paste("i",i,sep="")
  cell_pool = total_cell_counts
  new_data = list()
  for (r in 1:nrow(compiled_clust_df_rg[,1:12])) {
    clone_size = sum(compiled_clust_df_rg[,1:12][r, -1])  # Total size of the current clone
    result = create_random_clone(cell_pool, clone_size)
    #print(result)
    cell_pool = result$updated_pool
    new_row = c(Clone_IDs = compiled_clust_df_rg$Clone_IDs[r], result$distribution)
    new_data[[r]]=new_row}
  i_random_df = do.call(rbind, new_data)
  i_random_df = as.data.frame(i_random_df)
  colnames(i_random_df) <- c("Clone_IDs", names(total_cell_counts))
  i_random_df[2:12] <- lapply(i_random_df[2:12], as.integer)
  upset_notations = c()
  for(j in 1:nrow(i_random_df)) {
    upset_notations=c(upset_notations, paste(colnames(i_random_df[j,2:12] %>% select(where(~ any(. != 0)))), collapse="&"))
  }
  randomized_df_rg[dim(randomized_df_rg)[2]+1] = upset_notations
  colnames(randomized_df_rg)[dim(randomized_df_rg)[2]] = new_colname
  print(i)
}
head(randomized_df_rg)

# Convert bootstraps to long format
boot_long <- randomized_df_rg %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "clone_type")

# Count occurrences of each clone type across all iterations
boot_counts <- boot_long %>%
  group_by(clone_type) %>%
  summarise(mean_count = n() / length(unique(iteration)), .groups = "drop")


actual_rg_clone_comps <- as.data.frame(table(randomized_df_rg$upset))
colnames(actual_rg_clone_comps) <- c("clone_type", "actual_count")

# Merge with bootstrapped mean counts
rg_clone_comparison_df <- actual_rg_clone_comps %>%
  right_join(boot_counts, by = "clone_type") %>%
  replace_na(list(mean_count = 0))  # Replace NA with 0 if a clone type was not observed in bootstraps

rg_clone_comparison_df

#keep just the "real" clone types instead?
rg_clone_comparison_df_filt = actual_rg_clone_comps %>%
  left_join(boot_counts, by = "clone_type") %>%
  replace_na(list(mean_count = 0))  # Replace NA with 0 if a clone type was not observed in bootstraps
rg_clone_comparison_df_filt


## everything done above was not quite right
## start from line 15963 from the full code block
rg_clone_comparison_df
randomized_df_rg[1:5,1000:1004]
randomized_df_rg$trg_brg_clone[randomized_df_rg$Clone_IDs %in% trg_clones_revision] = "tRG"
randomized_df_rg$trg_brg_clone[randomized_df_rg$Clone_IDs %in% brg_clones_revision] = "bRG"
dim(randomized_df_rg)
#[1]  572 1006
table(randomized_df_rg$trg_brg_clone)
#bRG tRG 
#468 104

dim(randomized_df_rg)
randomized_df_rg[1:5,1:5]

#compress DF and get 95% CIs for each clone type
rg_fdr_df = data.frame(Var1=unique(randomized_df_rg$upset))
rg_fdr_df$Var1 = factor(rg_fdr_df$Var1)
for (i in c(1:1000)){
  colname=paste("i",i,sep="")
  rg_fdr_df = left_join(rg_fdr_df, data.frame(table(randomized_df_rg[,i+4])))
  colnames(rg_fdr_df)[i+1] = colname
}
rg_fdr_df[1:5,1:10]
rg_fdr_df = left_join(rg_fdr_df, data.frame(table(randomized_df_rg$upset)))
colnames(rg_fdr_df)[colnames(rg_fdr_df) == "Freq"] = "observed_rg_clone_type"
colnames(rg_fdr_df)[colnames(rg_fdr_df) == "Var1"] = "rg_clone_type"
rg_fdr_df[, 1:5]
rg_fdr_df[is.na(rg_fdr_df)] = 0

bootstrap_rg <- rg_fdr_df %>%
  select(rg_clone_type, starts_with("i")) %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "count") %>%
  group_by(rg_clone_type) %>%
  summarise(
    lower_95CI = quantile(count, 0.025),  # Lower 2.5 percentile
    upper_95CI = quantile(count, 0.975),  # Upper 97.5 percentile
    .groups = 'drop'
  )
bootstrap_rg

for (r in c(1:length(bootstrap_rg$rg_clone_type))){rg_fdr_df$mean[r] = mean(as.numeric(rg_fdr_df[r, 2:1001]))}
rg_fdr_df$mean
rg_fdr_df[,c(1,1000:1003)]
rg_fdr_df <- rg_fdr_df %>%
  left_join(bootstrap_rg, by = "rg_clone_type")
rg_fdr_df[,c(1,1000:1005)]


rg_fdr_df$rg_clone_type <- factor(rg_fdr_df$rg_clone_type, levels = rg_fdr_df$rg_clone_type[order(rg_fdr_df$mean, decreasing = TRUE)])


ggplot(rg_fdr_df, aes(x = rg_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = observed_rg_clone_type), color = "orange", size = 3) +
  labs(
    title = "Clone composition, tRG/bRG clones only",
    y = "Number of clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=45, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.title.x = element_blank(),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))

image = ggplot(rg_fdr_df, aes(x = rg_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = observed_rg_clone_type), color = "orange", size = 3) +
  labs(
    title = "Clone composition, tRG/bRG clones only",
    y = "Number of clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=45, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.title.x = element_blank(),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("all_rg_clones_distribution_with_actuals.svg", image, width=10, height=8)

rg_fdr_df[1:5,1000:1005]
dim(randomized_df_rg)
rg_fdr_df$mean_ratio = rg_fdr_df$mean/dim(randomized_df_rg)[1]
rg_fdr_df$lower_95CI_ratio = rg_fdr_df$lower_95CI/dim(randomized_df_rg)[1]
rg_fdr_df$upper_95CI_ratio = rg_fdr_df$upper_95CI/dim(randomized_df_rg)[1]

randomized_df_rg[1:5,1000:1005]
trg_actual_clone_comp = data.frame(table(randomized_df_rg[randomized_df_rg$trg_brg_clone == "tRG", 'upset']))
brg_actual_clone_comp = data.frame(table(randomized_df_rg[randomized_df_rg$trg_brg_clone == "bRG", 'upset']))
colnames(trg_actual_clone_comp)[colnames(trg_actual_clone_comp) == "Freq"] = "trg_counts"
colnames(brg_actual_clone_comp)[colnames(brg_actual_clone_comp) == "Freq"] = "brg_counts"
colnames(trg_actual_clone_comp)[colnames(trg_actual_clone_comp) == "Var1"] = "rg_clone_type"
colnames(brg_actual_clone_comp)[colnames(brg_actual_clone_comp) == "Var1"] = "rg_clone_type"

actual_trg_brg_clone_comp = full_join(trg_actual_clone_comp, brg_actual_clone_comp)
actual_trg_brg_clone_comp[is.na(actual_trg_brg_clone_comp)] = 0
actual_trg_brg_clone_comp$trg_ratio = actual_trg_brg_clone_comp$trg_counts/sum(actual_trg_brg_clone_comp$trg_counts)
actual_trg_brg_clone_comp$brg_ratio = actual_trg_brg_clone_comp$brg_counts/sum(actual_trg_brg_clone_comp$brg_counts)
actual_trg_brg_clone_comp

rg_fdr_df = left_join(rg_fdr_df, actual_trg_brg_clone_comp, by="rg_clone_type")
rg_fdr_df[,1005:1012]


ggplot(rg_fdr_df, aes(x = rg_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean_ratio, color = "random"), size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI_ratio, ymax = upper_95CI_ratio),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = trg_ratio, color = "tRG"), size = 3) +
  geom_point(aes(y = brg_ratio, color = "bRG"), size = 3) +
  scale_color_manual(name="Clone type", values = c("random"="darkgray", "tRG"=trg_col, "bRG"=brg_col)) + 
  labs(
    title = "Clone composition, tRG/bRG clones only",
    y = "Ratio of all tRG/bRG clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=45, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.title.x = element_blank(),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))

upset_for_rg_fdr_df <- rg_fdr_df[,c("rg_clone_type", "mean")] %>%
  mutate(mean = ceiling(mean)) %>%  # Round mean to the nearest integer
  uncount(mean)  # Expand rows based on mean
upset_for_rg_fdr_df
pdf('randomized_mean_rg_upset_plot.pdf', width=14, height=8)
upset(fromExpression(table(upset_for_rg_fdr_df$rg_clone_type)),
      nsets=40,
      nintersects=40,
      order.by = "freq",
      decreasing = T,
      point.size=5, line.size=1.5,
      text.scale = c(4,3,3,3,3.5,0), mb.ratio = c(0.5, 0.5))
dev.off()

ggplot(rg_fdr_df[rg_fdr_df$trg_counts > 1 | rg_fdr_df$brg_counts > 1,], aes(x = rg_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean_ratio, color = "random"), size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI_ratio, ymax = upper_95CI_ratio),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = trg_ratio, color = "tRG"), size = 3) +
  geom_point(aes(y = brg_ratio, color = "bRG"), size = 3) +
  scale_color_manual(name="Clone type", values = c("random"="darkgray", "tRG"=trg_col, "bRG"=brg_col)) + 
  labs(
    title = "Clone composition, tRG/bRG clones only",
    y = "Ratio of all tRG/bRG clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=45, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.title.x = element_blank(),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(rg_fdr_df[rg_fdr_df$trg_counts > 1 | rg_fdr_df$brg_counts > 1,], aes(x = rg_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean_ratio, color = "random"), size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI_ratio, ymax = upper_95CI_ratio),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = trg_ratio, color = "tRG"), size = 3) +
  geom_point(aes(y = brg_ratio, color = "bRG"), size = 3) +
  scale_color_manual(name="Clone type", values = c("random"="darkgray", "tRG"=trg_col, "bRG"=brg_col)) + 
  labs(
    title = "Clone composition, tRG/bRG clones only",
    y = "Ratio of all tRG/bRG clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=45, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.title.x = element_blank(),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("all_rg_clones_distribution_with_actuals_trg_brg_separated_filt1.svg", image, width=10, height=8)


## Fig ED10f -- Two cell clone types, young samples
## going to do the randomization for 2-cell clones
compiled_clust_df_2cell = compiled_clust_df_v1[compiled_clust_df_v1$Clone_size == 2,]
# define the clones more simply, and closer to how it will be done for the randomization
p_types <- c("RG", "EX_IPCs", "IN_IPCs")
n_types <- c("ENs", "IN_local", "IN_OB", "IN_CGE", "IN_MGE")
g_types <- c("OPCs", "Astrocytes", "Oligodendrocytes")
# Create a new dataframe with P, N, and G counts
compiled_clust_df_2cell$P_count <- rowSums(compiled_clust_df_2cell[, p_types])
compiled_clust_df_2cell$N_count <- rowSums(compiled_clust_df_2cell[, n_types])
compiled_clust_df_2cell$G_count <- rowSums(compiled_clust_df_2cell[, g_types])
head(compiled_clust_df_2cell)
compiled_clust_df_2cell$two_cell_clone_type[compiled_clust_df_2cell$P_count==2] = "PP"
compiled_clust_df_2cell$two_cell_clone_type[compiled_clust_df_2cell$N_count==2] = "NN"
compiled_clust_df_2cell$two_cell_clone_type[compiled_clust_df_2cell$G_count==2] = "GG"
compiled_clust_df_2cell$two_cell_clone_type[compiled_clust_df_2cell$P_count==1 & compiled_clust_df_2cell$N_count==1] = "PN"
compiled_clust_df_2cell$two_cell_clone_type[compiled_clust_df_2cell$P_count==1 & compiled_clust_df_2cell$G_count==1] = "PG"
compiled_clust_df_2cell$two_cell_clone_type[compiled_clust_df_2cell$N_count==1 & compiled_clust_df_2cell$G_count==1] = "NG"
head(compiled_clust_df_2cell, 10)
table(compiled_clust_df_2cell$two_cell_clone_type)
#  GG   NG   NN   PG   PN   PP 
# 230   82  868  317  598 1605

compiled_clust_df_2cell_young = compiled_clust_df_2cell[compiled_clust_df_2cell$age_binned == "<=GW20",]
head(compiled_clust_df_2cell_young)
randomized_df_2cell_young = compiled_clust_df_2cell_young %>% select(Clone_IDs, Clone_size, number_unique_celltypes, upset, two_cell_clone_type)
total_cell_counts <- compiled_clust_df_v1[compiled_clust_df_v1$Clone_size == 2 & compiled_clust_df_v1$age_binned == "<=GW20",1:12] %>%
  select(-Clone_IDs) %>%
  summarise(across(everything(), sum))
cell_pool = total_cell_counts

total_cell_counts

for (i in 1:1000){
  # set.seed(i*3)
  new_colname = paste("i",i,sep="")
  # Step 1: Create the character vector (the pool)
  pool = unlist(lapply(names(total_cell_counts), function(cell_type) rep(cell_type, total_cell_counts[[cell_type]])))
  # Step 2: Simpllify the celltypes
  pool[pool %in% c("RG", "EX_IPCs", "IN_IPCs")] = "P"
  pool[pool %in% c("ENs", "IN_local", "IN_OB", "IN_CGE", "IN_MGE")] = "N"
  pool[pool %in% c("OPCs", "Astrocytes", "Oligodendrocytes")] = "G"
  # Step 4: Shuffle pool to ensure random selection
  pool <- sample(pool)
  # Step 5: Make new 2-cell clones
  new_clones <- paste0(pool[seq(1, length(pool) - 1, by = 2)], 
                       pool[seq(2, length(pool), by = 2)])
  clone_indices <- cumsum(clone_sizes)  # End indices for each clone
  start_indices <- c(1, head(clone_indices, -1) + 1)  # Start indices
  # Step 7: Add to the dataframe
  randomized_df_2cell_young[dim(randomized_df_2cell_young)[2]+1] <- new_clones
  colnames(randomized_df_2cell_young)[dim(randomized_df_2cell_young)[2]] = new_colname
}
head(randomized_df_2cell_young)
randomized_df_2cell_young[1:5,1:20]
write.csv(randomized_df_2cell_young, "randomized_1k_two_cell_clones_P_N_G_young.csv")

table(randomized_df_2cell_young$i104)
#  GG  GN  GP  NG  NN  NP  PG  PN  PP 
#   4  34  70  43 297 495  57 487 948

randomized_df_2cell_young = data.frame(lapply(randomized_df_2cell_young, function(x) {gsub("GN", "NG", x)}))
randomized_df_2cell_young = data.frame(lapply(randomized_df_2cell_young, function(x) {gsub("GP", "PG", x)}))
randomized_df_2cell_young = data.frame(lapply(randomized_df_2cell_young, function(x) {gsub("NP", "PN", x)}))
table(randomized_df_2cell_young$i104)
# GG  NG  NN  PG  PN  PP 
#  4  77 297 127 982 948 

# Reshape the dataframe to long format
df_long_young <- randomized_df_2cell_young %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "clone_type")
head(df_long_young)

# Count occurrences of each clone type in each iteration
df_abundance_young <- df_long_young %>%
  group_by(iteration, clone_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(iteration) %>%  # Ensure normalization happens per iteration
  mutate(relative_abundance = count / sum(count)) %>%
  ungroup()  # Remove grouping for further processing

head(df_abundance_young)

# Compute actual observed frequencies of `two_cell_clone_type`
table(randomized_df_2cell_young$two_cell_clone_type)
#GG   NG   NN   PG   PN   PP 
#43   21  601  105  430 1235 
table(randomized_df_2cell_young$two_cell_clone_type)/sum(table(randomized_df_2cell_young$two_cell_clone_type))
df_actual_young = data.frame(table(randomized_df_2cell_young$two_cell_clone_type)/sum(table(randomized_df_2cell_young$two_cell_clone_type)))
colnames(df_actual_young) = c("two_cell_clone_type", "relative_abundance")
df_actual_young
df_abundance_young$clone_type = factor(df_abundance_young$clone_type, levels=c("PP", "NN", "GG", "PN", "PG", "NG"))

# Create boxplot with actual values overlaid
ggplot(df_abundance_young, aes(x = clone_type, y = relative_abundance)) +
  geom_boxplot(outlier.shape = NA, fill = "gray", alpha = 0.5) +  # Boxplot
  geom_point(data = df_actual_young, 
             aes(x = two_cell_clone_type, y = relative_abundance, color = "≤GW20"), 
             size = 3) +
  scale_color_manual(name = "Sample age", 
                     values = c("≤GW20" = young_col)) +
  theme_minimal() +
  labs(x = "Clone type", y = "Relative abundance", title = "Distribution of two cell clone types <=GW20") +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))

## get 95% CIs
two_cell_young_fdr_df = data.frame(Var1=c("PP", "NN", "GG", "PN", "PG", "NG"))
two_cell_young_fdr_df$Var1 = factor(two_cell_young_fdr_df$Var1)
for (i in c(1:1000)){
  colname=paste("i",i,sep="")
  two_cell_young_fdr_df = left_join(two_cell_young_fdr_df, data.frame(table(randomized_df_2cell_young[,i+5])))
  colnames(two_cell_young_fdr_df)[i+1] = colname
}
two_cell_young_fdr_df[1:5,1:10]
two_cell_young_fdr_df = left_join(two_cell_young_fdr_df, data.frame(table(randomized_df_2cell_young$two_cell_clone_type)))
colnames(two_cell_young_fdr_df)[colnames(two_cell_young_fdr_df) == "Freq"] = "observed_two_cell_clone_type"
colnames(two_cell_young_fdr_df)[colnames(two_cell_young_fdr_df) == "Var1"] = "two_cell_clone_type"
two_cell_young_fdr_df[, 1:5]
two_cell_young_fdr_df[is.na(two_cell_young_fdr_df)] = 0

bootstrap_two_cell_young <- two_cell_young_fdr_df %>%
  select(two_cell_clone_type, starts_with("i")) %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "count") %>%
  group_by(two_cell_clone_type) %>%
  summarise(
    lower_95CI = quantile(count, 0.025),  # Lower 2.5 percentile
    upper_95CI = quantile(count, 0.975),  # Upper 97.5 percentile
    .groups = 'drop'
  )
bootstrap_two_cell_young

for (r in c(1:length(c("PP", "NN", "GG", "PN", "PG", "NG")))){two_cell_young_fdr_df$mean[r] = mean(as.numeric(two_cell_young_fdr_df[r, 2:1001]))}
two_cell_young_fdr_df$mean
two_cell_young_fdr_df[,c(1,1000:1003)]
two_cell_young_fdr_df <- two_cell_young_fdr_df %>%
  left_join(bootstrap_two_cell_young, by = "two_cell_clone_type")
two_cell_young_fdr_df[,c(1,1000:1005)]

two_cell_young_fdr_df$two_cell_clone_type = factor(two_cell_young_fdr_df$two_cell_clone_type, levels=c("PP", "NN", "GG", "PN", "PG", "NG"))
ggplot(two_cell_young_fdr_df, aes(x = two_cell_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = observed_two_cell_clone_type), color = young_col, size = 3) +
  labs(
    title = "Two cell clone types, ≤GW20 samples",
    x = "Number of unique cell types",
    y = "Number of clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(two_cell_young_fdr_df, aes(x = two_cell_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = observed_two_cell_clone_type), color = young_col, size = 3) +
  labs(
    title = "Two cell clone types, ≤GW20 samples",
    x = "Number of unique cell types",
    y = "Number of clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("randomized_vs_actual_two_cell_clone_type_dist_young.svg", image, width=10, height=8)

observed_values_two_cell_young = two_cell_young_fdr_df$observed_two_cell_clone_type
p_values_two_sided_two_cell_young <- sapply(1:length(observed_values_two_cell_young), function(i) {
  (sum(abs(as.numeric(two_cell_young_fdr_df[i, 2:1001]) - mean(as.numeric(two_cell_young_fdr_df[i, 2:1001]))) >= abs(observed_values_two_cell_young[i] - mean(as.numeric(two_cell_young_fdr_df[i, 2:1001])))) + 1) / (1000 + 1)
})
p_values_two_sided_two_cell_young
p_values_two_sided_two_cell_young_bonf <- p.adjust(p_values_two_sided_two_cell_young, method = "bonferroni")
p_values_two_sided_two_cell_young_bonf
#[1] 0.005994006 0.005994006 0.005994006 0.005994006 0.005994006 0.005994006

## Fig ED10g -- Two cell clone types, old samples
## two cell clonetypes for old samples:
compiled_clust_df_2cell_old = compiled_clust_df_2cell[compiled_clust_df_2cell$age_binned == ">GW20",]
head(compiled_clust_df_2cell_old)
randomized_df_2cell_old = compiled_clust_df_2cell_old %>% select(Clone_IDs, Clone_size, number_unique_celltypes, upset, two_cell_clone_type)
total_cell_counts <- compiled_clust_df_v1[compiled_clust_df_v1$Clone_size == 2 & compiled_clust_df_v1$age_binned == ">GW20",1:12] %>%
  select(-Clone_IDs) %>%
  summarise(across(everything(), sum))
cell_pool = total_cell_counts

for (i in 1:1000){
  set.seed(i*3)
  new_colname = paste("i",i,sep="")
  # Step 1: Create the character vector (the pool)
  pool = unlist(lapply(names(total_cell_counts), function(cell_type) rep(cell_type, total_cell_counts[[cell_type]])))
  # Step 2: Simpllify the celltypes
  pool[pool %in% c("RG", "EX_IPCs", "IN_IPCs")] = "P"
  pool[pool %in% c("ENs", "IN_local", "IN_OB", "IN_CGE", "IN_MGE")] = "N"
  pool[pool %in% c("OPCs", "Astrocytes", "Oligodendrocytes")] = "G"
  # Step 4: Shuffle pool to ensure random selection
  pool <- sample(pool)
  # Step 5: Make new 2-cell clones
  new_clones <- paste0(pool[seq(1, length(pool) - 1, by = 2)], 
                       pool[seq(2, length(pool), by = 2)])
  clone_indices <- cumsum(clone_sizes)  # End indices for each clone
  start_indices <- c(1, head(clone_indices, -1) + 1)  # Start indices
  # Step 7: Add to the dataframe
  randomized_df_2cell_old[dim(randomized_df_2cell_old)[2]+1] <- new_clones
  colnames(randomized_df_2cell_old)[dim(randomized_df_2cell_old)[2]] = new_colname
}
head(randomized_df_2cell_old)
randomized_df_2cell_old[1:5,1:20]
write.csv(randomized_df_2cell_old, "randomized_1k_two_cell_clones_P_N_G_old.csv")

table(randomized_df_2cell_old$i104)
#  GG  GN  GP  NG  NN  NP  PG  PN  PP 
#. 87  86 133 106 118 182 148 153 252 

randomized_df_2cell_old = data.frame(lapply(randomized_df_2cell_old, function(x) {gsub("GN", "NG", x)}))
randomized_df_2cell_old = data.frame(lapply(randomized_df_2cell_old, function(x) {gsub("GP", "PG", x)}))
randomized_df_2cell_old = data.frame(lapply(randomized_df_2cell_old, function(x) {gsub("NP", "PN", x)}))
table(randomized_df_2cell_old$i104)
#  GG  NG  NN  PG  PN  PP 
#. 87 192 118 281 335 252 

# Reshape the dataframe to long format
df_long_old <- randomized_df_2cell_old %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "clone_type")
head(df_long_old)

# Count occurrences of each clone type in each iteration
df_abundance_old <- df_long_old %>%
  group_by(iteration, clone_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(iteration) %>%  # Ensure normalization happens per iteration
  mutate(relative_abundance = count / sum(count)) %>%
  ungroup()  # Remove grouping for further processing

head(df_abundance_old)

# Compute actual observed frequencies of `two_cell_clone_type`
table(randomized_df_2cell_old$two_cell_clone_type)
#GG  NG  NN  PG  PN  PP 
#187  61 267 212 168 370 
table(randomized_df_2cell_old$two_cell_clone_type)/sum(table(randomized_df_2cell_old$two_cell_clone_type))
df_actual_old = data.frame(table(randomized_df_2cell_old$two_cell_clone_type)/sum(table(randomized_df_2cell_old$two_cell_clone_type)))
colnames(df_actual_old) = c("two_cell_clone_type", "relative_abundance")
df_actual_old
df_abundance_old$clone_type = factor(df_abundance_old$clone_type, levels=c("PP", "NN", "GG", "PN", "PG", "NG"))

# Create boxplot with actual values overlaid
ggplot(df_abundance_old, aes(x = clone_type, y = relative_abundance)) +
  geom_boxplot(outlier.shape = NA, fill = "gray", alpha = 0.5) +  # Boxplot
  geom_point(data = df_actual_old, 
             aes(x = two_cell_clone_type, y = relative_abundance, color = ">GW20"), 
             size = 3) +
  scale_color_manual(name = "Sample age", 
                     values = c(">GW20" = old_col)) +
  theme_minimal() +
  labs(x = "Clone type", y = "Relative abundance", title = "Distribution of two cell clone types >GW20") +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))

## get 95% CIs
two_cell_old_fdr_df = data.frame(Var1=c("PP", "NN", "GG", "PN", "PG", "NG"))
two_cell_old_fdr_df$Var1 = factor(two_cell_old_fdr_df$Var1)
for (i in c(1:1000)){
  colname=paste("i",i,sep="")
  two_cell_old_fdr_df = left_join(two_cell_old_fdr_df, data.frame(table(randomized_df_2cell_old[,i+5])))
  colnames(two_cell_old_fdr_df)[i+1] = colname
}
two_cell_old_fdr_df[1:5,1:10]
two_cell_old_fdr_df = left_join(two_cell_old_fdr_df, data.frame(table(randomized_df_2cell_old$two_cell_clone_type)))
colnames(two_cell_old_fdr_df)[colnames(two_cell_old_fdr_df) == "Freq"] = "observed_two_cell_clone_type"
colnames(two_cell_old_fdr_df)[colnames(two_cell_old_fdr_df) == "Var1"] = "two_cell_clone_type"
two_cell_old_fdr_df[, 1:5]
two_cell_old_fdr_df[is.na(two_cell_old_fdr_df)] = 0

bootstrap_two_cell_old <- two_cell_old_fdr_df %>%
  select(two_cell_clone_type, starts_with("i")) %>%
  pivot_longer(cols = starts_with("i"), names_to = "iteration", values_to = "count") %>%
  group_by(two_cell_clone_type) %>%
  summarise(
    lower_95CI = quantile(count, 0.025),  # Lower 2.5 percentile
    upper_95CI = quantile(count, 0.975),  # Upper 97.5 percentile
    .groups = 'drop'
  )
bootstrap_two_cell_old

for (r in c(1:length(c("PP", "NN", "GG", "PN", "PG", "NG")))){two_cell_old_fdr_df$mean[r] = mean(as.numeric(two_cell_old_fdr_df[r, 2:1001]))}
two_cell_old_fdr_df$mean
two_cell_old_fdr_df[,c(1,1000:1003)]
two_cell_old_fdr_df <- two_cell_old_fdr_df %>%
  left_join(bootstrap_two_cell_old, by = "two_cell_clone_type")
two_cell_old_fdr_df[,c(1,1000:1005)]

two_cell_old_fdr_df$two_cell_clone_type = factor(two_cell_old_fdr_df$two_cell_clone_type, levels=c("PP", "NN", "GG", "PN", "PG", "NG"))
ggplot(two_cell_old_fdr_df, aes(x = two_cell_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = observed_two_cell_clone_type), color = old_col, size = 3) +
  labs(
    title = "Two cell clone types, >GW20 samples",
    x = "Number of unique cell types",
    y = "Number of clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
image = ggplot(two_cell_old_fdr_df, aes(x = two_cell_clone_type)) +
  # Plot the mean of the random distribution
  geom_point(aes(y = mean), color = "darkgray", size = 3, shape = 18) +
  # Plot the confidence intervals as a shaded area (or line)
  geom_errorbar(
    aes(ymin = lower_95CI, ymax = upper_95CI),
    width = 0.2,
    color = "gray",
    size = 1
  ) +
  # Plot the observed value
  geom_point(aes(y = observed_two_cell_clone_type), color = old_col, size = 3) +
  labs(
    title = "Two cell clone types, >GW20 samples",
    x = "Number of unique cell types",
    y = "Number of clones"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=15, hjust=1), axis.title = element_text(size=25), 
        legend.title = element_text(size=25), legend.text = element_text(size=20), 
        plot.title = element_text(size=30, hjust=0.5),
        plot.subtitle = element_text(size=15, hjust=0.5),
        axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1), axis.ticks=element_line(size=1))
ggsave("randomized_vs_actual_two_cell_clone_type_dist_old.svg", image, width=10, height=8)

observed_values_two_cell_old = two_cell_old_fdr_df$observed_two_cell_clone_type
p_values_two_sided_two_cell_old <- sapply(1:length(observed_values_two_cell_old), function(i) {
  (sum(abs(as.numeric(two_cell_old_fdr_df[i, 2:1001]) - mean(as.numeric(two_cell_old_fdr_df[i, 2:1001]))) >= abs(observed_values_two_cell_old[i] - mean(as.numeric(two_cell_old_fdr_df[i, 2:1001])))) + 1) / (1000 + 1)
})
p_values_two_sided_two_cell_old
p_values_two_sided_two_cell_old_bonf <- p.adjust(p_values_two_sided_two_cell_old, method = "bonferroni")
p_values_two_sided_two_cell_old_bonf
#[1] 0.005994006 0.005994006 0.005994006 0.005994006 0.005994006 0.005994006


## Fig ED14 -- EX marker gene expression in VZs (REVISION)
## ED14b-e -- Eze UMAP expression and highlighting of subclusters
library(data.table)
mat = fread('/media/chang/HDD-11/mgkeefe/240426_eze_data_ucsc/exprMatrix.tsv.gz')
meta = read.table('/media/chang/HDD-11/mgkeefe/240426_eze_data_ucsc/meta.tsv', header=T, sep='\t', as.is=T, row.names=1)
head(meta)
genes = mat[,1][[1]]
head(genes)
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
eze_data_v2 = CreateSeuratObject(counts = mat, project = "CS_data", meta.data=meta)
eze_data_v2
#An object of class Seurat 
#19368 features across 45156 samples within 1 assay 
#Active assay: RNA (19368 features, 0 variable features)
#1 layer present: counts
rm(mat)
rm(meta)

eze_data_v2@assays$RNA$counts
RowsNA<-names(which(rowSums(is.na(eze_data_v2@assays$RNA$counts))>0))
RowsNA
'%!in%' <- function(x,y)!('%in%'(x,y)) #this is a NOT IN function
RowsKEEP<-rownames(eze_data_v2)[rownames(eze_data_v2) %!in% RowsNA]
eze_data_v2<-subset(eze_data_v2,features=RowsKEEP)
# normalize data
eze_data_v2 <- NormalizeData(eze_data_v2, normalization.method = "LogNormalize", scale.factor = 10000)
# scale data based on normalization
all.genes <- rownames(eze_data_v2)
eze_data_v2 <- ScaleData(eze_data_v2, features = all.genes)
#vz_pipseq <- ScaleData(vz_pipseq, features = all.genes, vars.to.regress = "nCount_RNA")
# find variably expressed genes
eze_data_v2 <- FindVariableFeatures(eze_data_v2, selection.method = "vst", nfeatures = 2000)
# perform dimensionality reduction
eze_data_v2 <- RunPCA(eze_data_v2, features = VariableFeatures(object = eze_data_v2))
ElbowPlot(eze_data_v2)
# cluster cells
eze_data_v2 <- FindNeighbors(eze_data_v2, dims = 1:15)
eze_data_v2 <- FindClusters(eze_data_v2, resolution = .5)
# make graph
eze_data_v2 <- RunUMAP(eze_data_v2, dims = 1:15)
DimPlot(eze_data_v2, reduction = "umap",label=T)
DimPlot(eze_data_v2, group.by="orig.ident")

## Fig ED14b -- Eze seurat clusters
DimPlot(eze_data, label=T, pt.size=0.5) +
  ggtitle("First Timester -- 10x") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(eze_data, label=T, pt.size=0.5) +
  ggtitle("First Timester -- 10x") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="eze_data_dimplot.svg", plot=image, width=10, height=8)
ggsave(file="eze_data_dimplot.png", plot=image, width=10, height=8)

## Fig ED14c -- Eze sample ID
DimPlot(eze_data, label=F, pt.size=0.5, group.by="orig.ident") +
  ggtitle("First Timester -- 10x") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(eze_data, label=F, pt.size=0.5, group.by="orig.ident") +
  ggtitle("First Timester -- 10x") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="eze_data_dimplot_by_sample.svg", plot=image, width=10, height=8)
ggsave(file="eze_data_dimplot_by_sample.png", plot=image, width=10, height=8)

## Fig ED14d -- Eze marker gene UMAPS
FeaturePlot(eze_data, features=c("HES1"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(eze_data, features=c("HES1"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="eze_data_hes1.svg", plot=image, width=4, height=3.2)
ggsave(file="eze_data_hes1.png", plot=image, width=4, height=3.2)
FeaturePlot(eze_data, features=c("MKI67"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(eze_data, features=c("MKI67"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="eze_data_mki67.svg", plot=image, width=4, height=3.2)
ggsave(file="eze_data_mki67.png", plot=image, width=4, height=3.2)
FeaturePlot(eze_data, features=c("EOMES"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(eze_data, features=c("EOMES"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="eze_data_eomes.svg", plot=image, width=4, height=3.2)
ggsave(file="eze_data_eomes.png", plot=image, width=4, height=3.2)
FeaturePlot(eze_data, features=c("SLC17A7"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(eze_data, features=c("SLC17A7"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="eze_data_slc17a7.svg", plot=image, width=4, height=3.2)
ggsave(file="eze_data_slc17a7.png", plot=image, width=4, height=3.2)

## Fig ED14e -- Eze neuron highlights
DimPlot(eze_data, pt.size=0.5, label=T, cells.highlight = rownames(eze_data@meta.data)[eze_data@meta.data$seurat_clusters %in% c(4,8,10,12)]) +
  ggtitle("First Trimester -- 10x") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(eze_data, pt.size=0.5, label=T, cells.highlight = rownames(eze_data@meta.data)[eze_data@meta.data$seurat_clusters %in% c(4,8,10,12)]) +
  ggtitle("First Trimester -- 10x") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="eze_data_neuron_highlight.svg", plot=image, width=10, height=8)
ggsave(file="eze_data_neuron_highlight.png", plot=image, width=10, height=8)

## ED14f-i -- 10x midgestation VZ UMAP expression and highlighting of subclusters
subdissect_gw18_vz.data = Read10X('/media/chang/HDD-4/mgkeefe/240521_subdissect_data/GW18VZ/')
subdissect_gw18_vz = CreateSeuratObject(counts = subdissect_gw18_vz.data, project = "subdissect_gw18_vz", min.cells = 10, min.features = 200)
subdissect_gw18_vz
#An object of class Seurat 
#16215 features across 10040 samples within 1 assay 
subdissect_gw18_vz[["percent.mt"]] = PercentageFeatureSet(subdissect_gw18_vz, pattern="^MT-")
VlnPlot(subdissect_gw18_vz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
subdissect_gw18_vz <- subset(subdissect_gw18_vz, subset = 
                           nFeature_RNA > 50 & nFeature_RNA < 30000 &
                           nFeature_RNA > 500 & nCount_RNA < 30000 & percent.mt < 5)
VlnPlot(subdissect_gw18_vz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# normalize data
subdissect_gw18_vz <- NormalizeData(subdissect_gw18_vz, normalization.method = "LogNormalize", scale.factor = 10000)
# scale data based on normalization
all.genes <- rownames(subdissect_gw18_vz)
subdissect_gw18_vz <- ScaleData(subdissect_gw18_vz, features = all.genes)
# find variably expressed genes
subdissect_gw18_vz <- FindVariableFeatures(subdissect_gw18_vz, selection.method = "vst", nfeatures = 2000)
# perform dimensionality reduction
subdissect_gw18_vz <- RunPCA(subdissect_gw18_vz, features = VariableFeatures(object = subdissect_gw18_vz))
ElbowPlot(subdissect_gw18_vz)
# cluster cells
subdissect_gw18_vz <- FindNeighbors(subdissect_gw18_vz, dims = 1:10)
subdissect_gw18_vz <- FindClusters(subdissect_gw18_vz, resolution = .5)
# make graph
subdissect_gw18_vz <- RunUMAP(subdissect_gw18_vz, dims = 1:10)
DimPlot(subdissect_gw18_vz, reduction = "umap",label=T)
FeaturePlot(subdissect_gw18_vz, features=c("CRYAB", "HOPX", "EOMES", "NEUROD2"))
FeaturePlot(subdissect_gw18_vz, features=c("GAD2", "SCGN", "PAX6", "LHX6"))
FeaturePlot(subdissect_gw18_vz, features=c("TBR1", "NFIB", "SORCS1", "LMO3"))
FeaturePlot(subdissect_gw18_vz, features=c("nCount_RNA"), max.cutoff="q95")
FeaturePlot(subdissect_gw18_vz, features=c("percent.mt"), max.cutoff="q95")
FeaturePlot(subdissect_gw18_vz, features=c("RPS4Y1", "TTR", "RPS2", "ATP5D", "GORASP1", "ARCN1"), ncol=3)
#check DEGs between the first and second groups -- is it just cell stress, or what's happening?
DimPlot(subdissect_gw18_vz, reduction = "umap",label=T)
subdissect_gw18_vz = AddMetaData(subdissect_gw18_vz, "group1", col.name = "supercluster")
subdissect_gw18_vz@meta.data$supercluster[subdissect_gw18_vz@meta.data$seurat_clusters %in% c(3,8,9)] = "group2"
DimPlot(subdissect_gw18_vz, group.by="supercluster")
subdissect_gw18_vz = SetIdent(subdissect_gw18_vz, value = "supercluster")
#devtools::install_github('immunogenomics/presto')
vz_supercluster_markers = FindAllMarkers(subdissect_gw18_vz, min.diff.pct = 0.3, only.pos = T)
vz_supercluster_markers
FeaturePlot(subdissect_gw18_vz, features=c("PPP1R17"))
setwd('/media/chang/HDD-4/mgkeefe/R_work_in_progress_cajal/240514_resubmission_figs')
write.csv(vz_supercluster_markers, 'subdissect_gw18_vz_supercluster_markers.csv')


subdissect_gw21_vz.data = Read10X('/media/chang/HDD-4/mgkeefe/240521_subdissect_data/GW21VZ/')
subdissect_gw21_vz = CreateSeuratObject(counts = subdissect_gw21_vz.data, project = "subdissect_gw21_vz", min.cells = 10, min.features = 200)
subdissect_gw21_vz
#An object of class Seurat 
#16215 features across 10040 samples within 1 assay 
subdissect_gw21_vz[["percent.mt"]] = PercentageFeatureSet(subdissect_gw21_vz, pattern="^MT-")
subdissect_gw21_vz[["percent.ribo"]] = PercentageFeatureSet(subdissect_gw21_vz, pattern="^RPL")
VlnPlot(subdissect_gw21_vz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
subdissect_gw21_vz <- subset(subdissect_gw21_vz, subset = 
                           nFeature_RNA > 50 & nFeature_RNA < 30000 &
                           nFeature_RNA > 500 & nCount_RNA < 30000 & percent.mt < 5)
VlnPlot(subdissect_gw21_vz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# normalize data
subdissect_gw21_vz <- NormalizeData(subdissect_gw21_vz, normalization.method = "LogNormalize", scale.factor = 10000)
# scale data based on normalization
all.genes <- rownames(subdissect_gw21_vz)
subdissect_gw21_vz <- ScaleData(subdissect_gw21_vz, features = all.genes)
# find variably expressed genes
subdissect_gw21_vz <- FindVariableFeatures(subdissect_gw21_vz, selection.method = "vst", nfeatures = 2000)
# perform dimensionality reduction
subdissect_gw21_vz <- RunPCA(subdissect_gw21_vz, features = VariableFeatures(object = subdissect_gw21_vz))
ElbowPlot(subdissect_gw21_vz)
# cluster cells
subdissect_gw21_vz <- FindNeighbors(subdissect_gw21_vz, dims = 1:10)
subdissect_gw21_vz <- FindClusters(subdissect_gw21_vz, resolution = .5)
# make graph
subdissect_gw21_vz <- RunUMAP(subdissect_gw21_vz, dims = 1:10)
DimPlot(subdissect_gw21_vz, reduction = "umap",label=T)
FeaturePlot(subdissect_gw21_vz, features=c("CRYAB", "HOPX", "EOMES", "NEUROD2"))
FeaturePlot(subdissect_gw21_vz, features=c("TBR1"))
FeaturePlot(subdissect_gw21_vz, features=c("GAD2", "SCGN", "PAX6", "LHX6"))
FeaturePlot(subdissect_gw21_vz, features=c("NR4A2", "NFIB", "SORCS1", "LMO3"))
FeaturePlot(subdissect_gw21_vz, features=c("nCount_RNA"), max.cutoff="q98")
FeaturePlot(subdissect_gw21_vz, features=c("percent.mt"), max.cutoff="q98")
FeaturePlot(subdissect_gw21_vz, features=c("percent.ribo"), max.cutoff="q98")

subdissect_gw18_vz = SCTransform(subdissect_gw18_vz, vars.to.regress = "percent.mt", verbose=FALSE, vst.flavor = "v2")
subdissect_gw18_vz = RunPCA(subdissect_gw18_vz, npcs=30, verbose=FALSE)
subdissect_gw21_vz = SCTransform(subdissect_gw21_vz, vars.to.regress = "percent.mt", verbose=FALSE, vst.flavor = "v2")
subdissect_gw21_vz = RunPCA(subdissect_gw21_vz, npcs=30, verbose=FALSE)

subdissect_merge.list = list(
  subdissect_gw18_vz=subdissect_gw18_vz,
  subdissect_gw21_vz=subdissect_gw21_vz)
features <- SelectIntegrationFeatures(object.list = subdissect_merge.list, nfeatures = 3000)
subdissect_merge.list <- PrepSCTIntegration(object.list = subdissect_merge.list, anchor.features = features)
ls.anchors <- FindIntegrationAnchors(object.list = subdissect_merge.list, normalization.method = "SCT",
                                     anchor.features = features, reduction="rpca") #using rpca to speed up
subdissect_merge <- IntegrateData(anchorset = ls.anchors, normalization.method = "SCT")
subdissect_merge <- RunPCA(subdissect_merge, verbose = FALSE)
subdissect_merge <- RunUMAP(subdissect_merge, reduction = "pca", dims = 1:30, verbose = FALSE)
subdissect_merge <- FindNeighbors(subdissect_merge, reduction = "pca", dims = 1:30)
subdissect_merge <- FindClusters(subdissect_merge, resolution = 0.4)

DefaultAssay(subdissect_merge) = "RNA"
subdissect_merge <- NormalizeData(subdissect_merge, normalization.method = "LogNormalize", scale.factor = 10000)
DimPlot(subdissect_merge, label=T)
DimPlot(subdissect_merge, group.by="orig.ident")
FeaturePlot(subdissect_merge, features=c("MKI67", "VIM", "EOMES", "NEUROD2", "TBR1", "FOXP2", "GAD2", "LHX6", "NR2F2"))
## delete doublet/low quality cells in right clusters
subdissect_merge = subset(subdissect_merge, seurat_clusters %in% c(0,2,3,4,5,6,9,10,12,14,16,17,20))

## Fig ED14f -- 10x VZ UMAP
DimPlot(subdissect_merge_vz, label=T, pt.size=0.5) +
  ggtitle("Midgestation VZ (10x)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(subdissect_merge_vz, label=F, pt.size=0.5) +
  ggtitle("Midgestation VZ (10x)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="subdissect_merge_vz_dimplot.svg", plot=image, width=10, height=8)
ggsave(file="subdissect_merge_vz_dimplot.png", plot=image, width=10, height=8)

## Fig ED14g -- 10X VZ sample UMAP
DimPlot(subdissect_merge_vz, label=F, pt.size=0.5, group.by="sample_age") +
  ggtitle("Sample age (Gestational week)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(subdissect_merge_vz, label=F, pt.size=0.5, group.by="sample_age") +
  ggtitle("Sample age (Gestational week)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="subdissect_merge_vz_dimplot_by_sample.svg", plot=image, width=10, height=8)
ggsave(file="subdissect_merge_vz_dimplot_by_sample.png", plot=image, width=10, height=8)

## Fig ED14h -- 10X VZ marker genes
FeaturePlot(subdissect_merge_vz, features=c("HES1"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(subdissect_merge_vz, features=c("HES1"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="subdissect_merge_vz_hes1.svg", plot=image, width=4, height=3.2)
ggsave(file="subdissect_merge_vz_hes1.png", plot=image, width=4, height=3.2)
FeaturePlot(subdissect_merge_vz, features=c("MKI67"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(subdissect_merge_vz, features=c("MKI67"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="subdissect_merge_vz_mki67.svg", plot=image, width=4, height=3.2)
ggsave(file="subdissect_merge_vz_mki67.png", plot=image, width=4, height=3.2)
FeaturePlot(subdissect_merge_vz, features=c("EOMES"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(subdissect_merge_vz, features=c("EOMES"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="subdissect_merge_vz_eomes.svg", plot=image, width=4, height=3.2)
ggsave(file="subdissect_merge_vz_eomes.png", plot=image, width=4, height=3.2)
FeaturePlot(subdissect_merge_vz, features=c("SLC17A7"), order=T, min.cutoff="q10", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(subdissect_merge_vz, features=c("SLC17A7"), order=T, min.cutoff="q10", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="subdissect_merge_vz_slc17a7.svg", plot=image, width=4, height=3.2)
ggsave(file="subdissect_merge_vz_slc17a7.png", plot=image, width=4, height=3.2)

## Fig ED14i -- 10X VZ neuron highlights
DimPlot(subdissect_merge_vz, pt.size=0.5, label=T, cells.highlight = rownames(subdissect_merge_vz@meta.data)[subdissect_merge_vz@meta.data$seurat_clusters %in% c(0,6)]) +
  ggtitle("Midgestation VZ (10x)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(subdissect_merge_vz, pt.size=0.5, label=T, cells.highlight = rownames(subdissect_merge_vz@meta.data)[subdissect_merge_vz@meta.data$seurat_clusters %in% c(0,6)]) +
  ggtitle("Midgestation VZ (10x)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="subdissect_merge_vz_neuron_highlight.svg", plot=image, width=10, height=8)
ggsave(file="subdissect_merge_vz_neuron_highlight.png", plot=image, width=10, height=8)

## ED14j-m -- Pipseq midgestation VZ UMAP expression and highlighting of subclusters
vz_pipseq.data = Read10X('/media/chang/HDD-4/mgkeefe/240507_GW19VZ_pipseq/filtered_matrix/sensitivity_4')
vz_pipseq = CreateSeuratObject(counts = vz_pipseq.data, project = "VZ_pipseq", min.cells = 10, min.features = 200)
vz_pipseq
#An object of class Seurat 
#24339 features across 12982 samples within 1 assay 
vz_pipseq[["percent.mt"]] = PercentageFeatureSet(vz_pipseq, pattern="^MT-")
VlnPlot(vz_pipseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vz_pipseq <- subset(vz_pipseq, subset = 
                      nFeature_RNA > 50 & nFeature_RNA < 30000 &
                      nFeature_RNA > 800 & nCount_RNA < 20000 & percent.mt < 10)
VlnPlot(vz_pipseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# normalize data
vz_pipseq <- NormalizeData(vz_pipseq, normalization.method = "LogNormalize", scale.factor = 10000)
# scale data based on normalization
all.genes <- rownames(vz_pipseq)
vz_pipseq <- ScaleData(vz_pipseq, features = all.genes)
# find variably expressed genes
vz_pipseq <- FindVariableFeatures(vz_pipseq, selection.method = "vst", nfeatures = 2000)
# perform dimensionality reduction
vz_pipseq <- RunPCA(vz_pipseq, features = VariableFeatures(object = vz_pipseq))
ElbowPlot(vz_pipseq)
# cluster cells
vz_pipseq <- FindNeighbors(vz_pipseq, dims = 1:10)
vz_pipseq <- FindClusters(vz_pipseq, resolution = .5)
# make graph
vz_pipseq <- RunUMAP(vz_pipseq, dims = 1:10)
DimPlot(vz_pipseq, reduction = "umap",label=T)
FeaturePlot(vz_pipseq, features=c("CRYAB", "HOPX", "EOMES", "NEUROD2"))
FeaturePlot(vz_pipseq, features=c("GAD2", "SCGN", "PAX6", "LHX6"))
FeaturePlot(vz_pipseq, features=c("nCount_RNA"), max.cutoff="q98")

## Fig ED14j -- Pipseq VZ seurat clusters
DimPlot(vz_pipseq, label=T, pt.size=0.5) +
  ggtitle("GW19 VZ -- Pipseq") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(vz_pipseq, label=T, pt.size=0.5) +
  ggtitle("GW19 VZ -- Pipseq") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="pipseq_vz_only_dimplot.svg", plot=image, width=10, height=8)
ggsave(file="pipseq_vz_only_dimplot.png", plot=image, width=10, height=8)

## Fig ED14l -- Pipseq VZ sample
DimPlot(vz_pipseq, label=F, pt.size=0.5, group.by="orig.ident") +
  ggtitle("Sample age (Gestational week)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(vz_pipseq, label=F, pt.size=0.5, group.by="orig.ident") +
  ggtitle("Sample age (Gestational week)") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="pipseq_vz_dimplot_by_sample.svg", plot=image, width=10, height=8)
ggsave(file="pipseq_vz_dimplot_by_sample.png", plot=image, width=10, height=8)


## Fig ED14k -- Pipseq VZ marker genes
FeaturePlot(vz_pipseq, features=c("HES1"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(vz_pipseq, features=c("HES1"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="pipseq_vz_only_hes1.svg", plot=image, width=4, height=3.2)
ggsave(file="pipseq_vz_only_hes1.png", plot=image, width=4, height=3.2)
FeaturePlot(vz_pipseq, features=c("MKI67"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(vz_pipseq, features=c("MKI67"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="pipseq_vz_only_mki67.svg", plot=image, width=4, height=3.2)
ggsave(file="pipseq_vz_only_mki67.png", plot=image, width=4, height=3.2)
FeaturePlot(vz_pipseq, features=c("EOMES"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(vz_pipseq, features=c("EOMES"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="pipseq_vz_only_eomes.svg", plot=image, width=4, height=3.2)
ggsave(file="pipseq_vz_only_eomes.png", plot=image, width=4, height=3.2)
FeaturePlot(vz_pipseq, features=c("SLC17A7"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = FeaturePlot(vz_pipseq, features=c("SLC17A7"), order=T, min.cutoff="q25", pt.size=0.5) +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="pipseq_vz_only_slc17a7.svg", plot=image, width=4, height=3.2)
ggsave(file="pipseq_vz_only_slc17a7.png", plot=image, width=4, height=3.2)

## Fig ED14m -- VZ pipseq neurons highlight
DimPlot(vz_pipseq, pt.size=0.5, label=T, cells.highlight = rownames(vz_pipseq@meta.data)[vz_pipseq@meta.data$seurat_clusters %in% c(0)]) +
  ggtitle("GW19 VZ -- Pipseq") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
image = DimPlot(vz_pipseq, pt.size=0.5, label=T, cells.highlight = rownames(vz_pipseq@meta.data)[vz_pipseq@meta.data$seurat_clusters %in% c(0)]) +
  ggtitle("GW19 VZ -- Pipseq") +
  theme(axis.line =element_blank(), legend.position="right", axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust=0.5, size=30), legend.text = element_text(size=20))
ggsave(file="pipseq_vz_neuron_highlight.svg", plot=image, width=10, height=8)
ggsave(file="pipseq_vz_neuron_highlight.png", plot=image, width=10, height=8)

## ED14n-p -- heatmap of EX gene expression in each sample
## Fig ED14n -- Eze heatmap
minimal_sp_panel_v2_ctgf = c('EOMES', 'PPP1R17', 'NEUROD1', 'NEUROD2', 'NEUROD6', 'NFIA', 'NFIB', 'LMO3', 'TLE4', 'FOXP2', 'NR4A2', 'DPP10', 'EPB41L4A', 'ZFHX3', 'TOX', 'NEGR1', 'CDH10', 'TMEM178A', 'CRYM', 'SORCS1',  'TRPM3', 'TBR1', 'HS3ST4',  'HS3ST5', 'OCA2', 'SEMA3E', 'CTGF', 'ST18')
DoHeatmap(subset(eze_data, seurat_clusters %in% c(4,8,10,12)), features=minimal_sp_panel_v2_ctgf) +
  theme(legend.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        axis.title = element_text(size=20), 
        axis.text.y = element_text(size=15))
ggsave('eze_data_clust4_8_10_12_panel_order_v2.png', device = "png", width=8, height=6)
ggsave('eze_data_clust4_8_10_12_panel_order_v2.svg', device = "svg", width=8, height=6)

## Fig ED14o -- 10x VZ heatmap
minimal_sp_panel_v2_ctgf = c('EOMES', 'PPP1R17', 'NEUROD1', 'NEUROD2', 'NEUROD6', 'NFIA', 'NFIB', 'LMO3', 'TLE4', 'FOXP2', 'NR4A2', 'DPP10', 'EPB41L4A', 'ZFHX3', 'TOX', 'NEGR1', 'CDH10', 'TMEM178A', 'CRYM', 'SORCS1',  'TRPM3', 'TBR1', 'HS3ST4',  'HS3ST5', 'OCA2', 'SEMA3E', 'CTGF', 'ST18')
DoHeatmap(subset(subdissect_merge_vz, seurat_clusters %in% c(0,6)), features=minimal_sp_panel_v2_ctgf) +
  theme(legend.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        axis.title = element_text(size=20), 
        axis.text.y = element_text(size=15))
ggsave('subdissect_merge_vz_clust4_8_10_12_panel_order_v2.png', device = "png", width=8, height=6)
ggsave('subdissect_merge_vz_clust4_8_10_12_panel_order_v2.svg', device = "svg", width=8, height=6)

## Fig ED14p -- Pipseq VZ heatmap
minimal_sp_panel_v2 = c('EOMES', 'PPP1R17', 'NEUROD1', 'NEUROD2', 'NEUROD6', 'NFIA', 'NFIB', 'LMO3', 'TLE4', 'FOXP2', 'NR4A2', 'DPP10', 'EPB41L4A', 'ZFHX3', 'TOX', 'NEGR1', 'CDH10', 'TMEM178A', 'CRYM', 'SORCS1',  'TRPM3', 'TBR1', 'HS3ST4',  'HS3ST5', 'OCA2', 'SEMA3E', 'CCN2', 'ST18')
DoHeatmap(subset(vz_pipseq, seurat_clusters %in% c(0)), features=minimal_sp_panel_v2) +
  theme(legend.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        axis.title = element_text(size=20), 
        axis.text.y = element_text(size=15))
ggsave('vz_pipseq_clust0_panel_order_v2.png', device = "png", width=8, height=6)
DoHeatmap(subset(vz_pipseq, seurat_clusters %in% c(0)), features=minimal_sp_panel_v2) +
  theme(legend.text = element_text(size=15), 
        legend.title = element_text(size=20), 
        axis.title = element_text(size=20), 
        axis.text.y = element_text(size=15))
ggsave('vz_pipseq_clust0_panel_order_v2.svg', device = "svg", width=8, height=6)
require(ggplot2)
require(data.table)
library(dplyr)
require(project.init)
project.init2("artemether")


#########
color_all<- c("#4D4D4E","#1CBCC2", "#4FB763", "#904198")
color_FOX<- c("#4D4D4E","#1CBCC2")


fig4_human<- as.data.frame(fread(dirout("FIG_03_DiffEXPR/", "MetaData_Export_human.tsv")))

fig4_mouse<- as.data.frame(fread(dirout("FIG_03_DiffEXPR/", "MetaData_Export_mouse.tsv")))

#FIGURE 4a- Violin plot for beta cells, human INS and mouse Ins1 and Ins2 with FoxOi

fig4a_human_beta <- fig4_human %>%
  filter(celltype2 == "Beta" & treatment%in%c("FoxOi", "DMSO"))

fig4a_mouse_beta <- fig4_mouse %>%
  filter(celltype2 == "Beta" & treatment%in%c("FoxOi", "DMSO"))

ggplot(fig4a_human_beta,aes(x=replicate, y=INS, color=treatment, fill=treatment)) + geom_violin()+   
  scale_color_manual(values = color_FOX) + scale_fill_manual(values = color_FOX)+ theme_bw() +ylim(8,14)+
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")+
  labs(y = "Log10 INS") +   labs(x = "Replicate")

ggplot(fig4a_mouse_beta,aes(x=replicate, y=Ins1, color=treatment, fill=treatment)) + geom_violin()+   
  scale_color_manual(values = color_FOX) + scale_fill_manual(values = color_FOX)+ theme_bw() +ylim(8,14)+
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")+
  labs(y = "Log10 INS") +   labs(x = "Replicate")

ggplot(fig4a_mouse_beta,aes(x=replicate, y=Ins2, color=treatment, fill=treatment)) + geom_violin()+   
  scale_color_manual(values = color_FOX) + scale_fill_manual(values = color_FOX)+ theme_bw() +ylim(8,14)+
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")+
  labs(y = "Log10 INS") +   labs(x = "Replicate")


#FIGURE 4C- bubble heatmap for key genes in alpha and beta cell identity

library(gridExtra)

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
}


fig4c<- as.data.frame(fread(dirout("FIG_03_DiffEXPR/", "DiffGenes_byReplicate_AlphaBetaOnly.tsv")))

genesbetafoxo <- c("NEUROD1","INS/Ins2", "Ins1", "ISL1", "NKX6-1", "NKX2-2", "FOXA2", 
                   "MAFA", "FOXO3", "FOXO1", "LY96", "CRYBA2", "C2CD4A", "PPP1CC", "PDX1",
                   "STX4", "YWHAB", "RHOQ", "PAX6", "INSR", "IGF2")

genesalphafoxo <- c("GCG","ARX", "INS/Ins2", "Ins1", "GC", "NEUROD1", "TTR", "PAX6",
                    "PCSK2", "MAFB", "TM4SF4", "DDIT3", "ID1", "NPTX2")

fig4c_alpha <- fig4c %>%
  filter(celltype == "Alpha" & treatment =="FoxOi"& human%in%genesalphafoxo)

fig4c_beta <- fig4c %>%
  filter(celltype == "Alpha" & treatment =="FoxOi"& human%in%genesbetafoxo)

fig4c_alpha<-ggplot(fig4c_alpha, aes(x=replicate, y=human, color=logFC, size=pmin(10, -log10(qval)))) + 
  geom_point() + facet_grid(~org , scales="free", space="free") +
  scale_color_gradient2(name="logFC", high="indianred", low="navyblue") + 
  theme_bw(16) + 
  scale_size_continuous(name=expression(-log[10](qval))) +
  ylab("") +  facet_grid(~org)

fig4c_beta<-ggplot(fig4c_beta, aes(x=replicate, y=human, color=logFC, size=pmin(10, -log10(qval)))) + 
  geom_point() + facet_grid(~org , scales="free", space="free") +
  scale_color_gradient2(name="logFC", high="indianred", low="navyblue") + 
  theme_bw(16) + 
  scale_size_continuous(name=expression(-log[10](qval))) +
  ylab("") +  facet_grid(~org)

grid_arrange_shared_legend (fig4c_beta, fig4c_alpha, nrow=2)

#FIGURE 5A - Fraction of alpha cells that express insulin. Use same files as fig4, but filter for alpha cells: 

color_5a<- c("#904198", "#4D4D4E","#1CBCC2")
#Filter treatments and samples
fig5a_human <- fig4_human %>%
  filter(celltype2 == "Alpha" & treatment%in%c("FoxOi", "DMSO", "Artemether"))

fig5a_mouse <- fig4_mouse %>%
  filter(celltype2 == "Alpha" & treatment%in%c("FoxOi", "DMSO", "Artemether"))
 #plot distribution

fig5a_human<-ggplot(fig5a_human, aes(x=INS, color=treatment))+ stat_ecdf()

fig5a_mouse_Ins1<-ggplot(fig5a_mouse, aes(x=Ins1, color=treatment))+
  stat_ecdf()

fig5a_mouse_Ins2<-ggplot(fig5a_mouse, aes(x=Ins2, color=treatment))+
  stat_ecdf()

#inverse plot
fig5a_human + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log INS") + xlim(0,15) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
scale_color_manual(values = color_5a) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")

fig5a_mouse_Ins1 + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log Ins1") + xlim(0,15) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+  
  scale_color_manual(values = color_5a) + theme_bw() +
  theme(text = element_text(size=19), axis.text.x = element_text(size= 16),axis.text.y = element_text(size= 16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")

fig5a_mouse_Ins2 + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log Ins2") + xlim(0,15) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+  
  scale_color_manual(values = color_5a) + theme_bw() +
  theme(text = element_text(size=19), axis.text.x = element_text(size= 16),axis.text.y = element_text(size= 16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")


#FIGure 5C -Bubble heatmap for key genes in alpha Ins+ cells vs alpha Ins- cells.

human_markers<-as.data.frame(fread(dirout("36_01_AlphaInsHighVsLow/", "humanMarkers.tsv")))
mouse_markers<-as.data.frame(fread(dirout("36_01_AlphaInsHighVsLow/", "mouseMarkers.tsv")))

#homology file to convert gene mouse names into human
hom<-as.data.frame(fread(dirout("FIG_03_DiffEXPR/", "Homology_unique.tsv")))

#merge homology file with mouse markers, to get human names
mouse_homology<-merge(hom, mouse_markers, by="gene_mouse")

mouse_names <- names(mouse_homology) %in% c("gene_mouse")
mouse_human_names<- mouse_homology[!mouse_names]

#add human and mouse markers, to have them in the same plot
mouse_human_markers<-rbind(human_markers, mouse_human_names)


#alpha and beta marker genes
genesalphains<- c("GCG", "IAPP", "DLK1","PAX4","NGN3", "ABCC8","PDX1", "MAFA", "NKX6-1", "NKX2-2", "ARX", "TTR", "IGF1R", "SLC2A2", "FOXO1")

fig5c <- mouse_human_markers %>%
  filter(gene_human%in%genesalphains)

ggplot(fig5c, aes(x=group, y=gene_human, color=logFC, size=pmin(9, -log10(p_val)))) + 
  geom_point() + facet_grid(~org , scales="free", space="free") +
  scale_color_gradient2(name="logFC", high="indianred", low="navyblue") + 
  theme_bw(16) + 
  scale_size_continuous(name=expression(-log[10](p_val))) +
  ylab("")

#FIGURE 5D, same data as figure 4A and 5A, but for beta
color_5d<- c("#904198", "#4D4D4E","#1CBCC2")

fig5d_hbeta <- fig4_human %>%
  filter(celltype == "Beta" & treatment%in%c("FoxOi", "DMSO", "Artemether"))

fig5d_mbeta <- fig4_mouse %>%
  filter(celltype == "Beta" & treatment%in%c("FoxOi", "DMSO", "Artemether"))

ggplot(fig5d_hbeta, aes(x=INS, color=treatment))+ stat_ecdf(size=1.2)+ 
  labs(x = "Log INS") + xlim(5.5,13.5) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_5d) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")


ggplot(fig5d_mbeta, aes(x=Ins1, color=treatment))+ stat_ecdf(size=1.2)+ 
  labs(x = "Log Ins1") + xlim(5.5,13.5) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_5d) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")

ggplot(fig5d_mbeta, aes(x=Ins2, color=treatment))+ stat_ecdf(size=1.2)+ 
  labs(x = "Log Ins2") + xlim(5.5,13.5) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_5d) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")
#inverse plots outside of R




#FIGURE 5E- Bubble plot for expression of key beta cell genes and markers of de-differentiaiton in mouse and human beta cells, Sam data as figure 4c.

genesbetahuman <- c("GCG","UCN3", "HN1", "IMPDH2", "HES1", "NEUROD1","INS/Ins2", "ISL1", "NKX6-1", "NKX2-2", "FOXO3", "PDX1", "FOXO1", "SCG5", "IGF1R", "ABCC8", "GCGR", "TTR", "PCSK2", "MAFB", "INSR", "VAMP2", "CRYBA2", "CDCD4A", "PPP1CC", "YWHAB", "SLC2A2", "RFX6", "SPP1", "PTEN", "PTPN1")

pDat5e_fox <- fig4c %>%
  filter(human%in%genesbetahuman &celltype == "Beta"& treatment%in%c("FoxOi") & org%in%c("mouse", "human"))

pDat5e_art <- fig4c %>%
  filter(human%in%genesbetahuman &celltype == "Beta"& treatment%in%c("A10") & org%in%c("mouse", "human"))

pDat5e <- hierarch.ordering(pDat5e, "celltype", "group","logFC")

fig5efoxo<-ggplot(pDat5e_art, aes(x=replicate, y=human, color=logFC, size=pmin(10, -log10(qval)))) + 
  geom_point() + facet_grid(~org , scales="free", space="free") +
  scale_color_gradient2(name="logFC", high="indianred", low="navyblue") + 
  theme_bw(16) + 
  scale_size_continuous(name=expression(-log[10](qval))) +
  ylab("") +  xlab("A10") +facet_grid(~org)+theme(legend.position="left", legend.box = "vertical")

fig5eart<-ggplot(pDat5e_fox, aes(x=replicate, y=human, color=logFC, size=pmin(10, -log10(qval)))) + 
  geom_point() + facet_grid(~org , scales="free", space="free") +
  scale_color_gradient2(name="logFC", high="indianred", low="navyblue") + 
  theme_bw(16) + 
  scale_size_continuous(name=expression(-log[10](qval))) +
  ylab("") +  facet_grid(~org)+theme(legend.position="left", legend.box = "vertical")

grid_arrange_shared_legend (fig5eart, fig5efoxo)


#FIGURES S7a- insulin expression in alpha cells, separated by replicates 

#same as figure 5a, but facet_grid by replicate to show different replicates

fig5a_human + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log Ins2") + xlim(0,15) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_5a) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")+
  facet_grid(~replicate)

fig5a_mouse_Ins1 + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log Ins2") + xlim(0,15) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_5a) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")+
  facet_grid(~replicate)

fig5a_mouse_Ins2 + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log Ins2") + xlim(0,15) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_5a) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")+
  facet_grid(~replicate)

#FIGURE S8- Insulin expression in alpha and beta cells in human donor 3
color_supp <- c("#904198","#4D4D4E")

fig_supp8_h3<- as.data.frame(fread(dirout("FIG_03_DiffEXPR/", "MetaData_Export_human3.tsv")))


#filter treatments and samples
figsupp8_h3_alpha_72 <- fig_supp8_h3 %>%
  filter(celltype2 == "Alpha" & time=="72" & treatment%in%c("FoxOi", "DMSO", "Artemether"))

figsupp8_h3_beta_72 <- fig_supp8_h3 %>%
  filter(celltype2 == "Beta" & time=="72" & treatment%in%c("FoxOi", "DMSO", "Artemether"))

figsupp8_h3_alpha_36 <- fig_supp8_h3 %>%
  filter(celltype2 == "Alpha" & time=="36" & treatment%in%c("FoxOi", "DMSO", "Artemether"))

figsupp8_h3_beta_36 <- fig_supp8_h3 %>%
  filter(celltype2 == "Beta" & time=="36" & treatment%in%c("FoxOi", "DMSO", "Artemether"))

#plot distributions and inverse plot
figsupp8_h3_alpha_72<-ggplot(figsupp8_h3_alpha_72, aes(x=INS, color=treatment))+
  stat_ecdf()
figsupp8_h3_alpha_72 + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log INS") + xlim(0,15) +
  theme(text = element_text(size=25), axis.text.x = element_text(size= 20),axis.text.y = element_text(size= 20))+
  theme(strip.text.y = element_blank())+  facet_grid(~replicate)+
  scale_color_manual(values = color_supp) + theme_bw() +
  theme(text = element_text(size=19), axis.text.x = element_text(size= 16),axis.text.y = element_text(size= 16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")

figsupp8_h3_alpha_36<-ggplot(figsupp8_h3_alpha_36, aes(x=INS, color=treatment))+
  stat_ecdf()
figsupp8_h3_alpha_36 + geom_line(aes(y = 1 - ..y..), stat='ecdf', size=1.2) + ylim(0,.17)+labs(x = "Log INS") + xlim(0,15) +
  theme(text = element_text(size=25), axis.text.x = element_text(size= 20),axis.text.y = element_text(size= 20))+
  theme(strip.text.y = element_blank())+  facet_grid(~replicate)+
  scale_color_manual(values = color_supp) + theme_bw() +
  theme(text = element_text(size=19), axis.text.x = element_text(size= 16),axis.text.y = element_text(size= 16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")

#BETA CELLS
 
ggplot(figsupp8_h3_beta_36, aes(x=INS, color=treatment))+ stat_ecdf(size=1.2)+ 
  labs(x = "Log Ins1") + xlim(5.5,13.5) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_supp) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")

ggplot(figsupp8_h3_beta_72, aes(x=INS, color=treatment))+ stat_ecdf(size=1.2)+ 
  labs(x = "Log Ins1") + xlim(5.5,13.5) +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(strip.text.y = element_blank())+ 
  scale_color_manual(values = color_supp) + theme_bw() +
  theme(text = element_text(size=27), axis.text.x = element_text(size= 22),axis.text.y = element_text(size= 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="none")

#invert outside of R

#############################################################
# Batch - corrected t-SNE 
###############################################################
library(splatter)
library(scater)
#Use this fixed version
set.seed(1)
p = newSplatParams(seed = 1234,
		   batchCells = rep(200, 4), 
		   batch.facLoc = 0.2, 
		   batch.facScale = 0.1,
		   group.prob = c(0.1, 0.2, 0.3, 0.4),
		   de.facLoc = 0.1,
		   de.facScale = 0.4)
sim = splatSimulate(p, method = "groups", verbose = FALSE)
sim = normalize(sim)
str(sim)
X =t(counts(sim))
X = scale(X,center = T,scale = F)
zBatch = model.matrix(~-1+factor(colData(sim)$Batch))
zCELL = model.matrix(~-1+factor(colData(sim)$Group))
#write.table(data.frame(X), file = "X.txt",row.names = F,col.names = F)
#write.table(data.frame(zBatch), file = "z.txt", row.names = F,col.names = F)
#write.table(data.frame(zCELL), file = "zcell.txt", row.names = F,col.names = F)

#+++++++++++++++++++++++++++++++++++++++++++++
#Takes some minutes
#+++++++++++++++++++++++++++++++++++++++++++++
# system("julia ./run_sim.jl ./X.txt ./z.txt")
# system("rm ./X.txt ./z.txt")
#+++++++++++++++++++++++++++++++++++++++++++++
rm(X);gc()
zB = match(colData(sim)$Batch, unique(colData(sim)$Batch)) ## batch
zG = match(colData(sim)$Group, unique(colData(sim)$Group)) ## group/cell-type


Y1     =  read.table("Y50U.txt")
Y1$ADJ =  "UNADJUSTED"
Y1$k   =  ""
Y1$Batches =  as.factor(zB)
Y1$CellType  =  as.factor(zG)

Y2     =  read.table("./tsne_Fastmnn.txt")
Y2$ADJ =  "ADJUSTED"
Y2$k   =  "MNN"
Y2$Batches =  as.factor(zB)
Y2$CellType  =  as.factor(zG)

Y3     =  read.table("./tsne_harmony.txt")
Y3$ADJ =  "ADJUSTED"
Y3$k   =  "HARMONY"
Y3$Batches =  as.factor(zB)
Y3$CellType  =  as.factor(zG)

Y4     =  read.table("./Y10C.txt")
Y4$ADJ =  "ADJUSTED"
Y4$k   =  "BCTSNE"
Y4$Batches =  as.factor(zB)
Y4$CellType  =  as.factor(zG)


cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df = rbind(Y2,Y3,Y4)

df$ADJ = ordered(df$ADJ, levels = c("UNADJUSTED", "ADJUSTED"))
pl_adj = ggplot(df) + geom_point(aes(V1,V2, color=CellType, shape = Batches ),size=1) + 
	facet_wrap(~k, scales = 'free')+
	theme_bw()+
  theme(legend.position="bottom", legend.box = "horizontal",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title= element_blank(),
       legend.title = element_text(size = 40),
	legend.text = element_text(size = 40),
	legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
	strip.text=element_text(size=40,face='bold'))+ 
  guides(color = guide_legend(override.aes = list(size = 10)), shape = guide_legend(override.aes = list(size = 10))) + 
  scale_color_manual(values=(cbbPalette[1:5]))

pl_u = ggplot(Y1) + geom_point(aes(V1,V2, color=CellType, shape = Batches ),size=1) + 
	facet_wrap(~ADJ, scales = 'free')+
	theme_bw()+
  theme(legend.position="none", legend.box = "horizontal",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title= element_blank(),
       legend.title = element_text(size = 40),
	legend.text = element_text(size = 40),
	legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
	strip.text=element_text(size=40,face='bold'))+ 
  guides(color = guide_legend(override.aes = list(size = 10)), shape = guide_legend(override.aes = list(size = 10))) + 
  scale_color_manual(values=(cbbPalette[1:5]))
pl_u

(ly = matrix(c(NA,1,NA,2,2,2),byrow = T,ncol=3))
pl_comb = gridExtra::grid.arrange(pl_u,pl_adj,layout_matrix = ly)

ggsave(pl_comb, filename = paste0(path,"sim.png"),width = 20,height = 13)


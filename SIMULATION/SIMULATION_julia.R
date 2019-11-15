############################################################
# Batch - corrected t-SNE 
# This script reproduces the results of the paper
###############################################################
library(splatter)
library(scater)
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
cor(zBatch)

# Save the data to disk

write.table(data.frame(X), file = "X.txt",row.names = F,col.names = F)
write.table(data.frame(zBatch), file = "z.txt", row.names = F,col.names = F)

#+++++++++++++++++++++++++++++++++++++++++++++
#Takes a while
#+++++++++++++++++++++++++++++++++++++++++++++
system("julia ./run_sim.jl ./X.txt ./z.txt")
# Free space
system("rm ./X.txt ./z.txt")
#+++++++++++++++++++++++++++++++++++++++++++++
rm(X);gc()
zB = match(colData(sim)$Batch, unique(colData(sim)$Batch)) ## batch
zG = match(colData(sim)$Group, unique(colData(sim)$Group)) ## group/cell-type


Y30C     =  read.table("Y30C.txt")
Y30C$ADJ =  "ADJUSTED"
Y30C$k   =  "K = 30"
Y30C$Batches =  as.factor(zB)
Y30C$CellType  =  as.factor(zG)

Y30U     =  read.table("Y30U.txt")
Y30U$ADJ =  "UNADJUSTED"
Y30U$k   =  "K = 30"
Y30U$Batches =  as.factor(zB)
Y30U$CellType  =  as.factor(zG)


Y50C     =  read.table("Y50C.txt")
Y50C$ADJ =  "ADJUSTED"
Y50C$k   =  "K = 50"
Y50C$Batches =  as.factor(zB)
Y50C$CellType  =  as.factor(zG)

Y50U     =  read.table("Y50U.txt")
Y50U$ADJ =  "UNADJUSTED"
Y50U$k   =  "K = 50"
Y50U$Batches =  as.factor(zB)
Y50U$CellType  =  as.factor(zG)

Y100U     =  read.table("Y10U.txt")
Y100U$ADJ =  "UNADJUSTED"
Y100U$k   =  "K = 10"
Y100U$Batches =  as.factor(zB)
Y100U$CellType  =  as.factor(zG)

Y100C     =  read.table("Y10C.txt")
Y100C$ADJ =  "ADJUSTED"
Y100C$k   =  "K = 10"
Y100C$Batches =  as.factor(zB)
Y100C$CellType  =  as.factor(zG)

cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df = rbind(Y30C,Y30U,Y50C,Y50U, Y100C, Y100U)
require(ggplot2)
df$ADJ = ordered(df$ADJ, levels = c("UNADJUSTED", "ADJUSTED"))
ggplot(df) + geom_point(aes(V1,V2, color=CellType, shape = Batches ),size=1) + 
	facet_wrap(~ADJ+k, scales = 'free')+
	theme_bw()+
  theme(legend.position="bottom", legend.box = "horizontal",
	axis.text = element_blank(),
	axis.ticks = element_blank(),
	axis.title= element_blank(),
       legend.title = element_text(size = 25),
	legend.text = element_text(size = 25),
	legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
	strip.text=element_text(size=20,face='bold'))+ 
  guides(color = guide_legend(override.aes = list(size = 10)), shape = guide_legend(override.aes = list(size = 10))) + 
  scale_color_manual(values=(cbbPalette[1:5]))

path="./"
ggsave(filename = paste0(path,"sim.pdf"),width = 20,height = 11)

## CLEAN directory 
system("rm *.txt")
## Print session info for the readme
sessionInfo()

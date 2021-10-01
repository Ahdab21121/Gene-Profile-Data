EM_BMLN=read.table("C:/Users/ahmed/OneDrive/Documents/R/EM_BMLN.csv",header=TRUE,row.names = 1, sep="\t")
DE_BMLN=read.table("C:/Users/ahmed/OneDrive/Documents/R/DE_BMLN.csv",header=TRUE,row.names = 1, sep="\t")
EM_BMPB=read.table("C:/Users/ahmed/OneDrive/Documents/R/EM_BMPB.csv",header=TRUE,row.names = 1, sep="\t")
DE_BMPB=read.table("C:/Users/ahmed/OneDrive/Documents/R/DE_BMPB.csv",header=TRUE,row.names = 1, sep="\t")
DE_PBLN=read.table("C:/Users/ahmed/OneDrive/Documents/R/DE_PBLN.csv",header=TRUE,row.names = 1, sep="\t")
EM_PBLN=read.table("C:/Users/ahmed/OneDrive/Documents/R/EM_PBLN.csv",header=TRUE,row.names = 1, sep="\t")
sig_genes=read.table("C:/Users/ahmed/OneDrive/Documents/R/sig genes.csv",header=TRUE,row.names = 1, sep="\t")
#Create sig genes vector
row.names(sig.genes)=sig.genes[,"SYMBOL"]
sig_genes=row.names(sig.genes)
#merge into universal EM 
EM=merge(EM_BMLN,EM_BMPB, by.x = 0, by.y=0)
row.names(EM)= EM[,1]
EM=EM[2:37]

write.csv(EM, file ="C:/Users/ahmed/OneDrive/Documents/R/EMALL1.csv")
#import tidied up BG
#Merging BG and EM for PCA
master_temp=merge(EM,BG_tidy,by.x = 0,by.y = 1)
master_BMLN=merge(master_temp,DE_BMLN, by.x = 1,by.y = 0)
master_BMPB=merge(master_temp,DE_BMPB, by.x = 1,by.y = 0)
master_PBLN=merge(master_temp,DE_PBLN, by.x = 1,by.y = 0)

master_BMLN$sig=as.factor(master_BMLN$P.ADJ<=0.05 )
master_BMPB$sig=as.factor(master_BMPB$P.ADJ<=0.05 )
master_PBLN$sig=as.factor(master_PBLN$P.ADJ<=0.05 )

master_BMLN_sig=subset(master_BMLN, sig==TRUE)
master_BMPB_sig=subset(master_BMPB, sig==TRUE)
master_PBLN_sig=subset(master_PBLN, sig==TRUE)
#Creating sig gene list for GSEA
row.names(master_temp)=master_temp[,"SYMBOL"]
master_temp=master_temp[sig_genes,]
master_temp=na.omit(master_temp)
row.names(master_temp)=master_temp[,1]
master_temp=master_temp[2:37]


write.csv(master_temp, file ="C:/Users/ahmed/OneDrive/Documents/R/GSEA of sig.csv")

#Symbols as row names and leaving only expression values
row.names(master_temp)=master_temp[,"SYMBOL"]
row.names(master_BMLN)=master_BMLN[,"SYMBOL"]
row.names(master_BMPB)=master_BMPB[,"SYMBOL"]
row.names(master_PBLN)=master_PBLN[,"SYMBOL"]
em_symbols=master_temp
#renaming columns to micro-environment 
#get the sig genes from the vogel dataset
sample_groups = c(rep("marrow",12), rep("node",12), rep("blood",12))
names(em_symbols)=sample_groups
em_symbols=na.omit(em_symbols)
em_scaled=na.omit(data.frame(t(scale(t(em_symbols)))))
em_scaled=em_scaled[sig_genes,]
em_scaled=na.omit(em_scaled)
em_scaled_marrow=em_scaled[1:12]
em_scaled_node=em_scaled[13:24]
em_scaled_blood=em_scaled[25:36]
master_BMLN=master_BMLN[sig_genes,]
master_BMPB=master_BMPB[sig_genes,]
master_PBLN=master_PBLN[sig_genes,]

#Make A PCA plot
#Convert to numeric matrix
# Take a look at expression values, very different from tutorial...
as.matrix(sapply(em_scaled,as.numeric)) #had to use "em_symbols =" otherwise i get an error that its a list when i try to only use finite numbers 
pca=prcomp(t(em_scaled))
pca_coordinates= data.frame(pca$x)
#Plot the PCA
library(ggplot2)
pca= ggplot(pca_coordinates,aes(x=PC1,y=PC2,colour=sample_groups))+
  geom_point(size=4)+ geom_text(aes(label=sample_groups), size=0, position= position_nudge(y=1))+
  scale_colour_manual(values=c("Yellow","Red","Blue"), labels=c("blood","marrow","node"))+ 
  labs(title= "PCA plot for PCA 1 vs PCA2", x="PC1",y="PC2", color="SAMPLE GROUP")+theme_bw()
pca
#Add variance 
vars = apply(pca$x, 2, var)#both PCA variances are exactly the same ... should that happen ?
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100
x_axis_label = paste("PC1", " (",prop_x,"%)",sep="")
y_axis_label = paste("PC2", " (",prop_y,"%)",sep="")
#replot with variance
pca1= ggplot(pca_coordinates,aes(x=PC1,y=PC2,colour=sample_groups))+
  geom_point(size=4)+ # Not sure which position nudge i should use 
  scale_colour_manual(values=c("Yellow","Red","Blue"), labels=c("blood","marrow","node"))+ 
  labs(title= "PCA plot for Blood, Lymph Nodes & Bone Marrow in CLL", x=x_axis_label,y=y_axis_label, color="SAMPLE GROUP")+theme_bw()
pca1
png("C:/Users/ahmed/OneDrive/Documents/R/PCA of gene signatures.png", height = 1300, width=1500)
print(pca1)
dev.off()
#Installing pathway analysis packages
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("RDAVIDWebService")


#Pathway Analysis
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)


master_temp=master_temp[sig_genes,]
master_temp=na.omit(master_temp)
master_BMLN=na.omit(master_BMLN)
master_BMPB=na.omit(master_BMPB)
master_PBLN=na.omit(master_PBLN)
write.csv(master_temp, file ="C:/Users/ahmed/OneDrive/Documents/R/DAVID gene list.csv")

write.csv(master_BMLN, file ="C:/Users/ahmed/OneDrive/Documents/R/BMLN KEGG list.csv")

#GSEA
original_gene_list = master_PBLN$LOG2FOLD
names(original_gene_list) = row.names(master_PBLN) 
gene_list = na.omit(original_gene_list)   #remove the NAs
gene_list = sort(gene_list, decreasing = TRUE)

gse = gseGO(gene_list,ont ="BP", keyType = "SYMBOL", nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, OrgDb = org.Hs.eg.db,pAdjustMethod = "none")
gse

dotplot(gse, showCategory=10, split=".sign") 
ridgeplot(gse)
emapplot(gse) 
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, colorEdge = TRUE, showCategory = 3) 
barplot(gse)


# Creating a heatmap
install.packages("reshape2")
library(reshape2)
library(ggplot2)
em_sig=em_scaled_node[1:41,]
hm.matrix=melt(em_sig)
hm.matrix= as.matrix(em_sig)
hm.matrix=melt(hm.matrix)
colours=c("blue","red")
heatmap= ggplot(hm.matrix, aes(x=Var2,y=Var1, fill=value))+
  geom_tile()+
  scale_fill_gradientn(colours = colorRampPalette(colours)(100))+
  theme(axis.text.x=element_text(angle=45,hjust=1)) 
heatmap
#Clustering the heatmap 
install.packages("amap")
library(amap)
#Y axis
hm.matrix=as.matrix(em_sig)
y.dist=Dist(hm.matrix, method="spearman")
y.cluster= hclust(y.dist, method="average")
y.dd= as.dendrogram(y.cluster)
y.dd.reorder= reorder(y.dd,0,FUN="average")
y.order=order.dendrogram(y.dd.reorder)
hm.matrix_clustered= hm.matrix[y.order,]
hm.matrix_clustered= melt(hm.matrix_clustered)
#X axis
hm.matrix=as.matrix(t(hm.matrix))
x.dist= Dist(hm.matrix, method = "spearman")
x.cluster= hclust(x.dist, method = "average")
x.dd= as.dendrogram(x.cluster)
x.dd.reorder=reorder(x.dd,0, FUN="average")
x.order=order.dendrogram(x.dd.reorder)
hm.matrix_clustered=hm.matrix[x.order,]
hm.matrix_clustered=melt(hm.matrix_clustered)
heatmap_sig = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours= colorRampPalette(colours)(100))+
  ylab("genes")+xlab("samples")+theme(axis.ticks=element_blank(), legend.title = element_blank(),legend.spacing.x = unit(1.0, 'cm'), legend.text = element_text(size = 18))+
  theme(axis.text.x=element_text(angle=45,hjust=1))+theme(plot.title = element_text(face = "bold", size = 30))+
  theme(axis.text.y = element_text(size=18))+theme(axis.text.x = element_blank(),axis.title.x = element_text(size = 28),axis.title.y = element_text(size = 28),legend.key.size = unit(4,'line'))+labs(title = "HeatMap of FOXO1 Gene Signature in node Samples")
heatmap_sig
png("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/heatmap_yclustnode.png", height = 800, width=1800)
print(heatmap_sig)
dev.off()

foxosig= em_scaled["PIK3CA",]
foxosig1=em_scaled[""]

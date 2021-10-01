BG_KMH2=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/BG_KMH2.csv",header=TRUE,row.names = 1, sep="\t")
BG_L428=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/BG_L428.csv",header=TRUE,row.names = 1, sep="\t")
BG_L1236=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/BG_L1236.csv",header=TRUE,row.names = 1, sep="\t")
BG_UHO=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/BG_UHO.csv",header=TRUE,row.names = 1, sep="\t")
BG_SUPHD=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/BG_SUPHD.csv",header=TRUE,row.names = 1, sep="\t")


EM_KMH2=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/EM_KMH2.csv",header=TRUE,row.names = 1, sep="\t")
EM_L428=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/EM_L428.csv",header=TRUE,row.names = 1, sep="\t")
EM_L1236=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/EM_L1236.csv",header=TRUE,row.names = 1, sep="\t")
EM_UHO=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/EM_UHO.csv",header=TRUE,row.names = 1, sep="\t")
EM_SUPHD=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/EM_SUPHD.csv",header=TRUE,row.names = 1, sep="\t")

DE_KMH2=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/DE_KMH2.csv",header=TRUE,row.names = 1, sep="\t")
DE_L428=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/DE_L428.csv",header=TRUE,row.names = 1, sep="\t")
DE_L1236=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/DE_L1236.csv",header=TRUE,row.names = 1, sep="\t")
DE_UHO=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/DE_UHO.csv",header=TRUE,row.names = 1, sep="\t")
DE_SUPHD=read.table("C:/Users/ahmed/OneDrive/Documents/Vogel Dataset/DE_SUPHD.csv",header=TRUE,row.names = 1, sep="\t")


#Creating master table for KvL data set
Merge_tempKL=merge(EM_KvL,BG_KvL, by.x=0, by.y=0)
masterKL=merge(Merge_tempKL,DE_KvL,by.x=1,by.y=0)

masterKL=na.omit(masterKL)
row.names(masterKL)=masterKL[,10]

#Creating master table for KvS data set
Merge_tempKS=merge(EM_KvS,BG_KvS, by.x=0, by.y=0)
masterKS=merge(Merge_tempKS,DE_KvS,by.x=1,by.y=0)

#Creating master table for L4 vs U data set
Merge_tempL4U=merge(EM_L4vU,BG_L4vU, by.x=0, by.y=0)
masterL4U=merge(Merge_tempL4U,DE_L4vU,by.x=1,by.y=0)

#Creating master table for all data sets
master1=merge(masterKL,masterKS, by.x = 1,by.y = 1)
master=merge(master1,masterL4U, by.x=1, by.y=1)
master=master[,-26]
row.names(master)=master[,"SYMBOL"]
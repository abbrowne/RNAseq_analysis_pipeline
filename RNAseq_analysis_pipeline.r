##Load necessary libraries and functions
library(DESeq)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(rgl)
library(edgeR)
library(EDASeq)
library(corrplot)
library(WGCNA)
library(limma)
library(sva)
library(variancePartition)
library(statmod)
library(calibrate)
library(doParallel)
library(pca3d)
library(reshape2)


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
}

corr_comparison <- function(input_matrix,input_meta=NEU_meta,corr_input_name=input_name){
centered_input <- t(scale(t(input_matrix),scale=F)) ## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC <- prcomp(centered_input)
topPC<- PC$rotation[,1:5]
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topvar <- varexp[1:5]
colnames(topPC) <- paste("PC_MATS2\n",colnames(topPC)," (",signif(100*topvar[1:5],2),"%)",sep="")
#covariate_df <- data.frame(Dx=as.factor(input_meta$Dx),Seq_batch=as.factor(input_meta$Seq_batch),Sex=as.factor(input_meta$Sex),Diff_batch=as.factor(input_meta$Diff_batch), Individuals=as.factor(input_meta$Individuals),Age=as.numeric(input_meta$Age),Families=as.factor(input_meta$Families),RIN=as.numeric(input_meta$RIN),Total_reads=as.numeric(input_meta$Total_reads), Del_size=as.numeric(input_meta$Del_size), Timepoint=as.numeric(input_meta$Timepoint))
covariate_df <- data.frame(Astrocytes=as.numeric(input_meta$Astrocytes), Endothelial_cells=as.numeric(input_meta$Endothelial_cells), Fetal_quiescent_neurons=as.numeric(input_meta$Fetal_quiescent_neurons), Fetal_replicating_neurons=as.numeric(input_meta$Fetal_replicating_neurons), Microglia=as.numeric(input_meta$Microglia), Neurons=as.numeric(input_meta$Neurons), Oligodendrocytes=as.numeric(input_meta$Oligodendrocytes), Dx=as.factor(input_meta$Dx),Sex=as.factor(input_meta$Sex),Diff_batch=as.factor(input_meta$Diff_batch), Individuals=as.factor(input_meta$Individuals),Age=as.numeric(input_meta$Age),Families=as.factor(input_meta$Families),RIN=as.numeric(input_meta$RIN),Total_reads=as.numeric(input_meta$Total_reads), Del_size=as.numeric(input_meta$Del_size), Timepoint=as.numeric(input_meta$Timepoint))
col_cond <- input_meta$sample_colors ## colors
topPC2 <- topPC
colnames(topPC2) <- c("PC1","PC2","PC3","PC4","PC5")
form <- ~ PC1 + PC2 + PC3 + PC4 + PC5 + Dx + Timepoint + Sex + Diff_batch + Total_reads + Age + Families + RIN + Individuals + Del_size + Astrocytes + Endothelial_cells + Fetal_quiescent_neurons + Fetal_replicating_neurons + Microglia + Neurons + Oligodendrocytes
cor_matrix <- canCorPairs(form,cbind(topPC2,covariate_df))
hM <- format(round(cor_matrix, 2))
pdf_name <- paste(corr_input_name,"PCA to covariate correlations.pdf",sep="_")
pdf(pdf_name)
plotCorrMatrix(cor_matrix,sort=FALSE,cellnote=hM,notecex=0.5,notecol="black")
pairs(cbind(topPC,covariate_df),col=col_cond,pch=19,cex.labels=0.5,upper.panel = panel.cor,main=paste(corr_input_name,"PCA to covariate correlations",sep=" "))
pca3d(topPC, show.ellipses=TRUE,ellipse.ci=0.95, show.plane=FALSE)
#PCA correlation plot
dev.off()
}

norm_comparison <- function(norm_input_name,comparison_set,comparison_matrix_t,comp_meta,comparison_matrix,name=NULL){
pdf_name=paste(norm_input_name,"norm_comparison.pdf",sep="_")
pdf(pdf_name)
boxplot(comparison_set,col=comp_meta$sample_colors, las=2, cex=0.5, cex.axis=0.5, main=paste(norm_input_name,name,sep="_"), outline=FALSE, ylab="log Expression",ylim=c(-2,4))
legend(10,4,fill = c("red","blue","pink","lightblue"), cex=0.6, legend = c("CTRL_F","CTRL_M","PMS_F","PMS_M"), bty='n', horiz=T,x.intersp=0.1)
hist(log(counts(comparison_set)), main=paste(norm_input_name,name,sep="_"), xlab="log Expression")
plotRLE(comparison_set,ylim=c(-1,1), cex.axis=0.5, las=2,col=comp_meta$sample_colors, outline=FALSE, ylim="RLE", main=paste(norm_input_name,name,sep="_"))
legend(10,1,fill = c("red","blue","pink","lightblue"), cex=0.6, legend = c("CTRL_F","CTRL_M","PMS_F","PMS_M"), bty='n', horiz=T,x.intersp=0.1)
plotMDS(comparison_matrix,cex=0.8,col=comp_meta$sample_colors,main=paste(norm_input_name,name,sep="_"),ylim=c(-3, 3), xlim=c(-4, 4))
legend(-3,3,fill = c("red","blue","pink","lightblue"), cex=0.6, legend = c("CTRL_F","CTRL_M","PMS_F","PMS_M"), bty='n', horiz=T,x.intersp=0.1)
plotMDS(comparison_matrix,cex=0.8,col=comp_meta$sample_colors,main=paste(norm_input_name,name,sep="_"),ylim=c(-3, 3), xlim=c(-4, 4))
legend(-3,3,fill = c("red","blue","pink","lightblue"), cex=0.6, legend = c("CTRL_F","CTRL_M","PMS_F","PMS_M"), bty='n', horiz=T,x.intersp=0.1)
plot(hclust(dist(comparison_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
dev.off()
}

full_pipeline <- function(input_name="FULL_NEU",small_only=FALSE,run_vp=FALSE,pipeline_dir="YOUR_WORKING_DIR", voom_design="~Diff_batch", lmreg_design="Dx+RIN+Age+Diff_batch+Families+Sex+Timepoint+Endothelial_cells",count_file="FULL_NEU_COUNTS.csv", meta_file="FULL_NEU_META.csv", remove_ERCCs=TRUE, exclude_samples=TRUE,NEU_groups_to_use="A"){

	setwd(pipeline_dir)

	NEU_counts <- NULL
	NEU_counts <- read.csv(count_file,header=TRUE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
	NEU_counts <- NEU_counts[,order(colnames(NEU_counts))]
	NEU_counts <- NEU_counts[order(rownames(NEU_counts)),]
	
	##Option to remove ERCCs
	if(remove_ERCCs){
		print(NEU_counts[56632:nrow(NEU_counts),1:5])
		NEU_counts <- NEU_counts[1:56632,]
	}
	
	##Meta loading and sorting
	NEU_meta <- NULL
	NEU_meta <- read.csv(meta_file,header=TRUE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
	NEU_meta <- NEU_meta[order(rownames(NEU_meta)),]
	
	##Remove outliers
	if(exclude_samples){
		NEU_counts <- NEU_counts[,NEU_meta$exclude=="no"]
		NEU_meta <- NEU_meta[NEU_meta$exclude=="no",]
	}
	
	##Subset to small deletion samples
	if(small_only){
		NEU_counts <- NEU_counts[,NEU_meta$Del_class==("Small") | NEU_meta$Del_class==("No_del")]
		NEU_meta <- NEU_meta[NEU_meta$Del_class==("Small") | NEU_meta$Del_class==("No_del"),]
	}
	
	##Subset to group of interest
	groups_to_use <- strsplit(NEU_groups_to_use,split="")
	groups_to_use <- groups_to_use[[1]]
	NEU_counts <- NEU_counts[,NEU_meta$Seq_batch %in% groups_to_use]
	NEU_meta <- NEU_meta[NEU_meta$Seq_batch %in% groups_to_use,]
	
	##Print sample summary
	print(strsplit(paste('The input matrix contains expression counts for ',nrow(NEU_counts),' genes in a total of ',ncol(NEU_counts),' samples\n',sep=""),"\n")[[1]])
	print(strsplit(paste('A total of ',sum(colnames(NEU_counts) %in% rownames(NEU_meta)),' out of ',ncol(NEU_counts),' total samples were found in the included metadata',sep=""),"\n")[[1]])
	print(strsplit(paste('The metadata contains a total of ',nrow(NEU_meta),' samples',sep=""),"\n")[[1]])
	print(strsplit(paste('The metadata consists of ',sum(NEU_meta$Dx == "PMS"),' samples derived from PMS patients and ',sum(NEU_meta$Dx == "CTRL"),' samples derived from their healthy siblings',sep=""),"\n")[[1]])
	print(strsplit(paste('A total of ',sum(NEU_meta$Sex == "M"),' samples were derived from males and ',sum(NEU_meta$Sex == "F"),' were from females',sep=""),"\n")[[1]])
	print(strsplit(paste('The samples came from ',length(unique(NEU_meta$Individuals)),' individuals from a total of ',length(unique(NEU_meta$Families)),' families',sep=""),"\n")[[1]])
	print(strsplit(paste('A total of ',sum(NEU_meta$Del_class == "Small"),' samples with a small deletion and ',sum(NEU_meta$Del_class == "Large"),' samples with a large deletion',sep=""),"\n")[[1]])
	print(strsplit(paste('There are ',sum(NEU_meta[!duplicated(NEU_meta$Individuals),"Del_class"] == "No_del"),' control siblings, ',sum(NEU_meta[!duplicated(NEU_meta$Individuals),"Del_class"] == "Small"),' patients with small deletions, and ',sum(NEU_meta[!duplicated(NEU_meta$Individuals),"Del_class"] == "Large"),' patients with large deletions',sep=""),"\n")[[1]])
	print(strsplit(paste('A total of ',sum(NEU_meta$Del_class == "Small"),' samples have a small deletion and ',sum(NEU_meta$Del_class == "Large"),' samples have a large deletion',sep=""),"\n")[[1]])
	
	##Convert covariates to factors
	NEU_meta$Dx = as.factor(NEU_meta$Dx)
	NEU_meta$Sex = as.factor(NEU_meta$Sex)
	NEU_meta$Diff_batch = as.factor(NEU_meta$Diff_batch)
	NEU_meta$Seq_batch = as.factor(NEU_meta$Seq_batch)
	NEU_meta$Families = as.factor(NEU_meta$Families)
	NEU_meta$Individuals = as.factor(NEU_meta$Individuals)
	NEU_meta$Del_class = as.factor(NEU_meta$Del_class)
	
	##Filtering of lowly expressed genes
	filtered_counts <- NULL
	DGE_counts <- NULL
	DGE_counts <- DGEList(counts=as.matrix(NEU_counts), genes=rownames(NEU_counts))
	tokeep <- NULL
	tokeep <- rowSums(cpm(DGE_counts)>1) >= ncol(NEU_counts)/2 # require 1 cpm in at least 10 samples
	filtered_counts <- DGE_counts[tokeep,]
	print(dim (filtered_counts))
	print(strsplit(paste('A total of ',nrow(filtered_counts),' out of ',nrow(NEU_counts),' genes remain after filtering for those with a CPM > 1 in at least half of the samples',sep=""),"\n")[[1]])
	
	##Voom normalization
	voom_design1 <- model.matrix(as.formula(voom_design),NEU_meta)
	set <- newSeqExpressionSet(as.matrix(filtered_counts),phenoData=data.frame(voom_design1,row.names=colnames(filtered_counts)))
	filtered_t <- t(filtered_counts$counts)
	NF_counts <- calcNormFactors(filtered_counts)
	pdf(file=paste(input_name,"_voom_disp.pdf",sep=""))
	VST <- voom(NF_counts, design = voom_design1, plot=TRUE)
	dev.off()
	voom_matrix1 <- NULL
	voom_matrix1 <- cbind(VST$genes, VST$E)
	rownames(voom_matrix1) <- voom_matrix1[,1]
	voom_matrix1 <- voom_matrix1[,2:ncol(voom_matrix1)]
	voom_matrix_t <- t(voom_matrix1)
	
	##Code for running variance partition to identify which covariates explain the most variance
	if(run_vp){
		form <- ~ (1|Dx) + (1|Sex) + Age + (1|Families) + (1|Individuals) + (1|Diff_batch) + RIN + Timepoint + Endothelial_cells
		geneExpr <- as.matrix(voom_matrix1)
		info <- NEU_meta
		cl <- makeCluster(4)
		registerDoParallel(cl)
		varPart <- fitExtractVarPartModel( geneExpr, form, info )
		vp <- sortCols( varPart )
		pdf("VP_test5.pdf")
		plotVarPart( vp )
		dev.off()
	}

	##Create voom expression set
	set_voom <- NULL
	set_voom <- newSeqExpressionSet(as.matrix(voom_matrix1),phenoData=data.frame(voom_design1,row.names=colnames(voom_matrix1)))
	corr_comparison(input_matrix=voom_matrix1,input_meta=NEU_meta,corr_input_name=paste(input_name,sep="_"))
	norm_comparison(norm_input_name=paste(input_name,sep="_"),name="voom",comparison_set=set_voom,comparison_matrix_t=voom_matrix_t,comparison_matrix=voom_matrix1,comp_meta=NEU_meta)
	
	##Correction for duplicate samples using duplicateCorrelation
	corfit1 <- duplicateCorrelation(voom_matrix1,design=voom_design1,block=NEU_meta$Individuals)
	pdf(file=paste(input_name,"corDup_disp.pdf",sep="_"))
	voom_corDup1 <- voom(NF_counts,design=voom_design1,block=NEU_meta$Individuals,correlation=corfit1$consensus, plot=TRUE)
	dev.off()
	corDup_matrix1 <- NULL
	corDup_matrix1 <- cbind(voom_corDup1$genes,voom_corDup1$E)
	rownames(corDup_matrix1) <- corDup_matrix1[,1]
	corDup_matrix1 <- corDup_matrix1[,2:ncol(corDup_matrix1)]
	corfit2 <- duplicateCorrelation(corDup_matrix1,voom_design1, block=NEU_meta$Individuals)
	corr_comparison(input_matrix=corDup_matrix1,input_meta=NEU_meta,corr_input_name=paste(input_name,"corDup1_corr",sep="_"))
	corDup_matrix_t <- t(corDup_matrix1)
	set_corDup <- NULL
	set_corDup <- newSeqExpressionSet(as.matrix(corDup_matrix1),phenoData=data.frame(voom_design1,row.names=colnames(corDup_matrix1)))
	
	##Linear modeling of covariates
	lmreg_cols <- length(strsplit(lmreg_design,"\\+")[[1]])+2
	lmreg_design3 <- paste("~0+",lmreg_design,sep="")
	lmreg_design1 <- model.matrix(as.formula(lmreg_design3),NEU_meta,ref="CTRL")
	lmFit_design <- NULL
	lmFit_design <- lmreg_design1
	lmFit_model <- lmFit(as.matrix(corDup_matrix1), lmFit_design, block=NEU_meta$Individuals,correlation=corfit2$consensus)
	
	##Differential expression analysis
	contr_matrix <- makeContrasts(PMS_vs_CTRL=DxPMS - DxCTRL,levels=lmFit_design)
    contr_fit <- NULL
	contr_fit <- contrasts.fit(lmFit_model,contr_matrix)	
	eFit <- eBayes(contr_fit)
	lmFit_results <- NULL
    lmFit_results <- topTable(eFit,coef=1, adjust="BH", number=nrow(corDup_matrix1))
	write.csv(lmFit_results,file=paste(input_name,"_regressed_PMS_DEX.csv",sep=""),quote=FALSE)
	

}

full_pipeline(pipeline_dir='.',NEU_groups_to_use="ABC")




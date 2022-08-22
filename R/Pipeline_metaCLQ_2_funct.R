

require(ggplot2)
require(reshape2)
require(RColorBrewer)
require(metafor)
require(ggrepel)



MetaCLQ <- function(ID,nTile,nNeig,bandwidth,metadata,dir_output1){

	dir.create(file.path(dir_output1,paste0("nTile",nTile,"_","nNeig",nNeig,"_","Band",bandwidth),"MetaCLQ",ID),recursive=TRUE)
	dir_output <- file.path(dir_output1,paste0("nTile",nTile,"_","nNeig",nNeig,"_","Band",bandwidth),"MetaCLQ",ID)
	dir_input <- file.path(dir_output1,paste0("nTile",nTile,"_","nNeig",nNeig,"_","Band",bandwidth))


	###################### Reading Data

	allIDs <- list.files(file.path(dir_input))
	allIDs <- allIDs[grep("CLQ_realANDresampled",allIDs)]
	CLQsub_all <- c()
	for(i in allIDs){
		CLQsub <- read.table(file.path(dir_input,i),sep="\t",row.names=1,header=TRUE)
		CLQsub_all <- rbind(CLQsub_all,CLQsub)
	}
	CLQsub_all$Pair <- paste0(CLQsub_all[,3],"-",CLQsub_all[,4])
	CLQsub_all$Sample <- gsub(".csv","",CLQsub_all$Sample)

	sample_GR1 <- metadata[which(metadata$Group == 0),]$Sample
	sample_GR2 <- metadata[which(metadata$Group == 1),]$Sample


	################## Running METACLQ
	###### Idenfity binomial distributions and make them unimodal or zero colocalization
	screen=function(y){
	 ss=rep(NA, 3)
	 for(b in 1:3)
	     {ss[b]=kmeans(y, centers=b, iter.max=500)$tot.withinss}
	 reduction=c(ss[2]/ss[1], ss[3]/ss[2])
	 ynew=y
	 if(reduction[1]<reduction[2]*0.64)
	    {fit=kmeans(y, centers=2, iter.max=500)
	     n1=sum(fit$cluster==1)
	     n2=sum(fit$cluster==2)
	     if(n1>n2)   ynew=y[fit$cluster==1]
	     if(n2>n1)   ynew=y[fit$cluster==2]
	     }

	  return(list(CLQs.new=ynew,  flag=(reduction[1]<reduction[2]*0.64)))
	 }

	###### MetaCLQ

	data_all_1 <- c()
	data_all_2 <- c()
	cellid_1 <- c()
	cellid_2 <- c()
	bimod <- c()
	bimod_aftercorrec <- c()
	unimod <- c()
	alldist <- c()
	for(pairCell in unique(CLQsub_all$Pair)){
		CLQsubjPair <- CLQsub_all[CLQsub_all$Pair %in% pairCell,]

		sigma_GR1 <- c()
		sigma_GR2 <- c()
		CLQ_GR1 <- c()
		CLQ_GR2 <- c()
		CLQs.new_GR1 <- c()
		CLQs.new_GR2 <- c()
		for(sampl_1 in sample_GR1){
			CLQsubj_1 <- CLQsubjPair[CLQsubjPair$Sample %in% sampl_1,]
			CLQ_resampled_1 <- CLQsubj_1[-which(CLQsubj_1$Iteration==0),]
			if(nrow(CLQ_resampled_1) == 0){next}

			if(sum(CLQ_resampled_1[,5]) != 0 & sd(CLQ_resampled_1[,5]) > 0.1 & dim(table(CLQ_resampled_1[,5]))>2){
				fit=screen(CLQ_resampled_1[,5])
				if(sum(fit$CLQs.new) != 0 & sd(fit$CLQs.new) > 0.1 & dim(table(fit$CLQs.new))>2){
					fit2=screen(fit$CLQs.new)
						CLQs.new_GR1 <- rbind(CLQs.new_GR1,data.frame(Sample=sampl_1,CellPair=pairCell,CLQ=fit$CLQs.new,Flag=fit$flag,Flag2=fit2$flag))
					}else{
						CLQs.new_GR1 <- rbind(CLQs.new_GR1,data.frame(Sample=sampl_1,CellPair=pairCell,CLQ=fit$CLQs.new,Flag=fit$flag,Flag2="ZeroColoc"))
				}
			}else{
				CLQs.new_GR1 <- rbind(CLQs.new_GR1,data.frame(Sample=sampl_1,CellPair=pairCell,CLQ=CLQ_resampled_1[,5],Flag="ZeroColoc",Flag2="ZeroColoc"))
			}

			cellid_1 <- c(cellid_1,pairCell)

			sigma1=sd(CLQs.new_GR1$CLQ) #instead of var
			sigma_GR1 <- c(sigma_GR1,sigma1)

			CLQ_real <- CLQsubj_1[which(CLQsubj_1$Iteration==0),]
			CLQ_GR1 <- c(CLQ_GR1,CLQ_real$value)
		}

		for(sampl_2 in sample_GR2){
			CLQsubj_2 <- CLQsubjPair[CLQsubjPair$Sample %in% sampl_2,]
			CLQ_resampled_2 <- CLQsubj_2[-which(CLQsubj_2$Iteration==0),]
			if(nrow(CLQ_resampled_2) == 0){next}

			if(sum(CLQ_resampled_2[,5]) != 0 & sd(CLQ_resampled_2[,5]) > 0.1 & dim(table(CLQ_resampled_2[,5]))>2){
				fit=screen(CLQ_resampled_2[,5])
				if(sum(fit$CLQs.new) != 0 & sd(fit$CLQs.new) > 0.1 & dim(table(fit$CLQs.new))>2){
					fit2=screen(fit$CLQs.new)
						CLQs.new_GR2 <- rbind(CLQs.new_GR2,data.frame(Sample=sampl_2,CellPair=pairCell,CLQ=fit$CLQs.new,Flag=fit$flag,Flag2=fit2$flag))
					}else{
						CLQs.new_GR2 <- rbind(CLQs.new_GR2,data.frame(Sample=sampl_2,CellPair=pairCell,CLQ=fit$CLQs.new,Flag=fit$flag,Flag2="ZeroColoc"))
				}
			}else{
				CLQs.new_GR2 <- rbind(CLQs.new_GR2,data.frame(Sample=sampl_2,CellPair=pairCell,CLQ=CLQ_resampled_2[,5],Flag="ZeroColoc",Flag2="ZeroColoc"))
			}

			cellid_2 <- c(cellid_2,pairCell)

			sigma2=sd(CLQs.new_GR2$CLQ)
			sigma_GR2 <- c(sigma_GR2,sigma2)

			CLQ_real <- CLQsubj_2[which(CLQsubj_2$Iteration==0),]
			CLQ_GR2 <- c(CLQ_GR2,CLQ_real$value)
		}

		bimod <- c(bimod,length(which(CLQs.new_GR1$Flag == TRUE)))
		bimod_aftercorrec <- c(bimod_aftercorrec,length(which(CLQs.new_GR1$Flag2 == TRUE)))
		unimod <- c(unimod,length(which(CLQs.new_GR1$Flag == FALSE)))
		alldist <- rbind(alldist,CLQs.new_GR1,CLQs.new_GR2)

		#if(sum(CLQ_GR1,sigma_GR1) == 0){
		if(sum(CLQ_GR1) == 0){
			paircell_data_1 <- c(0,0,0,0,0,0)
		}else{
			rma_GR1 <- rma.uni(yi=CLQ_GR1, sei=sigma_GR1, method="SJ") #sei: the square root of the sampling variances
			paircell_data_1 <- coef(summary(rma_GR1))
		}

		#if(sum(CLQ_GR2,sigma_GR2) == 0){
		if(sum(CLQ_GR2) == 0){
			paircell_data_2 <- c(0,0,0,0,0,0)
		}else{
			rma_GR2 <- rma.uni(yi=CLQ_GR2, sei=sigma_GR2, method="SJ")
			paircell_data_2 <- coef(summary(rma_GR2))
		}

		data_all_1 <- rbind(data_all_1,paircell_data_1)
		data_all_2 <- rbind(data_all_2,paircell_data_2)

	}

	#print(paste0("Total number of bimodal distributions: ",sum(bimod)))
	#print(paste0("Total number of bimodal distributions after correction: ",sum(bimod_aftercorrec)))
	#print(paste0("Total number of unimodal distributions: ",sum(unimod)))
	#print(paste0("Percentage of unimodal distributions: ",sum(unimod)/sum(bimod+unimod)*100))



	rownames(data_all_1) <- unique(CLQsub_all$Pair)
	rownames(data_all_2) <- unique(CLQsub_all$Pair)
	write.table(data_all_1,file.path(dir_output,"Meta_analysis_GR1.txt"),sep="\t",row.names=TRUE,col.names=NA)
	write.table(data_all_2,file.path(dir_output,"Meta_analysis_GR2.txt"),sep="\t",row.names=TRUE,col.names=NA)



	## Generating zscore
	identical(rownames(data_all_1),rownames(data_all_2))
	zscore <- (data_all_2$estimate - data_all_1$estimate)/sqrt((data_all_2$se)^2 + (data_all_1$se)^2) #Z=(est GR2-  est GR1)/ sqrt((se GR2)^2 +(se GR1)^2)
	#pval <- pnorm(zscore, mean = 0, sd = 1, lower.tail = TRUE)
	pval <- (1-pnorm(abs(zscore), mean = 0, sd = 1, lower.tail = TRUE)) * 2
	FDR <- p.adjust(pval, method="BH")
	Bonferroni <- p.adjust(pval, method="bonferroni")
	dat <- data.frame(CellPair=rownames(data_all_1),Zscore=zscore,Pvalue=pval,FDR=FDR,Bonferroni=Bonferroni)
	dat <- dat[order(-abs(dat$Zscore)),]
	dat <- dat[order(dat$Pvalue),]
	write.table(dat,file.path(dir_output,"Output_Zscores.txt"),sep="\t",row.names=TRUE,col.names=NA)

	dat2 <- dat[-grep("Unknown|general",dat$CellPair),]
	dat2$FDR <- p.adjust(dat2$Pvalue, method="BH")
	write.table(dat2,file.path(dir_output,"Output_Zscores_filtered.txt"),sep="\t",row.names=TRUE,col.names=NA)



	## Scatter plots per condition
	status_GR1 <- unique(metadata[which(metadata$Group == 0),]$Status)
	status_GR2 <- unique(metadata[which(metadata$Group == 1),]$Status)

	zscore1 <- data_all_1$estimate/data_all_1$se
	zscore2 <- data_all_2$estimate/data_all_2$se
	dat2 <- data.frame(CellPair=rownames(data_all_1),GR1=zscore1,GR2=zscore2)

	scatter2 <- ggplot(dat2, aes(x=GR1, y=GR2)) +
	  geom_point() +
		theme_bw() +
	    theme(axis.title=element_text(size=12,face="bold"),
			plot.title = element_text(face = "bold",hjust = 0.5)) +
		ggtitle("Meta-CLQs") +
	  xlab(status_GR1)+
		ylab(status_GR2) +
		geom_label_repel(
					  aes(label = CellPair),
	                  box.padding   = 0.35,
	                  point.padding = 0.5,
	                  segment.color = 'grey50'
					  )
	ggsave(file.path(dir_output,"Zscorecomparison.png"),width=5,height=5)


	identical(rownames(data_all_2),rownames(data_all_1))
	dat <- dat[order(match(dat$CellPair,rownames(data_all_1))),]
	identical(dat$CellPair,rownames(data_all_1))

	FC <- data.frame(CellPair=dat$CellPair,Pvalue=dat$Pvalue,Zscore=dat$Zscore,log2FC=(log2(data_all_2$estimate) - log2(data_all_1$estimate)))
	rownames(FC) <- dat$CellPair
	FC <- FC[!is.infinite(FC$log2F),]
	FC$Differential <- "NO"
	FC$Differential[FC$log2FC > 0.5 & FC$Pvalue < 0.05] <- "UP"
	FC$Differential[FC$log2FC < -0.5 & FC$Pvalue < 0.05] <- "DOWN"
	FC$delabel <- NA
	FC$delabel[FC$Differential != "NO"] <- rownames(FC)[FC$Differential != "NO"]

	  # plot adding up all layers we have seen so far
	  a <- ggplot(data=FC, aes(x=log2FC, y=-log10(Pvalue), col=Differential, label=delabel)) +
	          geom_point() +
	          theme_minimal() +
	          geom_text_repel() +
	          scale_color_manual(values=c("blue", "black","tan2")) +
	          geom_vline(xintercept=c(-0.5, 0.5), linetype="dotted", col="tomato3") +
	          geom_hline(yintercept=-log10(0.05), linetype="dotted", col="tomato3")
	  ggsave(file.path(dir_output,"VolcanoPlot.png"))

  return(list(MetaCLQ_1=data_all_1, MetaCLQ_2 = data_all_2, Zscores=dat2))


}





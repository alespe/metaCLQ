

###################

# Meta-CLQ
# Author: Almudena Espin-Perez, Stanford University

# @ CaseID
# @ nTile
# @ nNeig
# @ bandwidth
# @ dir_files
# @ dir_output


	###################### Hyperparameters

	# nTile =20
	# nNeig =10
	# bandwidth =100

	###################### Reading Data

	# data_codex <- read.table(file.path(dir_files,IDcase),sep=",",header=TRUE)
	# data_codexSub <- data.frame(Unique_id_all=rownames(data_codex),cell_type=data_codex$Final.cell.type,x=data_codex$X,y=data_codex$Y)
	# data_codexSub$cell_type[which(is.na(data_codexSub$cell_type))] <- "Unknown"
	# print(paste("Running metaCLQ for sample",IDcase,"with",nrow(data_codexSub),"number of cells"))

	# dir.create(file.path(dir_output,paste0("nTile",nTile,"_","nNeig",nNeig,"_","Band",bandwidth)),showWarnings = FALSE)
	# dir_output <- file.path(dir_output,paste0("nTile",nTile,"_","nNeig",nNeig,"_","Band",bandwidth))

  # Sample1 <- SpatialBootstap("A",nTile,nNeig,bandwidth,data_codexSub,dir_output) # output BootstraptCLQ and CLQs
  # Sample2 <- SpatialBootstap("B",nTile,nNeig,bandwidth,data_codexSub,dir_output) # output BootstraptCLQ and CLQs

	# metadata <- read.table(metadata_file,sep="\t",header=TRUE)
  # results <- MetaCLQ(ID,nTile,nNeig,bandwidth,metadata,dir_output) # output MetaCLQ group 1, MetaCLQ group 2 and Zscores

######################################################## RUNNING CLQs

require(spdep)
require(ggplot2)
require(reshape2)
require(RColorBrewer)
require(dplyr)
#source("CLQ_functions.R")

SpatialBootstrap <- function(IDcase,nTile,nNeig,bandwidth,data_codexSub,dir_output){
  cat("\n")
  print(paste("Running spatial bootstrap for sample ",IDcase))

	set.seed(1)

	###################### SpatialBootstrap

	## This scripts aims to break the images into tiles, and reassambled them into an image where the location for each tile has changed

	#################### CLQs
	mat_all <- c()
	ncells_all <- c()
	CLQsub_all <- c()

	  ## Indentify the min and max for x and y coordinates
		max_x <- max(data_codexSub$x)
		min_x <- min(data_codexSub$x)
		max_y <- max(data_codexSub$y)
		min_y <- min(data_codexSub$y)

	  ## Calculate the tile size if we divide the image into nTile tiles
		size_x <- round((max_x-min_x)/nTile)
		size_y <- round((max_y-min_y)/nTile)

	  ## Create individual tiles with coordinates comprised into the tile size.
		tile0 <- list()
		for(tilex in 1:nTile){
		## Subset the cells that are within the range for an invididual tile in the x coordinates
		limit_x_min <- min_x + (tilex-1)*size_x
		limit_x_max <- min_x + tilex*size_x
		data_codexSub_coordx <- data_codexSub[which(data_codexSub$x < limit_x_max & data_codexSub$x >= limit_x_min),]

			for(tiley in 1:nTile){
			## Same with y coordinates. data_codexSub_coordy is one individual tile.
			limit_y_min <- min_y + (tiley-1)*size_y
			limit_y_max <- min_y + tiley*size_y
			data_codexSub_coordy <- data_codexSub_coordx[which(data_codexSub_coordx$y < limit_y_max & data_codexSub_coordx$y >= limit_y_min),]

			ncells <- c(IDcase,
			limit_x_min,
			limit_x_max,
			limit_y_min,
			limit_y_max,
			dim(data_codexSub_coordy)[1],
			nrow(data_codexSub))
			ncells_all <- rbind(ncells_all,ncells)

			## Make all coordinates as they would be if it would be the first tile in x and y axes
			new_tile <- data_codexSub_coordy
			new_tile$x <- new_tile$x - (tilex-1)*size_x
			new_tile$y <- new_tile$y - (tiley-1)*size_y

			## Put tiles in a list. The name of the tile in the list corresponds to the location of tiles in the x and y axes in the original image
			tile0[[paste0(tilex,"-",tiley)]] <- new_tile
			}
		}

		##  Divide non empty from empty tiles so in the next loop only non empty get ressambled while tiles with zero cells stays on the same location
		emptyTiles <- tile0[which(lapply(tile0,nrow) == 0)]
		tilenon0 <- tile0[-which(lapply(tile0,nrow) == 0)]

		##  Resampling with replacement, 1000 times
		all_ranges <- c()
		tiles_order <- c()
		Ncells <- nrow(data_codexSub)
		pb = txtProgressBar(min = 0, max = 1000, initial = 0, style = 3)
		for(it in 1:1000){
		  setTxtProgressBar(pb,it)
			tileReshuffle <- list()
			idstiles <- names(tilenon0)
			## Sampling ids with replacement
			tileIDResampl <- sample(idstiles, replace=TRUE)
			tiles_order <- rbind(tiles_order,c(IDcase,tileIDResampl))
			## List all the selected tiles with a random order into the new object tileReshuffle
			for(reshuf in 1:length(tileIDResampl)){
				tileReshuffle[[reshuf]] <- tilenon0[[tileIDResampl[reshuf]]]
			}
			## Assign location regarding position in the image to the tiles from the list (give same name as original list, meaning that 1st element in the list will be the tile in the bottom x and y coords, etc.)
			names(tileReshuffle) <- idstiles

			## Add empty tiles with the real location
			tileReshuffle <- do.call(c, list(emptyTiles, tileReshuffle))

			## Reassemble the image with the selected tiles
			for(i in 1:length(tileReshuffle)){
				## For each tile, adjust the coordinates to their location in the original image according to the reshuffled location
				tileLoc <- names(tileReshuffle)[i]

				addxcoords <- size_x * (as.numeric(gsub("-.*","",tileLoc))-1)
				tileReshuffle[[i]]$x <- tileReshuffle[[i]]$x + addxcoords

				addycoords <- size_y * (as.numeric(gsub(".*-","",tileLoc))-1)
				tileReshuffle[[i]]$y <- tileReshuffle[[i]]$y + addycoords
			}

			## Convert the list into dataframe so it looks like the original dataframe with the correct tile positions
			tileReshuffle_table <- do.call(rbind.data.frame, tileReshuffle)
			Ncells <- c(Ncells,nrow(tileReshuffle_table))

			## Calculate CLQ:

			  ################################### CLQ calculations per tile

			  cell_type_assignment <- tileReshuffle_table$cell_type
			  coords <- data.frame(x=tileReshuffle_table$x,y=tileReshuffle_table$y)
			  row.names(coords) <- rownames(tileReshuffle_table)

			  ##################### Identifying Neighbors
			  nb_list <- KNN_neighbors(as.matrix(coords),number_of_neighbors=nNeig)

			  for(i in length(nb_list)){
				neigX <- nb_list[[i]]
				neigXCoor <- coords[neigX,]
				xrange <- max(neigXCoor$x)-min(neigXCoor$x)
				yrange <- max(neigXCoor$y)-min(neigXCoor$y)
				all_ranges <- rbind(all_ranges,c(IDcase,it,nNeig,xrange,yrange))
			  }

				#################### Physical distance
			  all_cell_nb_in_bandwidth <- tryCatch(dnearneigh(as.matrix(coords), 0, bandwidth, longlat = NULL),error=function(e){cat("ERROR :",conditionMessage(e), "\n");skip_to_next <<- TRUE})

			  int_nb_list <- list()
				for(i in 1:length(nb_list)){
				  int_nb <- intersect(nb_list[[i]],all_cell_nb_in_bandwidth[[i]])
					int_nb_list <- append(int_nb_list,list(int_nb))
				}

			  ###################### Building Empty CLQ Matrix
			  cell_type_num <- unique(cell_type_assignment)
			  cell_types <- unique(cell_type_assignment)
			  CLQ_matrix <- matrix(0L, nrow = length(cell_type_num), ncol = length(cell_type_num))
			  row.names(CLQ_matrix) <- cell_types
			  colnames(CLQ_matrix) <- cell_types

			  ###################### Estimating CLQ
			  for(k in 1:length(cell_type_num)){
				for(j in 1:length(cell_type_num)){
				  cell_A <- cell_type_num[k]
				  cell_B <- cell_type_num[j]

				  cell_B_name <- cell_types[j]
				  CLQ_result <- get_CLQ(cell_type_assignment,coords,int_nb_list,cell_A,cell_B,
										number_of_neighbor=nNeig,cell_B_name)
				  CLQ_matrix[k,j] <- CLQ_result
				}
			  }
			  CLQsub <- melt(CLQ_matrix) #Var1 are row names (cell A) and Var2 column names (cell B). CLQ measures how much the cell A is attracted to B (B=target cell type, A=neighboring cells types) or in other words degree of cell type B co-locales with cell type A : CLQ_ab <- (Cab/length(cell_A_indices))/(length(cell_B_indices)/(total_cell-1)). CLQAB = 2 means that A is twice as likely to have B as its nearest neighbor as would be expected by chance
			  CLQsub <- data.frame(Sample=IDcase,Iteration=it,CLQsub)
			  CLQsub_all <- rbind(CLQsub_all,CLQsub)

			  close(pb)
		}


	write.table(CLQsub_all,file.path(dir_output,paste0(IDcase,"_CLQ_resampled.txt")),sep="\t",row.names=TRUE,col.names=NA)

	colnames(tiles_order) <- c("Sample",idstiles)
	write.table(tiles_order,file.path(dir_output,paste0(IDcase,"tiles_order.txt")),sep="\t",row.names=TRUE,col.names=NA)

	# Overview of each tile: range of coordinates and number of cells included in each tile
	colnames(ncells_all) <- c("Sample",
				"Tile x coords start",
				"Tile x coords end",
				"Tile y coords start",
				"Tile y coords end",
				"Number of cells in the tile",
				"Total number of cells original image"
				)
	ncells_all <- as.data.frame(ncells_all)
	write.table(ncells_all,file.path(dir_output,paste0(IDcase,"_OverviewCellTiles.txt")),sep="\t",col.names=TRUE)

	# Same overview but in the form of a figure
	for(subj in unique(ncells_all$Sample)){
		ncells_all_subj <- ncells_all[ncells_all$Sample %in% subj,]
		a <- ggplot(as.data.frame(ncells_all_subj), aes(x=ncells_all_subj[,2], y=ncells_all_subj[,4])) +
		  #geom_point() +
		  geom_text(label=ncells_all_subj[,6]) +
		  theme_void() +
		  labs(title="N cells per tile",x ="x", y = "y")
		ggsave(file.path(dir_output,paste0("OverviewCellTiles_",subj,".png")),width=4,height=4)
	}

	############################################## Calculating CLQ from original image


	  ## Select a subject and select columns of interest
	  cell_type_assignment <- data_codexSub$cell_type
	  coords <- data.frame(x=data_codexSub$x,y=data_codexSub$y)
	  row.names(coords) <- data_codexSub[,1]

	  ####################### Identifying Neighbors
	  nb_list <- KNN_neighbors(as.matrix(coords),number_of_neighbors=nNeig)

		  for(i in length(nb_list)){
			neigX <- nb_list[[i]]
			neigXCoor <- coords[neigX,]
			xrange <- max(neigXCoor$x)-min(neigXCoor$x)
			yrange <- max(neigXCoor$y)-min(neigXCoor$y)
			all_ranges <- rbind(all_ranges,c(IDcase,0,nNeig,xrange,yrange))
		  }

		#################### Physical distance
	  all_cell_nb_in_bandwidth <- tryCatch(dnearneigh(as.matrix(coords), 0, bandwidth, longlat = NULL),error=function(e){cat("ERROR :",conditionMessage(e), "\n");skip_to_next <<- TRUE})
	  #print(paste("Average number cells when imposing physical distance limit:", mean(lengths(all_cell_nb_in_bandwidth))))

	  int_nb_list <- list()
		for(i in 1:length(nb_list)){
		  int_nb <- intersect(nb_list[[i]],all_cell_nb_in_bandwidth[[i]])
			int_nb_list <- append(int_nb_list,list(int_nb))
		}

	  #print(paste("Average number of neighbor cells after imposing physical distance limit:", mean(lengths(int_nb_list))))

	  NcellsNeig <- data.frame(CaseID = IDcase,nTile = nTile,nNeig=nNeig,bandwidth=bandwidth,
	  		TotalCells=length(all_cell_nb_in_bandwidth),Av.Physical=mean(lengths(all_cell_nb_in_bandwidth)),Av.NeighPhysical=mean(lengths(int_nb_list)))
		write.table(NcellsNeig,file.path(dir_output,paste0(IDcase,"_NcellsNeig.txt")),sep="\t")

	  ######################## Building Empty CLQ Matrix
	  cell_type_num <- unique(cell_type_assignment)
	  cell_types <- unique(cell_type_assignment)
	  CLQ_matrix <- matrix(0L, nrow = length(cell_type_num), ncol = length(cell_type_num))
	  row.names(CLQ_matrix) <- cell_types
	  colnames(CLQ_matrix) <- cell_types

	  ######################## Estimating CLQ
	  for(k in 1:length(cell_type_num)){
	    for(j in 1:length(cell_type_num)){
	      cell_A <- cell_type_num[k]
	      cell_B <- cell_type_num[j]

	      cell_B_name <- cell_types[j]
	      CLQ_result <- get_CLQ(cell_type_assignment,coords,int_nb_list,cell_A,cell_B,
	                            number_of_neighbor=nNeig,cell_B_name)
	      CLQ_matrix[k,j] <- CLQ_result
	      #CLQmax_matrix[k,j] <- CLQ_result[[1]]
	    }
	  }
	  filename <- paste0(IDcase,"_CLQ.txt")
	  write.table(CLQ_matrix,file.path(dir_output,filename),sep="\t",row.names=TRUE,col.names=NA)
	  CLQsub <- melt(CLQ_matrix)
	  CLQsub <- data.frame(Sample=IDcase,Iteration=0,CLQsub)
	  CLQsub_all <- rbind(CLQsub,CLQsub_all)

	write.table(CLQsub_all,file.path(dir_output,paste0(IDcase,"_CLQ_realANDresampled.txt")),sep="\t",row.names=TRUE,col.names=NA)


	############################### Range of coordinates for each iteraction across Nneig
	colnames(all_ranges) <- c("Sample","Iteration","Nneigh","Xrange","Yrange")
	write.table(all_ranges,file.path(dir_output,paste0(IDcase,"_RangeCoords.txt")),sep="\t",col.names=TRUE)

	############################### Difference between mean bootstrap and real CLQ

		CLQsub_all$Pair <- paste(CLQsub_all[,3],CLQsub_all[,4])
		CLQsub_all <- CLQsub_all[,c(2,6,5)]
		CLQ_subj_CLQ <- CLQsub_all[which(CLQsub_all$Iteration == 0),]
		CLQ_subj_boot <- CLQsub_all[which(CLQsub_all$Iteration != 0),]
		CLQ_subj_boot <- CLQ_subj_boot %>%
		  group_by(Pair) %>%
		  summarise_at(vars(value), list(name = mean))
		CLQ_subj_boot <- as.data.frame(CLQ_subj_boot)
		colnames(CLQ_subj_boot) <- c("Pair","value_boot")

		CLQ_subj_boot <- CLQ_subj_boot[match(CLQ_subj_CLQ$Pair,CLQ_subj_boot$Pair),]
		CLQ_subj_CLQ <- merge(CLQ_subj_CLQ,CLQ_subj_boot,by="Pair")
		#CLQ_subj_CLQ$value[CLQ_subj_CLQ$value == 0] <- 0.000000001
		perc <- CLQ_subj_CLQ$value_boot/CLQ_subj_CLQ$value
		perc <- data.frame(Sample=subj,Pair=CLQ_subj_CLQ$Pair,perc=perc)
	write.table(perc,file.path(dir_output,paste0(IDcase,"_Difference_realANDresampled.txt")),sep="\t",col.names=TRUE)

	p<-ggplot(data=perc, aes(x=Pair, y = perc)) +
		 geom_point() +
		 theme_void() +
		 theme(axis.text.x = element_text(size=7,angle = 90, vjust = 0.5, hjust=0.5))
	ggsave(file.path(dir_output,paste0(IDcase,"_Difference_realANDresampled.png")),p,height=5,width=30)


	############################### Difference number of cells between real image and bootstrapped images

	  Ncells_2 <- data.frame(CaseID = rep(IDcase,length(Ncells)),nTile = rep(nTile,length(Ncells)),nNeig=rep(nNeig,length(Ncells)),bandwidth=rep(bandwidth,length(Ncells)),
	  		Image=c("Real",rep("Bootstrap",(length(Ncells)-1))),TotalCells=Ncells)
		write.table(Ncells_2,file.path(dir_output,paste0(IDcase,"_TotalNcells.txt")),sep="\t",col.names=TRUE,row.names=FALSE)

  return(list(BootstrapCLQs=CLQsub_all, CLQ = CLQ_matrix))

}


#module load R/4.0.2
#module load gcc/10.1.0
#module load physics udunits gdal
#module load proj/4.9.3
#module load geos

######################################################## FUNCTIONS
### Get CLQ
get_CLQ <- function(cell_type_assignment,coords,nb_list,cell_A,cell_B,
                    number_of_neighbor,cell_B_name){
  ### Cell B is the cell type to check, cell A is attracted to B
  ### CLQ in this case is the degree to which A can be attached to B
  cell_A_indices <- which(cell_type_assignment == cell_A)
  cell_B_indices <- which(cell_type_assignment == cell_B)

  cell_A_nb_of_B <- find_cell_type_neighbors(cell_B,nb_list,cell_type_assignment,
                                             cell_A_indices,cell_B_name)
  total_cell <- length(which(cell_type_assignment!=0 & is.na(cell_type_assignment)==FALSE))
  Cab <- 0
  for(i in 1:length(cell_A_indices)){
    Cab <- Cab + length(cell_A_nb_of_B[i,1][[1]])
  }
  if(cell_A==cell_B){
    if(length(cell_B_indices)<=5 | length(cell_A_indices)<=5){
      CLQ_ab <- 0
    }else{
      CLQ_ab <- (Cab/length(cell_A_indices))/((length(cell_B_indices)-1)/(total_cell-1))
    }
  }else{
    if(length(cell_B_indices)<=5 | length(cell_A_indices)<=5){
      CLQ_ab <- 0
    }else{
      CLQ_ab <- (Cab/length(cell_A_indices))/(length(cell_B_indices)/(total_cell-1))
    }
  }
  return(CLQ_ab)
}
############find_cell_type_neighbors
find_cell_type_neighbors <- function(cell_B,nb_list,cell_type_assignment,
                                     cell_A_indices,cell_B_name){
  if(length(cell_A_indices)==0){

  }else{
    cell_type_nb <- matrix(rep(list(),length(cell_A_indices)),nrow=length(cell_A_indices),ncol=1)
    row.names(cell_type_nb) <- cell_A_indices
    colnames(cell_type_nb) <- cell_B_name
    for(j in 1:length(cell_A_indices)){
      #print(j)
      current_cell_ID <- cell_A_indices[j]
      neighbors <- nb_list[[current_cell_ID]]
      neighbor_types <- cell_type_assignment[neighbors]

      cell_type_neighbor_IDs <- neighbors[which(neighbor_types == cell_B)]
      cell_type_nb[j,1][[1]] <- cell_type_neighbor_IDs
    }
    return(cell_type_nb)
  }
}
############### KNN_neighbors
KNN_neighbors <- function(coords,number_of_neighbors){
  xxx <- knearneigh(coords,k=number_of_neighbors)
  nb_list <- list()
  for(i in 1:dim(xxx$nn)[1]){
    nb_list[[i]] <- xxx$nn[i,]
  }
  return(nb_list)
}


##############
GetNeighborInfo <- function(coords, number_of_neighbors = 5, bandwidth = 100) {
  print("Getting the nearest neighbors")
  nb_list <- knearneigh(coords, k = number_of_neighbors)$nn
  colnames(nb_list) <- paste0("neighbor", seq(1, number_of_neighbors, by = 1))

  print("Identifying neighboring cells within a defined circle bandwidth")
  all_cell_nb_in_bandwidth <- dnearneigh(coords, 0, bandwidth, longlat = NULL)

  print("Identify distances for all the cells within the circle bandwidth")
  cell_nb_dist <- nbdists(all_cell_nb_in_bandwidth, coords)
  return(list(nb_list, all_cell_nb_in_bandwidth, cell_nb_dist))
}


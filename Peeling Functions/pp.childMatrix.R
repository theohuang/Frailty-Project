#' Creates a Binary 3D Array Which Tells You Whose Parents are Who
#' 
#' Entry \code{(ijk)} = 1 if individual \code{k} is the offspring of individuals 
#' code{i} and code{j}
#' @param ped the pedigree
#' @details Author: Tom Madsen
#' 
#' Pre-processing function
#' 
#' This only includes individuals whose parents are both in the pedigree
#' @family peeling
pp.childMatrix=function(ped){
  # Initialze 3d array of zeros
  C = array(0, rep(nrow(ped), 3))
  # Iterate through all the individuals in `ped`
  for(i in 1:nrow(ped)){
    # If both parents of individual i are in the pedigree, assign them as the 
    # parents
    if(ped$FatherID[i] != 0 & ped$MotherID[i] != 0){
      C[ped$ID==ped$MotherID[i], ped$ID==ped$FatherID[i], i] = 1
      C[ped$ID==ped$FatherID[i], ped$ID==ped$MotherID[i], i] = 1
      # If only one parent is in the pedigree, assign them to both parents
    } else if(ped$FatherID[i] != 0) {
      C[ped$ID==ped$FatherID[i], ped$ID==ped$FatherID[i], i] = 1
    } else if (ped$MotherID[i] != 0) {
      C[ped$ID==ped$MotherID[i], ped$ID==ped$MotherID[i], i] = 1
    }
  }
  return(C)
}
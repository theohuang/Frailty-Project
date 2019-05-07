#' Creates a Binary Matrix Which Tells You Who's Mated with Whom
#' 
#' Entry \code{(i,j)} = 1 if individual \code{i} mated with individual \code{j}
#' @param ped the pedigree
#' @details Author: Tom Madsen
#' 
#' Pre-processing function
#' @family peeling
pp.edgeMatrix=function(ped){
  # Initialze matrix of zeros
  E = matrix(0, nrow=nrow(ped), ncol=nrow(ped))
  # Iterate through all the individuals in `ped`
  for(i in 1:nrow(ped)){
    # If both parents of individual i are in the pedigree, assign them as each 
    # other's mates. 
    if(ped$FatherID[i] != 0 & ped$MotherID[i] != 0){
      E[ped$ID==ped$FatherID[i], ped$ID==ped$MotherID[i]] = 1
      E[ped$ID==ped$MotherID[i], ped$ID==ped$FatherID[i]] = 1
    # If only one parent is in the pedigree, assign them as their own mate
    } else if (ped$FatherID[i] != 0){
      E[ped$ID==ped$FatherID[i], ped$ID==ped$FatherID[i]] = 1
    } else if (ped$MotherID[i] != 0) {
      E[ped$ID==ped$MotherID[i], ped$ID==ped$MotherID[i]] = 1
    }
  }
  return(E)
}
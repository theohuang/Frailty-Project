#' Calculates the Posterior Probability for a Particular Individual and 
#' Potential Genotype
#' 
#' The posterior probability is the probability of the genotype and all posterior  
#' information. Posterior information is information garnered via mates and 
#' offspring. See Fernando et al. for details
#' @param index the index (row) number of the individual who anterior probability 
#' is being computed
#' @param mate the index of the individual's spouse
#' @param ped the pedigree
#' @param E edge matrix (binary matrix which tells you who's mated with whom)
#' @param C child matrix (binary 3D array which tells you whose parents are who)
#' @param prev a vector of genotype prevalences (see \code{pp.prevalences})
#' @param LIK a likelihood matrix
#' @param TR a 3D array which gives the probability of a genotype given a set of 
#' parental genotypes (see \code{pp.inheritanceProb})
#' @details Author: Tom Madsen
#' 
#' Helper function called by \code{pp.peeling_paring}
#' @family peeling
pp.sub = function(index, mate, ped, E, C, prev, LIK, TR){
  if(!any(is.na(P[ped$ID==index,ped$ID==mate,]))){
    return(P[ped$ID==index,ped$ID==mate,])
  }
  else{
    U = length(prev)
    
    B = matrix(1, nrow=U, ncol=U)
    children = ped$ID[C[ped$ID==index,ped$ID==mate,]==1]
    for(child in children){
      spouse = ped$ID[E[ped$ID==child,]==1]
      if(length(spouse) == 0){spouse = child}
      S = sweep(TR, 3, LIK[ped$ID==child, ] *
                  apply(sapply(spouse, pp.sub, 
                               index=child, ped=ped, E=E, C=C, prev=prev, LIK=LIK, TR=TR), 
                        1, prod), '*')
      B = B*apply(S, c(1,2), sum)
    }
    
    if (index != mate && length(children) > 0) { ### if child has only 1 parent in pedigree
      B = sweep(B, 1, LIK[ped$ID==mate,]*pp.ant(mate, ped, E, C, prev, LIK, TR),'*')
    } else { 
      B = sweep(B, 1, prev, '*')
    }
    
    
    # REVIEW THIS -- PRETTY SURE IT'S NOT IN THE RIGHT PLACE
    other_spouses = ped$ID[E[ped$ID==mate,]==1]
    other_spouses = other_spouses[!(other_spouses %in% c(index,mate))]
    if (length(other_spouses) > 0){
      for(spouse in other_spouses){ 
        B = sweep(B, 2, pp.sub(mate, spouse, ped, E, C, prev, LIK, TR),'*')
      }
    }
    
    P #here because of stupid buggy global variables
    P[ped$ID==index,ped$ID==mate,] <<- apply(B, 2, sum) # apply(B, 3, sum)
    return(P[ped$ID==index,ped$ID==mate,])
  }
}
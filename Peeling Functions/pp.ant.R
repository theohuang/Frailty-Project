#' Calculates the Anterior Probability for a Particular Individual and 
#' Potential Genotype
#' 
#' The anterior probability is the probability of the genotype and all anterior 
#' information. Anterior information is information garnered via blood relatives 
#' (rather than mates). See Fernando et al. for details. 
#' @param index the index (row) number of the individual who anterior probability 
#' is being computed
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
pp.ant = function(index, ped, E, C, prev, LIK, TR){
  if (!any(is.na(A[ped$ID==index,]))){
    return(A[ped$ID==index,])
  } 
  else{
    U = length(prev)
    m = ped$MotherID[ped$ID==index]
    f = ped$FatherID[ped$ID==index]
    
    #see Lorenzo's notes on the algorithm to understand this part - just using the formula for anterior probability
    B = matrix(1, nrow=U, ncol=U)
    sibs = ped$ID[C[ped$ID==m,ped$ID==f,]==1]
    sibs = sibs[sibs!=index]
    if(length(sibs) > 0){
      for(sib in sibs){
        spouse = ped$ID[E[ped$ID==sib,]==1]
        if(length(spouse) == 0){
          spouse = sib
        }
        S = sweep(TR, 3, LIK[ped$ID==sib, ] *
                    apply(sapply(spouse, pp.sub,
                                 index=sib, ped=ped, E=E, C=C, prev=prev, LIK=LIK, TR=TR),
                          1, prod), '*')
        B = B*apply(S, c(1,2), sum)
      }
    }
    
    if (m != 0){ ### Skip this step if there is no mother, to avoid NAs
      B = sweep(B, 1, LIK[ped$ID==m,]*pp.ant(m, ped, E, C, prev, LIK, TR), '*')
    } else if (f != 0) { 
      B = sweep(B, 1, prev, '*')
    }
    
    # Pull out all of the mother's spouses except for herself and the father
    mom_spouses = ped$ID[E[ped$ID==m,]==1]
    mom_spouses = mom_spouses[!(mom_spouses %in% c(m,f))]
    if (length(mom_spouses) > 0){
      for(spouse in mom_spouses){
        B = sweep(B, 1, pp.sub(m, spouse, ped, E, C, prev, LIK, TR),'*')
      }
    }
    
    if (f != 0) { ### Skip this step if there is no father, to avoid NAs
      B = sweep(B, 2, LIK[ped$ID==f,]*pp.ant(f, ped, E, C, prev, LIK, TR),'*')
    } else if (m != 0) { 
      B = sweep(B, 2, prev, '*')
    }
    
    # Pull out all of the father's spouses except for himself and the mother
    dad_spouses = ped$ID[E[ped$ID==f,]==1]
    dad_spouses = dad_spouses[!(dad_spouses %in% c(m,f))]
    if (length(dad_spouses) > 0){
      for(spouse in dad_spouses){
        B = sweep(B, 2, pp.sub(f, spouse, ped, E, C, prev, LIK, TR),'*')
      }
    }
    
    B = sweep(TR, c(1,2), B, '*')
    
    A #here because of stupid buggy global variables
    A[ped$ID==index,] <<- apply(B, 3, sum) 
    return(A[ped$ID==index,])
  }
}
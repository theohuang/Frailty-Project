#' Get the Allele Pairs for Each Gene with Mutation Restriction
#' 
#' Returns a matrix with all permissible genotypes under the mutation restriction
#' @param K total number of genes
#' @param T max number of mutations allowed
#' @details Author: Gang Peng
#' 
#' Adapted by: Tom Madsen
#' 
#' Pre-processing function
#' @family peeling
pp.possibleGenotypes = function(K, T){
  
  if(length(K) > 1){
    K = length(K)
  }
  
  ####Gene number for each gene
  G <- array(0, c(2,K))
  G[,1]=c(0,1)
  
  if(K > 1){
    for(i in 2:K){
      # all the genotypes with T-1 or fewer mutations before locus i
      GG=as.matrix(G[(G%*%rep(1,K))<T,]) 
      # if that's only one genotype, transpose it
      if((dim(GG)[2])==1){GG=t(GG)} 
      # for these genotypes, it's possible to have a mutation at locus i under the restriction
      GG[,i]=1 
      
      G=rbind(G,GG)
      # there's an implication here that you can have at most two mutations at each locus
      # each allele is either WT(=0) or mutant(>0); the numbering of mutant types is non-informative
    }
  }
  PG = c()
  total_muts = apply(G, 1, sum)
  for(i in 0:T){
    PG = rbind(PG, G[total_muts==i,])
  }
  return(PG)
}
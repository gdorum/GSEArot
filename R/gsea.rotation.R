#' GSEA rotation
#' 
#' GSEA with rotation test
#' @param X Gene expression matrix with samples as rows and genes as columns.
#' @param y Optional vector of 0/1 indicating sample phenotype.
#' @param S Matrix indicating gene set membership with genes as rows and gene sets as columns. 1 indicates that the gene is a member, 0 indicates not a member.
#' @param nrot The number of rotations to perform. Default is 500.
#' @param fun Name of function for computing a gene wise test statistic.
#' @return  \item{p.values }{A p-value for each gene set}
#' \item{q.values }{An FDR q-value for each gene set}
#' @references 
#' Dorum, G., Snipen, L., Solheim, M. and Sabo (2009) Rotation Testing in Gene Set Enrichment
#' Analysis for Small Direct Comparison Experiments. \emph{Statistical Applications in Genetics
#'  and Molecular Biology}, \bold{8}(1), article 34.
#'  
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#' Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.
#' and Mesirov, J. P (2005) Gene set enrichment analysis: A knowledge-based
#' approach for interpreting genome-wide expression profiles, \emph{PNAS}, \bold{102},
#' 15545-15550.
#' @author Guro Dorum
#' @export
gsea.rotation <-
function(X,y=NULL,S,nrot=500,fun=t90th) {
    
    nset <- ncol(S)
    
    stat <- fun(X,y)       
        
    #Calculate an enrichment score for each gene set
    escore <- numeric(nset)
    for( i in 1:nset ) {
        escore[i] <- es(stat,S[,i])
    }
    
    #Calculate a null distribution for each gene set
    es.null <- matrix(nrow=nrot,ncol=nset)
    for(j in 1:nrot) {
    
        if( j%%100 == 0 ) cat(paste(j,"rotations completed... \n"))
        
       #Rotation method depends on the data type
        if( is.null(y) ) {
            Xrot <- rotation(X, method=1) 
        } else Xrot <- rotation(X, method=2)
                
        #Compute test statistics and enrichment scores for rotated data
        testRot <- fun(Xrot,y) 
        for(i in 1:nset) {   
            es.null[j,i] <- es(testRot,S[,i])
        }
    }
    
    #Estimating p-values and normalised enrichment scores
    pvals <- numeric(nset)
    NES <- numeric(nset)
    NES.null <- matrix(nrow=nrot,ncol=nset)
    for( i in 1:nset) {
        res <- significance(escore[i], es.null[,i])
        pvals[i] <- res$p.value            
        NES[i] <- res$NES                  
        NES.null[,i] <- res$NESnull        
    }
    
    #Estimating FDR q-values
    qvals <- sapply(NES, FUN=fdr, NESobs=NES, NESnull=NES.null)
    
    return(list(p.values=pvals,q.values=qvals))
}

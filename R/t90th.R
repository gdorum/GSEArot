#' Moderated gene wise t-statistics.
#'
#' Computes gene wise t-statistic with variance stabilised by adding the 90th percentile.
#' 
#' @param X Gene expression matrix with samples as rows and genes as columns
#' @param y Optional vector of 0/1 indicating sample phenotype
#' @return Vector with one smoothed t-value for each gene
#' @details Computes a gene wise t-statistic with variance stabilised by adding the 90th percentile 
#' of all the genes' standard deviations to a gene's standard devitation (Efron et al. 2001).
#' If y is given, a two-sample test is performed, otherwise a one-sample test is performed.
#' For a two-sample t-test the approach of Pan et al. (2003) is used.
#' @references Efron, B., Tibshirani, R., Storey, J. D. and Tusher, V. (2001) 
#' Empirical Bayes analysis of a microarray experiment, 
#' \emph{J Amer Statist Assoc}, 
#' \bold{96}, 1151-1160. 
#' Pan,W., Lin,J. and Le,C.T. (2003)
#' A mixture model approach to detecting differentially expressed genes with microarray data, 
#' \emph{Funct Integr Genomics}, 
#' \bold{3}, 117-124.
#' @author Guro Dorum
#' @importFrom stats quantile sd var
#' @export
t90th <- function( X, y=NULL ) {
    
    if( is.null(y) ){
        #One-sample test
        n <- ncol(X)
        m <- apply(X,2,FUN=mean,na.rm=TRUE)
        std <- apply(X,2,FUN=sd,na.rm=TRUE)
        percStd <- quantile(std,probs=0.9)
        stat <- m / ( (std + percStd) / sqrt(2*n) )
    } else {   
        #Two-sample test 
        n1 <- sum(y==0)
        n2 <- sum(y)
        m1 <- apply(X[y==0,],2,FUN=mean,na.rm=TRUE)
        m2 <- apply(X[y==1,],2,FUN=mean,na.rm=TRUE)
        v1 <- apply(X[y==0,],2,FUN=var,na.rm=TRUE)
        v2 <- apply(X[y==1,],2,FUN=var,na.rm=TRUE)
        std <- sqrt(v1/n1 + v2/n2)
        percStd <- quantile(std,probs=0.9)
        stat <- (m1 - m2) / (std + percStd)
    }
        
    stat
}

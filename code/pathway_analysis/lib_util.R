


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Numeric Helper Functions
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#' Compute the sum of values grouoed into integer bins
#' Similar to tabulate(), but compute the sum of values instead of counting
#' @param bin integer vector defining the bins
#' @param nbins number of bin (length of the output vector)
#' @param w the weight vector (values to sum)
#' @return a numeric vector of size nbins which sum the values of w into each bin
tabulate_sum <- function(bin,nbins=max(1,bin,na.rm=TRUE),w=rep_len(1L,length(bin))) {
  as.vector(Matrix::sparseMatrix(i=bin,j=rep_len(1L,length(bin)),x=w,dims = c(nbins,1L)))
}


#' Compute columns means accross rows 
#' Similar to rowsum, but compute the mean instead
#' @param M a (n x m) matrix
#' @param f a factor grouping matrix rows
#' @param ... additional values are passed to rowsum
#' @return a numeric matrix of means
rowmean <- function(x,group,...) {
  group <- as.factor(group)
  m <- rowsum(x,group,...)
  n <- rowsum(0+!is.na(x),group,...)
  x <- m / n
  x[levels(group),]
}


#' Compute margin difference between two conditions 
#' @param A numeric matrix
#' @param B numeric matrix of the reference values
#' @param q tolerance
#' @return a numeric vector of differences
rowMargin <- function(A,B,q=1L) {
  pmax(Biobase::rowQ(A,q) + Biobase::rowQ(-B,q),0) - pmax(Biobase::rowQ(B,q) + Biobase::rowQ(-A,q),0)
}


clip <- function(m,low=-1,high=1) {
  stopifnot(low<=high)
  m[m<low] <- low
  m[m>high] <- high
  m
}




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Strings Processing Functions
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#' Trim trailing white spaces in a string 
str_chomp <- function(txt) {sub("\\s*$","", txt)}

#' Insert new lines into long strings
str_break <- function(txt,len=20) {
  pat <- sprintf('(.{1,%d})(\\s|$)',len)
  str_chomp(gsub(pat, '\\1\n', txt))
}







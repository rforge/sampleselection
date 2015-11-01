
### intervals: extract the intervals by observation from the model
### return: a matrix with columns ub, lb, the corresponding interval boundaries
intervals <- function(object) 
   UseMethod("intervals")

intervals.intReg <- function(object) {
   r <- model.response(model.frame(object))
   if(is.null(dim(r))) {
                           # return a vector -> factor if interval codes
      b <- boundaries(object)
      names(b) <- NULL
      LB <- b[as.integer(r)]
      UB <- b[as.integer(r) + 1]
      return(cbind(LB=LB, UB=UB))
   }
   ## observation-specific boundaries -> already in the right form
   colnames(r) <- c("LB", "UB")
                           # may have different names in the user-supplied
                           # original data
   return(r)
}

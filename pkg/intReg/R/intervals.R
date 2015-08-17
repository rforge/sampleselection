
### intervals: extract the intervals by observation from the model
### return: a matrix with columns ub, lb, the corresponding interval boundaries
intervals <- function(object) 
   UseMethod("intervals")

intervals.intReg <- function(object) {
   r <- model.response(model.frame(object))
   if(is.null(dim(r))) {
                           # return a vector -> factor if interval codes
      lb <- boundaries(object)[as.integer(r)]
      ub <- boundaries(object)[as.integer(r) + 1]
      return(cbind(lb, ub))
   }
   ## observation-specific boundaries -> already in the right form
   return(r)
}

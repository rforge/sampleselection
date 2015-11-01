
### intervalObs: return TRUE/FALSE depending on if the observation was treated
### as an interval or point observation

intervalObs <- function(object)
   UseMethod("intervalObs")

intervalObs.intReg <- function(object) {
   object$param$intervalObs
}


### Boundaries: return the set of interval boundaries used in the model
### Note that these may be overlap in case of observation-specific boundaries
intervalObs <- function(object)
   UseMethod("intervalObs")

intervalObs.intReg <- function(object) {
   object$param$intervalObs
}

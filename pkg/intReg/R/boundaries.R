
### Boundaries: return the set of interval boundaries used in the model
### Note that these may be overlap in case of observation-specific boundaries
boundaries <- function(object)
   UseMethod("boundaries")

boundaries.intReg <- function(object) {
   object$param$boundaries
}


### disturbancies: return the disturbancies type of the interval regression model

disturbances <- function(object)
   UseMethod("disturbances")

disturbances.intReg <- function(object) {
   object$method
}

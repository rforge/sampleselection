
### disturbancies: return the disturbancies type of the interval regression model

disturbancies <- function(object)
   UseMethod("disturbancies")

disturbancies.intReg <- function(object) {
   object$method
}


### minIntervalWidth: fetch the minimal interval width data from intReg object

setGeneric("minIntervalWidth",
           function(object, ...) standardGeneric("minIntervalWidth")
           )

setMethod("minIntervalWidth", "intReg",
          function(object, ...) object$param$minIntervalWidth)


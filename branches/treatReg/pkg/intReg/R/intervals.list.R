### fetch interval information from a list or data.frame

intervals.list <- function(object) {
   if(!is.null(object$yInt)) {
      if(!is.factor(object$yInt)) {
         cat("Component 'yInt' must be a factor with levels ",
             "[LB, UB] or (LB, UB] or [LB, UB) or (LB, UB)",
             fill=TRUE)
         stop("'yInt' is not a factor")
      }
      boundaries <- strsplit(levels(object$yInt), ",")
      if(!all(l <- sapply(boundaries, length) == 2)) {
         cat("Levels of 'yInt' must be in the form ",
             "[LB, UB] or (LB, UB] or [LB, UB) or (LB, UB)",
             fill=TRUE)
         cat("Levels that do not match:\n")
         print(levels(object$yInt)[l != 2])
         stop("wrong levels in 'yInt'")
      }
      LBv <- sapply(boundaries, "[", 1)
      UBv <- sapply(boundaries, "[", 2)
      LBv <- as.numeric(sub("^[[(]", "", LBv))
      UBv <- as.numeric(sub("[])]$", "", UBv))
      LB <- LBv[object$yInt]
      UB <- UBv[object$yInt]
   }
   else if(!is.null(object$LB) & !is.null(object$UB)) {
      if(length(object$LB) != length(object$UB)) {
         stop("'LB' and 'UB' must be of the same length for intervals")
      }
      LB <- object$LB
      UB <- object$UB
   }
   else {
      cat("Intervals can be submitted either as\n",
          "* a factor named 'yInt' with levels in the form\n",
          "  [lb, ub] or (lb, ub] or [lb, ub) or (lb, ub)\n",
          "  (this is the first choice)\n",
          "* or as variables called 'LB' and 'UB'\n",
          "  if 'yInt' not found\n")
      stop("No suitable intervals found in data")
   }
   return(cbind(LB=LB, UB=UB))
}

intervals.data.frame <- function(object)
   intervals.list(object)

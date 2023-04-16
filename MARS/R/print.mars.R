print.mars <- function(x,...){
  cat("\nCall:\n")
  print(x$call)
  cat("Coefficients:\n")
  print(coefficients(x))
  invisible(x)
}


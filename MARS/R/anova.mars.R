anova.mars <- function(mars_output) {
  # anova.lm(mars_output)
  anova(mars_output$formula)
}

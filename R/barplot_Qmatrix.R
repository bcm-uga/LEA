barplot.Qmatrix = function(height, sort.by.Q = TRUE, lab = FALSE, ...){
  Q = height
  rm(height)

  if (class(Q)[1] != "Qmatrix") {warning("Object Q not of class Qmatrix.")}

  if (sort.by.Q) {
    gr = apply(Q, MARGIN = 1, which.max)
    gm = max(gr)
    gr.o = order(sapply(1:gm, FUN = function(g) mean(Q[,g])))
    gr = sapply(gr, FUN = function(i) gr.o[i])
    or = order(gr)
    Qm = t(Q[or,])
    class(Qm) = "matrix"
    graphics::barplot(Qm,...)
    return(list(order = or))
    }
  else {
    Qm = t(Q)
    class(Qm) = "matrix"
    graphics::barplot(Qm,...)
    return(list(order = 1:nrow(Q)))
      }
}
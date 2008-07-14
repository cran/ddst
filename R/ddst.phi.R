`ddst.phi` <-
function(x, j, base = ddst.base.legendre) {
  .Internal(mean(base(x,j)))
}


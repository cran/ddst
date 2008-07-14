`ddst.exp.test` <-
function(x, base = ddst.base.legendre, c = 100, B=1000, compute.p = F, Dmax = 5, ...) {
 data(ddst.polynomial.fun)
 data(MMexp)

# method.name = as.character(substitute(base)) 
 # only Legendre is implemented yet
 base = ddst.base.legendre
 method.name = "ddst.base.legendre"

n = length(x)
if(n<5) 
	stop("length(x) should be at least 5")

er = .Internal(mean(x))
maxN = max(min(Dmax, n-2, 20),1)
coord = numeric(maxN)
u = numeric(maxN)
for (j in 1:maxN) { 
u[j] = ddst.phi(pexp(x,1/er), j, base)
coord[j] = t(u[1:j]) %*% MMexp[[j]] %*% u[1:j] * n
}

l = ddst.IIC(coord, n, c)
attr(l, "names") = "n. coord"
t = coord[l]
attr(t, "names") = "WT*"
result = list(statistic = t, parameter = l, method = "Data Driven Smooth Test for Expotentiality")
result$data.name = paste(paste(as.character(substitute(x)), collapse=""), ",   base: ", method.name, ",   c: ", c, sep="")
class(result) = "htest"
if (compute.p) {
tmp = numeric(B)
for (i in 1:B) {
y = rexp(length(x))
tmpC = ddst.exp.Nk(y, base, Dmax = Dmax, n=length(y))
l = ddst.IIC(tmpC, n, c)
tmp[i] = tmpC[l]
}
p.val = .Internal(mean(tmp > t))
result$p.value = p.val  
}
result
}


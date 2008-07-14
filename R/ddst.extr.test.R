`ddst.extr.test` <-
function(x, base = ddst.base.legendre, c = 100, B=1000, compute.p = F, Dmax = 5, ...) {
 data(ddst.polynomial.fun)
 data(MMextr12)
 data(MMextr)
 data(sM22)
 library(evd)

# method.name = as.character(substitute(base)) 
 # only Legendre is implemented yet
 base = ddst.base.legendre
 method.name = "ddst.base.legendre"

n   = length(x)
if(n<5) 
	stop("length(x) should be at least 5")
sx  = sort(x)
er2 = sum(sx*(1:n - n:1))/(n*(n-1)*log(2))
er1 = .Internal(mean(x))+0.5772156649*er2
maxN = max(min(Dmax, length(x)-2, 20),1)

u = NULL
for (j in 1:maxN) 
u[j] = ddst.phi(1-pgumbel(-x,-er1,er2), j, base)
coord = NULL
gg1=.Internal(mean(1-exp((x-er1)/er2)))
gg2=.Internal(mean((1-exp((x-er1)/er2))*(x-er1)/er2+1))
for (k in 1:Dmax) {
korekta = u[1:k] + t(MMextr12[[k]]) %*% sM22 %*% c(gg1,gg2)
    coord[k] = t(korekta) %*% MMextr[[k]] %*% korekta * n
}

l = ddst.IIC(coord, n, c)
attr(l, "names") = "n. coord"
t = coord[l]
attr(t, "names") = "WT*"
result = list(statistic = t, parameter = l, method = "Data Driven Smooth Test for Extreme Values")
result$data.name = paste(paste(as.character(substitute(x)), collapse=""), ",   base: ", method.name, ",   c: ", c, sep="")
class(result) = "htest"
if (compute.p) {
tmp = numeric(B)
for (i in 1:B) {
y = rnorm(n)
tmpC = ddst.extr.Nk(y, base, Dmax = Dmax, n=length(y))
l = ddst.IIC(tmpC, n, c)
tmp[i] = tmpC[l]
}
p.val = .Internal(mean(tmp > t))
result$p.value = p.val  
}
result
}


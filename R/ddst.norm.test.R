`ddst.norm.test` <-
function(x, base = ddst.base.legendre, c = 100, B=1000, compute.p = F, Dmax = 5, ...) {
 data(ddst.polynomial.fun)
 data(MMnorm12)
 data(MMnorm)

# method.name = as.character(substitute(base)) 
 # only Legendre is implemented yet
 base = ddst.base.legendre
 method.name = "ddst.base.legendre"

n = length(x)
if(n<5) 
	stop("length(x) should be at least 5")
er1 = .Internal(mean(x))
sx  = sort(x)
H   = qnorm((1:n - 3/8)/(n+1/4))
er2 = .Internal(mean((sx[-1] - sx[-n])/(H[-1] - H[-n])))
pp   = (x-er1)/er2
tmpp = c(.Internal(mean(pp)), (.Internal(mean(pp^2)) - 1)/2)

maxN = max(min(Dmax, length(x)-2, 20),1)
u = numeric(maxN)
for (j in 1:maxN) 
u[j] = ddst.phi(pnorm(x,er1,er2), j, base)

coord = numeric(Dmax)
for (k in 1:Dmax) {
  korekta = u[1:k] - t(MMnorm12[[k]]) %*% tmpp
  coord[k] = t(korekta) %*% MMnorm[[k]] %*% korekta * n
}

l = ddst.IIC(coord, n, c)
attr(l, "names") = "n. coord"
t = coord[l]
attr(t, "names") = "WT*"
result = list(statistic = t, parameter = l, method = "Data Driven Smooth Test for Normality")
result$data.name = paste(paste(as.character(substitute(x)), collapse=""), ",   base: ", method.name, ",   c: ", c, sep="")
class(result) = "htest"
if (compute.p) {
tmp = numeric(B)
for (i in 1:B) {
y = rnorm(n)
tmpC = ddst.norm.Nk(y, base, Dmax = Dmax, n=length(y))
l = ddst.IIC(tmpC, n, c)
tmp[i] = tmpC[l]
}
p.val = .Internal(mean(tmp > t))
result$p.value = p.val  
}
result
}


.First.lib <- function(lib, pkg)
{
     library.dynam("orsk", pkg, lib)
     require("BB")
     vers <- library(help=orsk)$info[[1]]
     vers <- vers[grep("Version:",vers)]
     vers <- rev(strsplit(vers," ")[[1]])[1]
     cat("Loaded orsk",vers,"\n")
}

orsk <-
function(x, y, a, al, au, level=0.95, method=c("optim","grid"), d=1e-4){
call <- match.call()
method <- match.arg(method)
if(x < 0 || y < 0) stop("x and y must be positive integers\n")
if(a < 0 || al < 0 || au < 0) stop("odds ratio or confidence intervals must be positive integers\n")
if(level <= 0 || level >= 1) stop("confidence level must be between 0 and 1\n")
if(d <= 0) stop("thereshold value must be positive\n")
q <- qnorm(1-(1-level)/2)
se = log(au/a)/q
if(method=="grid"){
z <- .Fortran("oddsratio",
     x=as.integer(x),
     y=as.integer(y),
     u=as.double(matrix(1, x-1, y-1)),
     w=as.double(matrix(0, x-1, y-1)),
     a=as.double(a),
     al=as.double(al),
     au=as.double(au),
     q=as.double(q),
     d=as.double(d),
     m=as.integer(0),
     wind=as.integer(matrix(0, nrow=(x-1)*(y-1), ncol=2)),
     package="orsk")
w2 <- matrix(z$w, byrow=FALSE, nrow=x-1, ncol=y-1)
windex <- matrix(z$wind, byrow=FALSE, nrow=(x-1)*(y-1), ncol=2)
windex <- subset(windex, windex[,1] > 0)
m <- dim(windex)[1]
if(m == 0) stop("The threshold value d is too small\n")
### estimated odds ratio
n01 <- windex[,1]
n11 <- windex[,2]
n10 <- y-n11
n00 <- x-n01
or.est <- n11*n00/(n10*n01)
s <- sqrt(1/windex[,1] + 1/windex[,2] + 1/(x-windex[,1])+1/(y-windex[,2]))
s <- exp(q*s)
res.lo <- or.est/s
res.up <- or.est*s
rr <- (windex[,2]/y)/(windex[,1]/x)
#standard error (SE) log relative risk
slrr <- sqrt(1/windex[,2]-1/y+1/windex[,1]-1/x)
#lower limit= the exponential of (log(Rel risk)-(1.96*SElogR))
#upper limit= the exponential of (log(Rel risk)+(1.96*SElogR)) 
rr.lo <- exp(log(rr) - q*slrr)
rr.up <- exp(log(rr) + q*slrr)
res <- cbind(windex[,1], x-windex[,1], windex[,2], y-windex[,2], w2[windex], or.est, res.lo, res.up, rr, rr.lo, rr.up)
colnames(res) <- c("ctr_yes","ctr_no","trt_yes","trt_no","SS","OR","OR_lower","OR_upper", "RR", "RR_lower","RR_upper") 
### sort R^2=w2
z1 <- .Fortran("sortw",
       x=as.integer(dim(res)[1]),
       w=as.double(as.vector(w2[windex])),
       ind=as.integer(rep(0, dim(res)[1])),
       package="orsk")
if(length(z1$ind) == 1)
res <- as.data.frame(res)
else
res <- as.data.frame(res[z1$ind,])
}
else{
fw <- function(p){
n11 <- p[2]
n01 <- p[1]
n10 <- y - n11; n00 <- x - n01;
v <- n11*n00/(n10*n01) ### compare to odds ratio
f <- v - a ### compare to odds ratio
s <- exp(q*sqrt(1/n11 + 1/n10 + 1/n01 + 1/n00))
if(n11 <= 0 || n10 <= 0 || n01 <= 0 || n00 <= 0)
print(c(n11, n10, n01, n00))  
tmp1 <- v/s          ### compute confidence interval lower or upper bound
return(f^2+(tmp1-al)^2)
}

set.seed(123)
tmp <- 1:max(2, ceiling(0.1*min(x-1, y-1))) # select number of random integer
p0 <- cbind(sample(y-1)[tmp], sample(x-1)[tmp])
ans <- multiStartNew(par=p0, fn=fw, lower=c(1,1), upper=c(x-1,y-1), action="optimize")
pmat <- ans$par[ans$conv, ]
pmat <- round(pmat)
pmat <- pmat[!duplicated(pmat),]
pmat <- matrix(pmat, ncol=2)
pmat <- cbind(pmat, apply(pmat, 1, fw))
ord1 <- order(pmat[, 3])
ans <- pmat[ord1, ]
ans <- matrix(ans, ncol=3)
n11 <- ans[,2]; n01 <- ans[,1];
n10 <- y - n11; n00 <- x - n01;
or.est <- n11*n00/(n10*n01) ### compare to odds ratio
s <- exp(q*sqrt (1/n11 + 1/n10 + 1/n01 + 1/n00)) 
 res.lo <- or.est / s 
 res.up <- or.est * s
rr <- (n11/y)/(n01/x)
slrr <- sqrt(1/n11-1/y+1/n01-1/x)
#lower limit= the exponential of (log(Rel risk)-(1.96*SElogR))
#upper limit= the exponential of (log(Rel risk)+(1.96*SElogR)) 
rr.lo <- exp(log(rr) - q*slrr)
rr.up <- exp(log(rr) + q*slrr)
res <- cbind(n01,n00, n11, n10, ans[,3], or.est, res.lo, res.up, rr, rr.lo, rr.up)
colnames(res) <- c("ctr_yes","ctr_no","trt_yes","trt_no","SS","OR","OR_lower","OR_upper", "RR", "RR_lower","RR_upper") 
tmp <- which(ans[,3] <= d)
if(length(tmp) == 0)
stop("The threshold value d is too small\n")
else{
res <- as.data.frame(res)[tmp,]
}
}
RET <- list(x=x, y=y, a=a, al=al, au=au, method=method, d=d, res=res) 
RET$call <- call
class(RET) <- "orsk"
return(RET)
}

#plot.orsk <- function(x, type=c("or","rr"), ...) {
#  type <- match.arg(type)
#if(type=="rr")
#hist(x$RR, xlab="relative risk", main="Histogram of relative risk")
#else
#hist(x$OR, xlab="odds ratio", main="Histogram of odds ratio")
#}

print.orsk <- function(x, ...) {
  cat("\n")
  cat("\t Converting odds ratio to relative risk\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Method: ", x$method, "\n")
  cat("Threshold value: ", x$d, "\n")
  cat("\n")
  invisible(x)
}

summary.orsk <- function(object, ...) {
  x <- object
  cat("\n")
  cat("\t Converting odds ratio to relative risk\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Method: ", x$method, "\n")
  cat("Threshold value: ", x$d, "\n")
  cat(paste("The reported odds ratio: ",x$a, ", confidence interval ", x$al, ", ", x$au, "\n", sep=""))
  cat("The estimated results. The calculated odds ratios and relative risks are for \n the scenarios created with different numbers of events in both control and \n treatment group that lead to comparable results for the reported odds ratio \n and confidence interval.\n")
  print(x$res, digits=3)
  cat("\n")
  invisible(x)
}


setwd("~/Dropbox-Personal/Students and Collaborators/Parijat Chakraborty/")
#
# Pareto combination test
#

Pareto.scale <- function(y){
  #
  # Convention: extremes are transformed into extremes
  #
  d = dim(y)[2];
  n = dim(y)[1];
  x = c();
  for (i in c(1:d)){
    x = cbind(x,(n+1)/rank(-y[,i],))
  }
  return(x)
}

Cauchy.scale <- function(y){
  #
  # Convention: extremes are transformed into extremes
  #
  d = dim(y)[2];
  n = dim(y)[1];
  x = c();
  for (i in c(1:d)){
    x = cbind(x,tan(pi*(rank(y[,i],)/(n+1)-1/2)))
  }
  return(x)
}


experiment.multivariate.t<- function(n=1e5, d = 10, nu=1,
                                     w = rep(1,d)/d,
                                     sigma.half=diag(rep(1,d)),
                                     t = n^(0.5)*c(1:100)/100,to.plot=TRUE
                                     ){
 #
 # Simulating multivariate t-distribution variates
 #
  Y = matrix(rnorm(d*n),n,d) %*% sigma.half;
  G= sqrt(rgamma(n,nu/2,rate=1/2)/nu);
  Y = Y/(outer(G,rep(1,d)));  
  
  X.pareto = Pareto.scale(Y);
  T.pareto = X.pareto %*% w;
  
  X.cauchy = Cauchy.scale(Y);
  T.cauchy =  X.cauchy %*% w;
  
  nt = length(t);
  P.pareto = c();
  P.cauchy = c();
  pareto.test = c();
  for (it in c(1:nt)){
    pareto.test =  c(pareto.test, mean(T.pareto>t[it]));
    P.pareto = c(P.pareto,
               mean(T.pareto>t[it])/mean(X.pareto[,1]>t[it]));
    P.cauchy = c(P.cauchy,
                 mean(T.cauchy>t[it])/mean(X.cauchy[,1]>t[it]));
    
  }
  
  if (to.plot){
    plot(t,P.pareto,type="l",ylim=c(-0.5,1.2),cex=0.75,lwd=2,
         main = paste0("Multivariate t-distribution: df=nu=",nu,
                       ", d=",d,
                       " n=",n),
         xlab="1/alpha",ylab="combined p-value/alpha");
    lines(t,P.cauchy,lty=2,col="red",lwd=2);
    Cauchy.limit = sqrt(sum((t(w)%*%sigma.half)^2))/
      sqrt(sum(sigma.half[,1]^2));
    abline(h=Cauchy.limit,col="green",lty=4)
    legend("bottomleft",cex=0.5, legend=c("Pareto comb test","Cauchy comb test",
                                          "Cauchy limit (valid for nu=1 only)"),
           col=c("black","red","green"),lwd=c(2,2,1),lty=c(1,1,4));
    abline(h=c(0,1),lty=3,col="blue")
  }
  return(list("t"=t, "P.pareto"=P.pareto,"P.cauchy"=P.cauchy,
              "pareto.test"=pareto.test))
}

save.some.plots <-function(){
  png(filename = "t-distribution-comparison%d.png");
  out = experiment.multivariate.t(n=10^6, nu=1);
  out = experiment.multivariate.t(n=10^6, nu=4);
  out = experiment.multivariate.t(n=10^6, nu=40);
  out = experiment.multivariate.t(n=10^6, nu=400);
  dev.off()
}

mvt.effect.of.dimension.and.nu <- function(n=10^6,
                                           dims = c(10,50,100,500),
                                           nus  = c(1,4,40,400)){
  library(data.table)
  results = data.table();
  for (d in dims){
    cat("\nDimension d=",d," nu =")
    for (nu in nus){
      cat(" ",nu)
      out = experiment.multivariate.t(n=n,d=d,nu=nu,to.plot = T);
      results = rbind(results,
                      data.table("Pareto"=out$P.pareto,
                                 "Cauchy"=out$P.cauchy,
                                 "t"=out$t,
                                 "nu"=nu,
                                 "d"=d))
    }
  }
  return(results)
}

save.some.results <- function(){
  res=mvt.effect.of.dimension.and.nu() 
  save(file="mvt_cauchy_pareto_effect_of_dim_and_nu.RData","res")
}

plotting.the.results <- function(fname="mvt_cauchy_pareto_effect_of_dim_and_nu.RData",
                                 outfile="mvt_cauchy_pareto_effect_of_dim_and_nu.pdf"){
 load(fname)
    dims = unique(res$d)
    nus = unique(res$nu)
    pdf(file = outfile,
        paper="special",
        width=length(dims)*3,
        height=length(nus)*3)
    par(mfrow=c(length(dims),length(nus)))
    for(d in dims){
      for (nu in nus){
        idx = which((res$d==d)&(res$nu==nu))
        boxplot(res[idx,c("Pareto","Cauchy")],
                main=paste("d=",d,"nu=",nu),xlab="",
                ylim = c(0,1.2));
        abline(h=1,col="blue",lty=2)
        }
    }
    dev.off()
}

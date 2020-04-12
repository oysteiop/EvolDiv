plot_GD=function(G=out$G, D=out$D, species="Species", xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL){
  
# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(D, Beta = g_ev)$e
  
d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(D, Beta = d_ev)$e
  
# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
#beta_g
r2_g = summary(mg)$r.squared
#r2_g
  
md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
#beta_d
r2_d = summary(md)$r.squared
#r2_d
  
# Plot
x11(width=5, height=5)

if(is.null(xmin)){xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))}
if(is.null(xmax)){xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))}
if(is.null(ymin)){ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))}
if(is.null(ymax)){ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))}
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main=species, las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(G)), log10(diag(D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

outlist=data.frame(c(betaG=beta_g, r2G=r2_g, betaD=beta_d, r2D=r2_d))
names(outlist)="value"

return(outlist)
}

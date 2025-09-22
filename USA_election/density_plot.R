# Packages
library(tidyverse)
library(here)
#-------------------------------------------------------------------------------
# Load 'dWeibullT' function
source(here("USA_election", "R", "funcoesWeibullT.R"))
#-------------------------------------------------------------------------------
# Densities
#-------------------------------------------------------------------------------
# Density1
x <- seq(-1, 1, length.out = 500)
y <- dWeibullT(y=x, phi=1.2, m=-0.7)
#nome2 <- paste("USA_election/figures/density1",".pdf",sep="")
#pdf(file = nome2, width = 5, height = 5,family = "Times")
#png("USA_election/figures/density1.png", width = 800, height = 600)
par(mar=c(2.8, 2.7, 1.1, 1))
par(mgp=c(1.7, 0.45, 0))
plot(x, y, type="l", lty=1,lwd=2, xlab="y",
     ylab="f(y)", ylim = c(0,4.8), col = 1)
curve(dWeibullT(y=x, phi=1.5, m=-0.7),col=2,lty=2,lwd=2,add=T)
curve(dWeibullT(y=x, phi=2.0, m=-0.7),col=3,lty=3,lwd=2,add=T)
curve(dWeibullT(y=x, phi=2.5, m=-0.7),col=4,lty=4,lwd=2,add=T)
curve(dWeibullT(y=x, phi=3.0, m=-0.7),col=8,lty=5,lwd=2,add=T)
curve(dWeibullT(y=x, phi=4.0, m=-0.7),col=6,lty=6,lwd=2,add=T)
legend(x=0., y=4.8, 
       c(bquote(plain(m)==-0.7~"," ~ plain(phi) == .(sprintf("%.1f", 1.2))),
         bquote(plain(m)==-0.7~"," ~ plain(phi) == .(sprintf("%.1f", 1.5))),
         bquote(plain(m)==-0.7~"," ~ plain(phi) == .(sprintf("%.1f", 2.0))),
         bquote(plain(m)==-0.7~"," ~ plain(phi) == .(sprintf("%.1f", 2.5))),
         bquote(plain(m)==-0.7~"," ~ plain(phi) == .(sprintf("%.1f", 3.0))),
         bquote(plain(m)==-0.7~"," ~ plain(phi) == .(sprintf("%.1f", 4.0)))),
       lty=c(1,2,3,4,5,6),lwd=c(1,1,1,1,1,1),bty="n",col=c(1,2,3,4,8,6))
#box();axis(1) 
#dev.off()
#-------------------------------------------------------------------------------
# Density2
x <- seq(-1, 1, length.out = 500)
y <- dWeibullT(y=x, phi=2.5, m=0.0)

#nome2 <- paste("USA_election/figures/density2",".pdf",sep="")
#pdf(file = nome2, width = 5, height = 5,family = "Times") 
#png("USA_election/figures/density2.png", width = 800, height = 600)
par(mar=c(2.8, 2.7, 1.1, 1))
par(mgp=c(1.7, 0.45, 0))
plot(x, y, type="l", lty=1,lwd=2, xlab="y",
     ylab="f(y)", ylim = c(0,4.8), col = 1)
curve(dWeibullT(y=x, phi=3.0, m=0.0),col=2,lty=2,lwd=2,add=T)
curve(dWeibullT(y=x, phi=4.5, m=0.0),col=3,lty=3,lwd=2,add=T)
curve(dWeibullT(y=x, phi=6.0, m=0.0),col=4,lty=4,lwd=2,add=T)
curve(dWeibullT(y=x, phi=10.0, m=0.0),col=8,lty=5,lwd=2,add=T)
curve(dWeibullT(y=x, phi=12.0, m=0.0),col=6,lty=6,lwd=2,add=T)
legend(x=0.1, y=4.6, 
       c(bquote(plain(m) == .(sprintf("%.1f", 0.0)) ~ "," ~ plain(phi) == .(sprintf("%.1f", 2.5))),
         bquote(plain(m) == .(sprintf("%.1f", 0.0)) ~ "," ~ plain(phi) == .(sprintf("%.1f", 3.0))),
         bquote(plain(m) == .(sprintf("%.1f", 0.0)) ~ "," ~ plain(phi) == .(sprintf("%.1f", 4.5))),
         bquote(plain(m) == .(sprintf("%.1f", 0.0)) ~ "," ~ plain(phi) == .(sprintf("%.1f", 6.0))),
         bquote(plain(m) == .(sprintf("%.1f", 0.0)) ~ "," ~ plain(phi) == .(sprintf("%.1f", 10.0))),
         bquote(plain(m) == .(sprintf("%.1f", 0.0)) ~ "," ~ plain(phi) == .(sprintf("%.1f", 12.0)))),
       lty=c(1,2,3,4,5,6),lwd=c(1,1,1,1,1,1),bty="n",col=c(1,2,3,4,8,6))
box();axis(1) 
#dev.off()

#-------------------------------------------------------------------------------
# Density3
x <- seq(-1, 1, length.out = 500)
y <- dWeibullT(y=x, phi=3, m=0.5)

#nome2 <- paste("USA_election/figures/density3",".pdf",sep="")
#pdf(file = nome2, width = 5, height = 5,family = "Times") 
#png("USA_election/figures/density3.png", width = 800, height = 600)
par(mar=c(2.8, 2.7, 1.1, 1))
par(mgp=c(1.7, 0.45, 0))
plot(x, y, type="l", lty=1,lwd=2, xlab="y",
     ylab="f(y)", ylim = c(0,3.7), col = 1)
curve(dWeibullT(y=x, phi=5.0, m=0.5),col=2,lty=2,lwd=2,add=T)
curve(dWeibullT(y=x, phi=7.5, m=0.5),col=3,lty=3,lwd=2,add=T)
curve(dWeibullT(y=x, phi=10.0, m=0.5),col=4,lty=4,lwd=2,add=T)
curve(dWeibullT(y=x, phi=12.5, m=0.5),col=8,lty=5,lwd=2,add=T)
curve(dWeibullT(y=x, phi=15.0, m=0.5),col=6,lty=6,lwd=2,add=T)
legend(x=-0.95, y=3.8, 
       c(bquote(plain(m)==0.5~"," ~ plain(phi) == .(sprintf("%.1f", 3.0))),
         bquote(plain(m)==0.5~"," ~ plain(phi) == .(sprintf("%.1f", 5.0))),
         bquote(plain(m)==0.5~"," ~ plain(phi) == .(sprintf("%.1f", 7.5))),
         bquote(plain(m)==0.5~"," ~ plain(phi) == .(sprintf("%.1f", 10.0))),
         bquote(plain(m)==0.5~"," ~ plain(phi) == .(sprintf("%.1f", 12.5))),
         bquote(plain(m)==0.5~"," ~plain(phi) == .(sprintf("%.1f", 15.0)))),
       lty=c(1,2,3,4,5,6),lwd=c(1,1,1,1,1,1),bty="n",col=c(1,2,3,4,8,6))
box();axis(1) 
#dev.off()

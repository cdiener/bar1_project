N.A = 6.0221415e23

fit = as.numeric( t(read.table("ver_fit.txt", header=0)) )*1e-24*N.A
data = as.numeric( t(read.table("ver_data.txt", header=0)) )*1e-24*N.A

ord = order(fit)
svg("ver.svg", width=5, height=8)
plot(rep(0,78), data, xlim=c(-0.5,1.5), xaxt="n", col="black", pch=19, 
	ylab=expression(paste(alpha-plain(factor)," [molecs/",mu,m^2,"]",sep=" ")), main="", xlab="")
points( rep(1,78), fit, col="blue", pch=17)
dev.off()



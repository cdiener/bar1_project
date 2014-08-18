############################################################
# alpha_fit.R
# This file is part of Mesh Tools
#
# Copyright (C) 2012 - Christian Diener
#
# Mesh Tools is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Mesh Tools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Mesh Tools; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, 
# Boston, MA  02110-1301  USA
############################################################

source("mesh_tools.R")

SIZE = 100

# Those are the parameters used for config
# if you redefine them after the include and before
# calling write_conf they will be ignored
EC50 = 30.41
NH = 1.62
FMAX = 1507
FMIN = 275.9

# Parameters for Bar1 activity
k0 = 0.2402
k1 = 0.5992


# Values used for mesh element sizes
# Good for fast evals: LCW = 10-16, LCM=0.3-0.5, LCS = 6-8
# Good for accurate: LCW = 6-8, LCM = 0.1-0.25, LCS = 0.5-2

ALPHA = 3
LCW = 6
LCM = 0.5
LCS = 2

# Cell smoothing
SMOOTH = 0.3

IMG = read.table("config.txt", header=T)

AREA = vector( length=nrow(IMG) )
kalpha = vector( length=sum(IMG$date_index==1))
#scale = vector(mode="list", length=N_C)

pvals = vector(length=nrow(IMG))
N_ALPHA = vector(length=sum(IMG$date_index<3))
N_A = vector(length=sum(IMG$date_index<3))

F.real = NULL
F.mod = NULL

# Read ares ratios beforehand to get a mean ratio
for(i in 1:nrow(IMG))
{
	d_idx = IMG$date_index[i]
	p_idx = IMG$pos_index[i]

	cells <- read.cellID( sprintf("out_%d_%d",d_idx,p_idx) , RFP=T, remove.irregular=F )
	AREA[i] = sum(cells$area)/IMG$width[i]^2
}

for(i in 1:nrow(IMG))
{
	PIXELS = IMG$pix[i]
	SIZE = IMG$width[i]
	PIX_PER_UM <- PIXELS/SIZE

	d_idx = IMG$date_index[i]
	p_idx = IMG$pos_index[i]

	cells <- read.cellID( sprintf("out_%d_%d",d_idx,p_idx) , RFP=T, remove.irregular=F )
	bl = read.boundary( sprintf("BF_D0%d_P0%d.tif_BOUND.txt", d_idx, p_idx) )

	MATAL <- cells$ID[cells$RFP>RFP_CUT]+1
	MATA <- cells$ID[-MATAL]
	scale = AREA[i]
	
	dmap <- generate.dmap.poly(cells$x, cells$y, cells, SIZE)

	surfsa <- interpolate.surface(cells, bl)
	surfs <- get.surfaces(cells$ID+1, surfsa, do.smooth=SMOOTH, nP=128)
	
	write( sprintf( "Working with %d MATa and %d MATalpha cells...", length(MATA), length(MATAL) ), file="" )
	write( sprintf( "Cells areas: %g +- %g (total %g%% of image)", mean(cells$area), sd(cells$area), 100*AREA[i] ), file="" )
	write( sprintf( "World radius: %g um", SIZE ), file="" )

	write.geo.poly(cells, surfs, filename=sprintf("cells_%d_%d.geo", d_idx, p_idx) )
	write("Building mesh...", file="")
	system( sprintf("gmsh cells_%d_%d.geo -algo front2d -o cells_%d_%d.msh -2 -v 0", d_idx, p_idx, d_idx, p_idx ) )
	write("Done.", file="")
	
	if( d_idx < 3)	
	{
		P_INIT = c(10.0, k0*scale, k1*scale, 0)
		write.conf(cells, filename=sprintf("conf_%d_%d.txt", d_idx, p_idx) )
		#if(p_idx[i]<2)#!file.exists(sprintf("sol_%d_%d.vtu", d_idx[i], p_idx[i] )) )
		#{		
		write("Fitting...", file="")
	
		system( sprintf("./fit_bfgs cells_%d_%d.msh conf_%d_%d.txt", d_idx, p_idx, d_idx, p_idx) )
		system( sprintf("mv fit_result.txt fit_%d_%d.txt", d_idx, p_idx ) )
		system( sprintf("mv solStat.vtu sol_%d_%d.vtu", d_idx, p_idx ) )
		#}
		kalpha[i] = read.table(sprintf("fit_%d_%d.txt",d_idx, p_idx), header=F)[1,1]
		write("Done!\n", file="")
	}
	else
	{		
		sigs = kalpha[1:(i-1)]
		sigs = sigs[pvals[1:(i-1)]<0.05]
		P_INIT = c(mean(sigs), k0*scale, k1*scale, 0)
		write.conf(cells, filename=sprintf("conf_%d_%d.txt", d_idx, p_idx) )
	
		write("Testing...", file="")
	
		a = system( sprintf("./stationary cells_%d_%d.msh conf_%d_%d.txt", d_idx, p_idx, d_idx, p_idx), intern=TRUE )
		a = a[(length(a)-length(MATA)+1):length(a)]
		a = as.numeric(a)

		out = file( sprintf("fit_%d_%d.txt", d_idx, p_idx ), open="w" )		
		sapply(P_INIT, write, file=out)		
		sapply(a, write, file=out)
		close(out)
		
		system( sprintf("mv solStat.vtu sol_%d_%d.vtu", d_idx, p_idx ) )
		write("Done!\n", file="")
	}
	
	# Calculate goodness of fit...
	conf = read.table(sprintf("conf_%d_%d.txt", d_idx, p_idx), header=F)
	fit = read.table(sprintf("fit_%d_%d.txt", d_idx, p_idx), header=F)

	N_ALPHA[i] = conf[1,1]
	N_A[i] = conf[2,1]
	A = fit[-(1:4),1]
	F.fit = A
	F.ref = conf[-(1:6),1]
	rss.ref = sum(F.ref^2)
	rss.fit = sum( (F.ref - F.fit)^2 )
	pvals[i] = 1-pf( (rss.ref-rss.fit)/1/rss.fit*(N_A[i]-1), 1, N_A[i]-1 )

	F.real = c(F.real, F.ref)
	F.mod = c(F.mod, F.fit)
}

save(kalpha, AREA, pvals, file="alpha.dat")



write(sprintf("%d out of %d fits significant.", sum(pvals[IMG$date_index<3]<0.05), length(pvals[IMG$date_index<3])), file="")
write(sprintf("Fitted %d MATalpha cells and %d MATa cells.", sum(N_ALPHA[pvals[IMG$date_index<3]<0.05]), sum(N_A[pvals[IMG$date_index<3]<0.05])), file="")
write(sprintf("KALPHA = %g +- %g", mean(kalpha[pvals[IMG$date_index<3]<0.05]), sd(kalpha[pvals[IMG$date_index<3]<0.05])), file="")
write(sprintf("Verification out of sample p-value: %g", pvals[IMG$date_index==3]), file="")

r_casy = 2.25
N_A = 6.02214129e23
jalpha.sig = kalpha[pvals[IMG$date_index<3]<0.05]
jalpha.tot = jalpha.sig*2*pi*r_casy
ksec = jalpha.tot*1e-9*1e-15*N_A

save(ksec, file="sec.dat")

pdf("alpha_sec.pdf", width=6, height=5)
dens = density(ksec, bw=150) 
plot(x=dens$x, 100*dens$y/max(dens$y), type="l", lwd=2,
	xlab=expression(paste(alpha,"-factor secretion [molecules/s per cell]")),
	ylab=" normalized density [%max]", main="")
rug(ksec, lwd=1)
dev.off()


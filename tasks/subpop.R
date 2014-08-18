############################################################
# subpop.R
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

# This file quantifies the effects on varying density of mixed cell populations
# on the gradients, maximum information content and coefficient of variation of the 
# alpha-factor distribution

require(entropy)
source("mesh_tools.R")

# These parameters have to be changed if your image has different measures!

SIZE = 200
WIDTH = 100

# Those are the parameters used for config
# if you redefine them after the include and before
# calling write_conf they will be ignored
EC50 = 30.41
NH = 1.62
FMAX = 1507
FMIN = 275.9

# Alpha secretion rate
kalpha = 101.589

# Parameters for Bar1 activity
k0 = 0.2402
k1 = 0.5992

# Values used for mesh element sizes
# Good for fast evals: LCW = 10-16, LCM=0.3-0.5, LCS = 6-8
# Good for accurate: LCW = 6-8, LCM = 0.1-0.25, LCS = 0.5-2

ALPHA = 3
LCW = 50
LCM = 1
LCS = 4

# Mean and sd of haploid cell radii
R_MEAN = 1.982
R_SD = 0.251

N_CELLS = c(100)
N_SMALL = 20
N_REP = 4
N_P = 512

OFFSCREEN = ""

get.cv = function(sol_mat)
{
	af = as.numeric(sol[sol>0])
	cv = sd(af)/mean(af)

	return(cv)
}

get.ent = function(sol_mat)
{
	af = as.numeric(sol[sol>0])
	g.af = cut(af, seq(0, max(af), length.out=ceiling(log(length(af),2)+1)), include.lowest=1)

	ent = entropy.shrink(table(g.af))	

	return(ent)
}

get.grad = function(sol_mat)
{
	grad = matrix(nrow = nrow(sol_mat), ncol = ncol(sol_mat))
	h = 2*ceiling(SIZE/1.5)/(N_P-1)

	for(i in 1:nrow(sol_mat))
		for(j in 1:ncol(sol_mat))
		{
			if(sol_mat[i,j] == 0)
			{ 
				grad[i,j] = NA
				next
			}
			
			if(j==ncol(sol_mat) || sol_mat[i,j+1]==0) dx = sol_mat[i,j] - sol_mat[i,j-1]
			else dx = sol_mat[i,j] - sol_mat[i,j+1]
			
			if(i==1 || sol_mat[i-1,j]==0) dy = sol_mat[i+1,j] - sol_mat[i,j]
			else dy = sol_mat[i,j] - sol_mat[i-1,j]

			grad[i,j] = sqrt( (dx/h)^2 + (dy/h)^2 )
		}

	mu = max(sol_mat)

	return(grad/mu)	
}

get.maxmin = function()
{
	
	get.stats = function(m)
	{
		i.min = which.min(m$alpha)
		i.max = which.max(m$alpha)
		
		d = sqrt( (m$x[i.min]-m$x[i.max])^2 + (m$y[i.min]-m$y[i.max])^2 )

		return( c(m$alpha[i.min], m$alpha[i.max], d) )
	} 
	
	surf_alpha = read.table("absorbed", header=T)
	
	res = by(surf_alpha, surf_alpha$CellId, get.stats)
	res = t( sapply(res, identity) )

	return(res )
}

wt = list(
ent = matrix(nrow=N_REP, ncol=length(N_CELLS)),
alpha = matrix(nrow=N_REP, ncol=length(N_CELLS)),
maxmin = vector( mode="list", length=length(N_CELLS)),
d.maxmin = vector( mode="list", length=length(N_CELLS))
)

bar1D = list(
ent = matrix(nrow=N_REP, ncol=length(N_CELLS)),
alpha = matrix(nrow=N_REP, ncol=length(N_CELLS)),
maxmin = vector( mode="list", length=length(N_CELLS)),
d.maxmin = vector( mode="list", length=length(N_CELLS))
)

for(ic in 1:length(N_CELLS))
{
	nc = N_CELLS[ic]
	
	for(nr in 1:N_REP)
	{
		write("Sampling cells...", file="")
		cells = data.frame( x=rep(0,nc+N_SMALL), y=rep(0,nc+N_SMALL), r=rnorm(nc+N_SMALL, R_MEAN, R_SD) )
		cells[1,1:2] = runif(2, R_MEAN, WIDTH-R_MEAN)
		for(i in 2:nc) cells[i, 1:2] = sample.valid(i, cells, lim=c(R_MEAN, WIDTH-R_MEAN))
		for(i in (nc+1):(nc+N_SMALL)) cells[i, 1:2] = sample.valid(i, cells, lim=c(-WIDTH+R_MEAN, -WIDTH/2+R_MEAN))
		cells = cells[c( (nc+N_SMALL/2+1):(nc+N_SMALL), 1:(nc+N_SMALL/2)), ]
		AREA = sum(pi*cells$r^2)/SIZE^2

		# For wt...
		P_WT = c(kalpha, k0*AREA, k1*AREA, 0)
		
		write("Building mesh...", file="")
		dmap <- generate.dmap.poly(cells$x, cells$y, cells, SIZE)
		write.geo(cells, filename=sprintf("cells_%d_%d.geo", nc, nr))
		system( sprintf("gmsh cells_%d_%d.geo -algo front2d -o cells_%d_%d.msh -2 -v 0", nc, nr, nc, nr) )

		conf = file(sprintf("wt_%d_%d.txt", nc, nr), open='w')
		write((nc+N_SMALL)/2, file=conf)
		write((nc+N_SMALL)/2, file=conf)
		sapply(P_WT, write, file=conf)
		close(conf)

		write("Solving for wt...", file="")
		system( sprintf("./stationary cells_%d_%d.msh wt_%d_%d.txt", nc, nr, nc, nr) )
		system( sprintf("mv solStat.vtu wt_%d_%d.vtu", nc, nr) )

		system( sprintf("%s ./grid_to_m.py wt_%d_%d.vtu alpha %d %d", OFFSCREEN, nc, nr, ceiling(SIZE/1.5), N_P) )

		write("Optimizing measures...", file="")
		sol = as.matrix( read.table(sprintf("alpha_%d.txt",ceiling(SIZE/1.5)), header=0 ) )
		wt$ent[nr,ic] = get.ent(sol)
		wt$alpha[nr,ic] = mean(sol[sol>1e-6])
		cell.stats = get.maxmin()
		if(nc<4) 
		{ 
			maxmin = cell.stats[2]/cell.stats[1]
			d.maxmin = cell.stats[3]
        }
		else 
		{
			maxmin = cell.stats[,2]/cell.stats[,1]
			d.maxmin = cell.stats[,3]
		}
		wt$maxmin[[ic]] = c(wt$maxmin[[ic]], maxmin)
		wt$d.maxmin[[ic]] = c(wt$d.maxmin[[ic]], d.maxmin)

		plot.profile(cells, NULL,  sprintf("alpha_%d.txt",ceiling(SIZE/1.5)), in.molec=1, out_file=sprintf("images/wt_%d_%d.tif", nc, nr), width=ceiling(SIZE/1.5))

		# For bar1D
		P_BAR1D = c(kalpha, 1e-4, 0, 0)

		conf = file(sprintf("bar1D_%d_%d.txt", nc, nr), open='w')
		write((nc+N_SMALL)/2, file=conf)
		write((nc+N_SMALL)/2, file=conf)
		sapply(P_BAR1D, write, file=conf)
		close(conf)

		write("Solving for Bar1D...", file="")
		system( sprintf("./stationary cells_%d_%d.msh bar1D_%d_%d.txt", nc, nr, nc, nr) )
		system( sprintf("mv solStat.vtu bar1D_%d_%d.vtu", nc, nr) )

		system( sprintf("%s ./grid_to_m.py bar1D_%d_%d.vtu alpha %d %d", OFFSCREEN, nc, nr, ceiling(SIZE/1.5), N_P) )

		write("Optimizing measures...", file="")
		sol = as.matrix( read.table(sprintf("alpha_%d.txt",ceiling(SIZE/1.5)), header=0 ) )
		bar1D$ent[nr,ic] = get.ent(sol)
		bar1D$alpha[nr,ic] = mean(sol[sol>1e-6])
		cell.stats = get.maxmin()
		if(nc<4) 
		{ 
			maxmin = cell.stats[2]/cell.stats[1]
			d.maxmin = cell.stats[3]
        }
		else 
		{
			maxmin = cell.stats[,2]/cell.stats[,1]
			d.maxmin = cell.stats[,3]
		}
		bar1D$maxmin[[ic]] = c(bar1D$maxmin[[ic]], maxmin)
		bar1D$d.maxmin[[ic]] = c(bar1D$d.maxmin[[ic]], d.maxmin)

		plot.profile(cells, NULL, sprintf("alpha_%d.txt",ceiling(SIZE/1.5)), in.molec=1, out_file=sprintf("images/bar1D_%d_%d.tif", nc, nr), width=ceiling(SIZE/1.5))
		
	}
}

save(wt, bar1D, file="stats.dat")

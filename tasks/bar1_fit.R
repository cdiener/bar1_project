### Analysis file for cell positions and mesh generation

GFP <- 0
RFP <- 1
PIXELS = 512
SIZE = 70.14
REC_TIME = 180

PIX_PER_UM <- PIXELS/SIZE
RFP_CUT <- 4

P_INIT = c(100, 2.8, 1.4, 0.028)
EC50 = 50
NH = 1.45
EC50WT = 236.8
NHWT = 0.842
FMAX = 38.31
FMIN = 9.8
A_C = 10^seq(1,3, length.out=16)

CONC = c()

H_HEAD = 
"#ifndef CELL_POS_INCLUDED\n
#define CELL_POS_INCLUDED\n
/******** Position header file ********/\n\n
" 

# Values used for mesh elemt sizes
# Good for fast evals: A = 3, W = 10-16, LCM=0.3-0.5, LCS = 6-8
# Good for accurate: A = 3, W = 6-8, LCM = 0.1-0.25, LCS = 0.5-2

ALPHA = 3
LCW = 16
LCM = 0.5
LCS = 6

GEO_HEAD ='' #sprintf("lc = %0.3f;\n", LC)

generate.dmap.poly = function(xy.matrix, cells)
{   
    dens <- KernSmooth::bkde2D(xy.matrix, bandwidth=rep(mean(cells$r),2), gridsize=c(128,128))

    mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2

    return(list(xbin=mkBreaks(dens$x1), ybin=mkBreaks(dens$x2), d=dens$fhat/max(dens$fhat) )) 
}

get.dens = function(x,y, dm=dmap) 
{
    xbin <- cut(x, dm$xbin, labels = FALSE)
    ybin <- cut(y, dm$ybin, labels = FALSE)
 
    return(dm$d[xbin,ybin])
}

point = function(idx,x,y, d=dmap) sprintf( 'Point(%i) = {%0.3f, %0.3f, 0.0, %0.5f};' , idx, x, y, max(exp(-get.dens(x,y,d)*ALPHA)*LCS, LCM) )
point.out = function(idx,x,y, lc=LCW) sprintf( 'Point(%i) = {%0.3f, %0.3f, 0.0, %0.5f};' , idx, x, y, lc )
line = function(idx,p1,p2) sprintf( 'Line(%i) = {%i, %i};', idx, p1, p2 )
spline = function(idx, plist) sprintf( 'Spline(%i) = {%s};', idx, paste(plist, collapse=', ') )
circle = function(idx, p1, p2, p3) sprintf( 'Circle(%i) = {%i, %i, %i};', idx, p1, p2, p3 )
loop = function(idx, l1, l2, l3, l4) sprintf( 'Line Loop(%i) = {%i, %i, %i, %i};', idx, l1, l2, l3, l4)


write.geo.poly = function(cells, surfs ,filename='cells_poly.geo', nP=128, rad=SIZE)
{
    geo <- file(filename, open='w')
    write(GEO_HEAD, file=geo)
    pi <- 0
    li <- nrow(cells)+1
    
    #MATAL <- cells$ID[cells$RFP>RFP_CUT]+1
    #MATA <- cells$ID[cells$RFP<=RFP_CUT]+1
    if(length(MATAL)>0) surfs <- surfs[c(MATAL, MATA)]
    else surfs <- surfs[MATA]
    
    
    geo.cell.poly = function(i, surfs, nP, out_file)
    {
        write('\n', file=out_file)
        coord <- get.surface(i, surfs, do.smooth=T, nP=nP)
        n <- length(coord$x)
        
        sapply(1:n, function(k) write( point(pi+k, coord$x[k], coord$y[k]), file=out_file) )
        
        write( spline(li+1, c(pi+1:n, pi+1)), file=out_file )
   
        write( sprintf( 'Line Loop(%i) = {%i};', li+2, li+1), file=out_file )
        write( sprintf( 'Physical Line(%i) = {%i};', i, li+1 ), file=out_file )
        
        pi <<- pi+n
        li <<- li+2
    }
    
    
    write( point.out(pi+1, 0, 0) ,file=geo )
    write( point.out(pi+2, 0, rad) ,file=geo )
    write( point.out(pi+3, rad, 0) ,file=geo )
    write( point.out(pi+4, 0, -rad) ,file=geo )
    write( point.out(pi+5, -rad, 0) ,file=geo )
    pi <- pi+5
    
    write( circle(li+1, 2, 1, 3), file=geo )
    write( circle(li+2, 3, 1, 4), file=geo )
    write( circle(li+3, 4, 1, 5), file=geo )
    write( circle(li+4, 5, 1, 2), file=geo )
        
    write( loop(li+5, li+1, li+2, li+3 ,li+4), file=geo )
    write( sprintf( 'Physical Line(%i) = {%i, %i, %i, %i};', 0, li+1, li+2, li+3, li+4 ), file=geo )
    li <- li+5
    start <- li
    
    dummy <- sapply(1:length(surfs), geo.cell.poly, surfs, nP, geo)
    
    loops <- paste( seq(start,li, by=2), collapse=', ')
    write( sprintf('\nRuled Surface (%i) = { %s };', li+1, loops), file=geo )
    write( sprintf('Physical Surface(%i) = { %i };', nrow(cells)+1, li+1), file=geo )
    close(geo)
}

h.cell = function(i, cells, out_file, vec, idx)
{
    write( sprintf('%s[%i][0] = %0.5f;', vec, i-1, cells$x[ idx[i] ]), file=out_file )
    write( sprintf('%s[%i][1] = %0.5f;', vec, i-1, cells$y[ idx[i] ]), file=out_file )
    write( sprintf('%s[%i][2] = %0.5f;', vec, i-1, cells$r[ idx[i] ]), file=out_file )
}

write.conf = function(cells, filename='conf.txt')
{
	out = file(filename, open='w')
	
	write(length(MATAL), file=out)
	write(length(MATA), file=out)
	
	sapply(P_INIT, write, file=out)
	F = cells$GFP[(cells$ID+1) %in% MATA]

	F[F<FMIN] = FMIN
	F[1:length(F)] = FMIN + FMAX*(P_INIT[1]^NHWT/(P_INIT[1]^NHWT + EC50WT^NHWT))
	alpha = EC50*( (F-FMIN)/(FMAX-F+FMIN) )^(1/NH)
	
	sapply(alpha, write, file=out)
	
	close(out)
}

### Functions for smooth irregular shapes

read.boundary = function(file)
{
    bound <- read.table(file, header=F)
    bound[,1] <- (bound[,1] - PIXELS/2)/PIX_PER_UM
    bound[,2] <- (PIXELS/2 - bound[,2])/PIX_PER_UM
    
    return( lapply(unique(bound[,3]), function(i) bound[bound[,3]==i,1:2]) )
}

assign.borders = function(rd, bl)
{
    rd.xy <- rd[seq(1,nrow(rd),by=2), c('xpos','ypos')]
    rd$border <- vector(length=nrow(rd))
    closest <- apply(rd.xy, 1, function(xy) 
                                            which.min( apply(bl$mids,1, function(bli) dist(rbind(xy,bli))) )
                    )
    closest <- unique(bl$xy[,3])[closest]
    
    rd$border <- rep(closest ,each=2)
    
    return( list(rd=rd,bl=bl) )
}

interpolate.surface = function(cells, borders)
{
    #Helper functions to interpolate a single surface first
    
    d.e = function(p1,p2){
        p1 <- as.numeric(p1)
        p2 <- as.numeric(p2)
        return( sqrt(sum((p1-p2)*(p1-p2))) )
    }
    angle = function(p, ref, mid)
    {
        pm <- p - mid
        refm <- ref - mid
        
        if(all(pm==refm)) return(0)
        else if(all(pm==-refm)) return(pi)
        else return(  acos(sum(pm*refm)/sqrt(sum(pm*pm))/sqrt(sum(refm*refm)) ) )
    }
    
    inter.single = function(bID, bor)
    {
        surf.points <- bor[[bID+1]]
        surf.mid <- colMeans(surf.points)


        t.p <- surf.points[surf.points[,2]>surf.mid[2],]
        l.p <- surf.points[surf.points[,2]<=surf.mid[2],]

        t.i <- which.min( t.p[ ,2] )
        t.s <- t.p[t.i, ]
        t.ccw <- t.s[1] > surf.mid[1]
        
        l.i <- which.max( l.p[ ,2] )
        l.s <- l.p[l.i, ]
        l.ccw <- l.s[1] < surf.mid[1]
        
        t.a <- apply(t.p, 1, angle, ref=t.s, mid=surf.mid)
        t.r <- apply(t.p, 1, d.e, p2=surf.mid)
        if(!t.ccw) t.a <- pi - t.a
        
        l.a <- apply(l.p, 1, angle, ref=l.s, mid=surf.mid)
        l.r <- apply(l.p, 1, d.e, p2=surf.mid)
        if(l.ccw) l.a <- l.a + pi else l.a <- 2*pi - l.a

        phi <- c(t.a,l.a)
        r <- c(t.r, l.r)
        ix <- sort(phi, index.return=T)$ix
        phi <- phi[ix]
        r <- r[ix]
        #r <- rep(10,length(phi))
        
        out <- list(mid=surf.mid, phi=phi, r=r )
        return(out)
    }   
    
    surf <- list()
    surf[cells$ID+1] <- lapply(cells$ID, inter.single, bor=borders)
    return(surf)
    
}

get.surface = function(idx, surfaces, do.smooth=T, nP = 256)
{
    if(do.smooth)
    {		
        n <- length( surfaces[[idx]]$phi ) 
        x <- surfaces[[idx]]$phi
        x <- c(x[-n]-2*pi,x,2*pi+x[-1])
        y <- surfaces[[idx]]$r
        y <- c(y[-n],y,y[-1])
        smo <- smooth.spline(x,y,spar=0.6)
        #smo <- loess(y ~ x, span=1/nP*6, degree=2)
        phi <- seq(2*pi/nP,2*pi,length.out=nP)
        r <- predict(smo, phi)$y-2*SIZE/PIXELS

		#Resolve olverlap problems
        #r <- predict(smo,newdata=data.frame(x=phi))
    }
    else
    {
        phi <- surfaces[[idx]]$phi
        r <- surfaces[[idx]]$r
    }
    
    return( list(x=surfaces[[idx]]$mid[1]+cos(phi)*r, y=surfaces[[idx]]$mid[2]+sin(phi)*r) )
}


#### MAIN PART ####

N_C = 8
N_R = 3

source("diff.R")
BASE_SCALE = solve(SIZE , REC_TIME)$value/solve(0,REC_TIME)$value
bar1 = vector(length=8)

for(i in 1:length(A_C))
{
	raw_data <- read.table( "out_3_1", header=T)
	raw_data$xpos <- (raw_data$xpos - PIXELS/2)/PIX_PER_UM
	raw_data$ypos <- (PIXELS/2 - raw_data$ypos)/PIX_PER_UM	
	#raw_data <- raw_data[raw_data$fft.stat<0.4, ]	

	bl <- read.boundary( "BF_P03_T01.tif_BOUND.txt" )	

	cells <- matrix(nrow=nrow(raw_data), ncol=5)
	colnames(cells) <- c('ID', 'x', 'y', 'r', 'GFP')

	attach(raw_data)

	cells[, 1] <- cellID
	cells[, 2] <- xpos
	cells[, 3] <- ypos
	cells[, 4] <- .25*( maj.axis + min.axis )/PIX_PER_UM
	cells[, 5] <- (f.tot + f.tot.m3)/(a.tot.m3 + a.tot)/f.bg - 1 

	cells <- data.frame(cells)
	detach(raw_data)

	MATAL <- NULL
	MATA <- cells$ID+1
	P_INIT[1] = A_C[i]
	P_INIT[4] = 0#BASE_SCALE

	#allb <- NULL
	#dummy <- sapply(bl, function(k) allb <<- rbind(allb,k))

	#dmap <- generate.dmap.poly(allb, cells)

	#surfs <- interpolate.surface(cells, bl)

	#write.geo.poly(cells, surfs, nP=64, filename=sprintf("cells_%d_%d.geo", i, j) )
	#write("Building mesh...", file="")
	#system( "gmsh cells_3_1.geo -algo front2d -o cells_3_1.msh -2 -v 0", i, j, i, j ) )
	#write("Done.", file="")

	write.conf(cells, filename=sprintf("conf_%d.txt", i) )

	cols = colorRampPalette(c("Turquoise","White","Red"))

	write( sprintf("Dirichlet: %g", P_INIT[4]), file="" )
	write("Fitting...", file="")
	system( sprintf("./fit_nm cells_3_1.msh conf_%d.txt", i) )
	system( sprintf("mv fit_result.txt fit_%d.txt", i ) )
	bar1[i] = read.table(sprintf("fit_%d.txt",i), header=F)[2,1]
	write("Done!\n", file="")
} 

alpha_loc = vector(length=length(A_C))

for(i in 1:length(A_C))
{
	alpha_loc[i] = read.table(sprintf("conf_%d.txt", i), header=F, skip=6)[1,1]
}

act = alpha_loc^NH/(alpha_loc^NH + EC50^NH)
model = lm(bar1 ~ act) 
data = list(x=act, x=bar1, alpha=A_C)

pdf("kinetics.pdf", width=4, height=3)
plot(alpha_loc, bar1, pch=19, type="b", log="x")
dev.off()

pdf("linear.pdf", width=4, height=3)
plot(act, bar1, pch=19)
abline(coef(model))
dev.off()

save(data, file="data.txt")

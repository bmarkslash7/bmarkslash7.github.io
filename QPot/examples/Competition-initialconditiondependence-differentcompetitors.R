require(QPot); require(phaseR); require(plot3D)

### Focus on Lotka_Volterra competition model
dxstring <- 'rx*x*(kx - x - alphaxy*y)*(1/kx)'
dystring <- 'ry*y*(ky - y - alphayx*x)*(1/ky)'
### First look at the situation where both competitors coexist
model.state.init.cond2 <- c(rx = 4, ry = 4, kx = 10, ky = 15, alphaxy = 1.5, alphayx = 2.1)
dx <- Model2String(dxstring, parms = model.state.init.cond2)
dy <- Model2String(dystring, parms = model.state.init.cond2)


########################################################################
### Step 1: Analyze the deterministic skeleton #########################
require(phaseR)
competition.model <- function(t, y, parameters){
	x <- y[1]
	y <- y[2]
	rx <- parameters["rx"]
	ry <- parameters["ry"]
	alphaxy <- parameters["alphaxy"]
	alphayx <- parameters["alphayx"]
	kx <- parameters["kx"]
	ky <- parameters["ky"]
	dx <- rx*x*(kx - x - alphaxy*y)*(1/kx)
	dy <- ry*y*(ky - y - alphayx*x)*(1/ky)
	list(c(dx,dy))
}
flowField(deriv = competition.model, x.lim = c(0, 20), y.lim = c(0, 20), parameters = model.state.init.cond2, 
		add = F, points = 30, col = "grey70", ann = F, arrow.head = 0.025, frac = 1.1, 
		xaxs = "i", yaxs = "i", las = 1)
supp.print <- nullclines(deriv = competition.model, x.lim = c(0, 20), y.lim = c(0, 20), 
		parameters = model.state.init.cond2, col = c("blue", "red"), points = 250)
traj <- trajectory(competition.model,y0=c(4,1),
				parameters = model.state.init.cond2 , t.end = 250, lwd = 1.5 , pch = 16)
traj <- trajectory(competition.model,y0=c(1,4),
				parameters = model.state.init.cond2 , t.end = 250, lwd = 1.5 , pch = 16)
traj <- trajectory(competition.model,y0=c(20,1),
				parameters = model.state.init.cond2 , t.end = 250, lwd = 1.5 , pch = 16)
traj <- trajectory(competition.model,y0=c(10,15),
				parameters = model.state.init.cond2 , t.end = 250, lwd = 1.5 , pch = 16)



########################################################################
### Step 2: Stochastic simulation ######################################
dxstring <- 'rx*x*(kx - x - alphaxy*y)*(1/kx)'
dystring <- 'ry*y*(ky - y - alphayx*x)*(1/ky)'
model.state.init.cond2 <- c(rx = 4, ry = 4, kx = 10, ky = 10, 
							alphaxy = 1.5, alphayx = 2.1)
dx <- Model2String(dxstring, parms = model.state.init.cond2)
dy <- Model2String(dystring, parms = model.state.init.cond2)

sigma = 0.1
timesteps = 1000
deltat = 0.1

# first simulate where x wins
init.cond.xwins <- c(x = 3, y = 5) # goes to x = 10, y = 0

timesimulation.xwins <- TSTraj (y0 = init.cond.xwins, time = timesteps, 
				deltat = deltat, x.rhs = dx, y.rhs = dy, sigma = sigma)
TSPlot(timesimulation.xwins, deltat = deltat)
TSPlot(timesimulation.xwins, deltat = model.deltat, dim = 2)
TSDensity(timesimulation.xwins, dim = 1)
TSDensity(timesimulation.xwins, dim = 2)

# then simulate where y wins
init.cond.ywins <- c(x = 5, y = 3) # goes to x = 0, y = 10
timesimulation.ywins <- TSTraj (y0 = init.cond.ywins, time = timesteps, 
				deltat = deltat, x.rhs = dx, y.rhs = dy, sigma = sigma)
TSPlot(timesimulation.ywins, deltat = deltat)

########################################################################
### Step 3: Local quasi-potential calculations #########################
bounds.x = c(-5,20)     # upper and lower limit of X
bounds.y = c(-5,20)     # upper and lower limit of Y
#if only use (-5,20), ywins does not overlap with xwins, leading to problems
step.number.x = 2000 #4000  # number of division between upper and lower limit
step.number.y = 2000 #4000  # number of division between upper and lower limit

# compute quasipotential around stable equilibrium where x wins
xinit = 10 
yinit = 0

init.cond.xwins <- QPotential(x.rhs = dx, x.start = xinit, x.bound = bounds.x, 
			x.num.steps = step.number.x, y.rhs = dy, y.start = yinit, 
			y.bound = bounds.y, y.num.steps = step.number.y)
			
QPContour(surface = init.cond.xwins, dens = c(1000, 1000), x.bound = bounds.x,
	y.bound = bounds.y, c.parm = 10)
dev.new()

# compute quasipotential around stable equilibrium where y wins
xinit = 0
yinit = 10

init.cond.ywins <- QPotential(x.rhs = dx, x.start = xinit, x.bound = bounds.x, 
			x.num.steps = step.number.x, y.rhs = dy, y.start = yinit, 
			y.bound = bounds.y, y.num.steps = step.number.y)
QPContour(surface = init.cond.ywins, dens = c(1000, 1000), x.bound = bounds.x,
	y.bound = bounds.y, c.parm = 10)

# What happens if you don't start at the equilibrium?  This:
#bounds.x = c(-5,20)     # upper and lower limit of X
#bounds.y = c(-5,20)     # upper and lower limit of Y
#xinit = 0
#yinit = 15
#init.cond.ywins.notateq <- QPotential(x.rhs = dx, x.start = xinit, x.bound = bounds.x, 
#			x.num.steps = step.number.x, y.rhs = dy, y.start = yinit, 
#			y.bound = bounds.y, y.num.steps = step.number.y)
#QPContour(surface = init.cond.ywins.notateq, dens = c(1000, 1000), x.bound = bounds.x,
#	y.bound = bounds.y, c.parm = 10)


########################################################################
### Step 4: Global quasi-potential calculation #########################
global.qp <- QPGlobal(local.surfaces = list(init.cond.xwins, init.cond.ywins),
		unstable.eq.x = c(2.323), unstable.eq.y = c(5.116), 
		x.bound = bounds.x, y.bound = bounds.y)
QPContour(surface = global.qp, dens = c(1000, 1000), x.bound = bounds.x,
	y.bound = bounds.y, c.parm = 10)
	
#if we want to know what the quasi-potential values at the stable eq are:
QPInterp(X = 0, Y = 10, x.bound = bounds.x, y.bound = bounds.y, surface = global.qp)
QPInterp(X = 10, Y = 0, x.bound = bounds.x, y.bound = bounds.y, surface = global.qp)


########################################################################
### Step 5: Global quasi-potential visualization #######################
# Use QPContour above or view in 3D, if you like that sort of stuff
library(plot3D)
library(viridis)
subsample.x = 200
subsample.y = 200

global.subsample <- global.qp[seq(1,step.number.x,step.number.x/subsample.x), seq(1,step.number.y,step.number.y/subsample.y) ]
global.qp2 <- global.subsample
maxvalue = 5
min(global.qp, na.rm = TRUE)
global.qp2[global.subsample > maxvalue] <- NA
# play with theta to see it from different points of view
persp3D(z = global.qp2, theta = 0, phi = 0, col = viridis(100, option = "C"), shade = 0.1, colkey = list(side = 4, length = 0.85), contour = T, zlim = c(-0.001, maxvalue), zlab = intToUtf8(0x03A6))

### Step 6: Vector field decomposition #################################

VDAll <- VecDecomAll(surface = global.qp, x.rhs = dx, 
				y.rhs = dy, x.bound = bounds.x, y.bound = bounds.y)
#should look familiar, deterministic skeleton
VecDecomPlot(x.field = VDAll[,,1], y.field = VDAll[,,2], dens = c(25, 25), 
			x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 20), 
			ylim = c(0, 20), arrow.type = "equal", tail.length = 0.25, 
			head.length = 0.025)

#gradient field
VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25, 25), 
	x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", 
	tail.length = 0.25, head.length = 0.025, xlim = c(0,20), 
	ylim = c(0,20))
VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25, 25), 
	x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", 
	tail.length = 0.25, head.length = 0.025, xlim = c(0,20), 
	ylim = c(0,20))

#remainder field
VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25, 25), 
	x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", 
	tail.length = 0.35, head.length = 0.025, xlim = c(0,20), 
	ylim = c(0,20))
VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25, 25), 
	x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", 
	tail.length = 0.35, head.length = 0.025, xlim = c(0,20), 
	ylim = c(0,20))








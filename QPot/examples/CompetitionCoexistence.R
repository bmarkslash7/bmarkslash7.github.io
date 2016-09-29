install.packages('QPot'); install.packages('plot3D')
install.packages('phaseR'); install.packages('deSolve')

########################################################################
########################################################################
### COEXISTENCE in Lotka-Volterra competition model
########################################################################
########################################################################

### First look at the situation where both competitors coexist
dxstring <- 'rx*x*(kx - x - alphaxy*y)*(1/kx)'
dystring <- 'ry*y*(ky - y - alphayx*x)*(1/ky)'

# these parameters lead to coexistence between two competitors
model.state.coexist <- c(rx = 4, ry = 4, kx = 10, ky = 10, 
							alphaxy = 0.5, alphayx = 0.5)

########################################################################
### Step 1: Analyze the deterministic skeleton #########################
# First we can work with model.state.coexist
# coexistence, stable eq at (6.667,6.667)
require(phaseR)
competition.model <- function(t, y, parameters){
	x <- y[1]							# initial state of species x
	y <- y[2]							# initial state of species y
	rx <- parameters["rx"]				# growth rate of species x
	ry <- parameters["ry"]				# growth rate of species y
	alphaxy <- parameters["alphaxy"]	# competition coefficient
	alphayx <- parameters["alphayx"]	# competition coefficient
	kx <- parameters["kx"]				# carrying capacity species x
	ky <- parameters["ky"]				# carrying capacity species y
	dx <- rx*x*(kx - x - alphaxy*y)*(1/kx)	# equation for species x
	dy <- ry*y*(ky - y - alphayx*x)*(1/ky)	# equation for species y
	list(c(dx,dy))
}

# this function will draw the vector field that describes how the 
# densities will change
flowField(deriv = competition.model, x.lim = c(0, 20), y.lim = c(0, 20), 
		parameters = model.state.coexist, add = F, points = 30)
	# x.lim and y.lim are ranges for x axis and y axis, respectively 
	# add = F makes new figure, points*points is the number of arrows 
	# include xaxs = 'i', yaxs = 'i' to get rid of white space in graph

# Draw the null clines, where dx/dt = 0 or dy/dt = 0
supp.print <- nullclines(deriv = competition.model, x.lim = c(0, 20), 
				y.lim = c(0, 20), parameters = model.state.coexist, 
				col = c("blue", "red"), points = 250)
	# species x is blue, species y is red
	# nullclines() estimates null cline, more points = better estimation 

#plot points on graph for stable equilbrium (black) and unstable (white) 
points(x = 0, y = 10 , pch = 21 , bg = "white" , cex = 1.5)
points(x = 10, y = 0 , pch = 21 , bg = "white" , cex = 1.5)
points(x = 6.667, y = 6.667 , pch = 21 , bg = "black" , cex = 1.5)
points(x = 0, y = 0 , pch = 21 , bg = "white" , cex = 1.5)

# draw in a trajectory at a random starting point
# can rerun multiple times to get more trajectories
# can change xstart and ystart to choose starting value
xstart = runif(1,0,20)	#1 random point between 0 and 20 from uniform dist 
ystart = runif(1,0,20)	# change this to a value to set start point
traj <- trajectory(competition.model,y0=c(xstart,ystart),
				parameters = model.state.coexist , t.end = 250)
	# to make line thicker, include lwd = 1.5, 
	# to make circle darker, include pch = 16,

########################################################################
### Step 2: Stochastic simulation ######################################
require(QPot) # loads QPot so that we can use the functions

dxstring <- 'rx*x*(kx - x - alphaxy*y)*(1/kx)'	# L-V competition
dystring <- 'ry*y*(ky - y - alphayx*x)*(1/ky)'	# from earlier

# Parameters included in step one, but in case you are skipping around
model.state.coexist <- c(rx = 4, ry = 4, kx = 10, ky = 10, 
							alphaxy = 0.5, alphayx = 0.5)

# QPot needs the equations with parameter values as strings
# Model2String puts the parameter values into our equations
dx <- Model2String(model = dxstring, parms = model.state.coexist)
dy <- Model2String(model = dystring, parms = model.state.coexist)

# We need to define some values for our simulation
init.cond <- c(x = 3, y = 10)	# initial state of simulation
sigma = 0.1		# level of stochasticity, sigma = 0 is deterministic model
timesteps = 1000				# number of time steps
deltat = 0.1					# size of time step

# Run a stochastic simulation, but TSTraj does not return a plot
timesimulation <- TSTraj(y0 = init.cond, time = timesteps, 
				deltat = deltat, x.rhs = dx, y.rhs = dy, sigma = sigma)

# plot our stochastic simulation
# plot time series and histogram
TSPlot(mat = timesimulation, deltat = deltat)

# plot the time series on a two-dimensional phase plane: dim = 2
TSPlot(mat = timesimulation, deltat = model.deltat, dim = 2) 

# plot only the histogram part of the first TSPlot 
TSDensity(mat = timesimulation, dim = 1)

# plot a histogram on a two-dimensional phase plane
TSDensity(mat = timesimulation, dim = 2)
	# can set xlim = c(0,10), ylim = c(0,10) to plot a larger area

########################################################################
### Step 3: Local quasi-potential calculations #########################
# First we need to set some limits
bounds.x = c(-5,15)     # upper and lower limit of x
bounds.y = c(-5,15)     # upper and lower limit of y

step.number.x = 2000    # number of divisions in bounds.x
step.number.y = 2000    # number of divisions in bounds.y

xinit = 6.667			# equilibrium x value
yinit = 6.667			# equilibrium y value

# Compute the quasi-potential using the ordered upwind method
coexist <- QPotential(x.rhs = dx, x.start = xinit, x.bound = bounds.x, 
			x.num.steps = step.number.x, y.rhs = dy, y.start = yinit, 
			y.bound = bounds.y, y.num.steps = step.number.y) 

# Use QPContour to look at the local quasi-potential 
QPContour(surface = coexist, dens = c(1000, 1000), x.bound = bounds.x,
	y.bound = bounds.y, c.parm = 0) # add countour lines with c.parm = 5

# in default case, when QPotential hits bounds.x or bounds.y, it stops
# can use bounce = 'b', which prevents QPotential from doing this
coexist.BOUNCE <- QPotential(x.rhs = dx, x.start = xinit, 
		x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = dy, 
		y.start = yinit, y.bound = bounds.y, y.num.steps = step.number.y, 
		bounce = 'b') 

# Using bounce = 'b' leads to a yellow visualization
# because some of the calculated quasi-potentials are too large
QPContour(surface = coexist.BOUNCE, dens = c(1000, 1000), x.bound = bounds.x,
	y.bound = bounds.y, c.parm = 10)

# To fix this, we are going to remove all quasi-potentials > upper limit
hist(coexist.BOUNCE) # plots a histogram so we can decide on upper limit
		#from the histogram and graphing other values, 25 looks good. 

# need to remove all values greater than 25
# could use coexist.BOUNCE, but if we mess up, we have to rerun QPotential
coexist.BOUNCE.REMOVED <- coexist.BOUNCE
coexist.BOUNCE.REMOVED[coexist.BOUNCE.REMOVED > 25] <- NA # remove funk 

# Plot the quasi-potential only where the quasi-potential < 100
QPContour(surface = coexist.BOUNCE.REMOVED, dens = c(1000, 1000), 
				x.bound = bounds.x, y.bound = bounds.y, c.parm = 10)
# dens = c(1000,1000) plots a subset of points
# c.parm draws contour lines

########################################################################
### Step 4: Global quasi-potential calculation #########################
# For this example, the local quasi-potential is global quasi-potential
# Step 4 Done!

########################################################################
### Step 5: Global quasi-potential visualization #######################
# Use QPContour above or view in 3D, if you like that sort of stuff

require(plot3D)
subsample.x = 200 # 3D plotting all the points takes forever
subsample.y = 200 # so we subsample and only plot 200*200 points

my.theta = 110				# rotation on 3D graph
my.phi = 0					# altitude of 3D graph
my.zlim = c(-0.001, 25) 	# range of z axis on plot, by making many 
							# graphs, we decided 25 looks good
							# also notice that 25 removed the funk
my.zlab = intToUtf8(0x03A6) # this allows us to plot the letter phi

# First, look at the quasi-potential with default algorithm behavior
# subsample so that plotting is fast
coexist.subsample <- coexist[
						seq(1,step.number.x,step.number.x/subsample.x), 
						seq(1,step.number.y,step.number.y/subsample.y) ]

# plot the 3D graph of the quasi-potential with default default behavior
persp3D(z = coexist.subsample, theta = my.theta, phi = my.phi, 
		contour = T, zlim = c(-0.001, 15), zlab = my.zlab)

# Now plot the quasipotential when we used bounce = 'b'
dev.new() # this opens a new window so that we can compare plots

# again, need to subsample from the bigger matrix to make graphing faster
coexist.BOUNCE.SUBSAMPLE <- coexist.BOUNCE.REMOVED[
							seq(1,step.number.x,step.number.x/subsample.x), 
							seq(1,step.number.y,step.number.y/subsample.y) ]

# plot the quasipotential when we used bounce = 'b'
persp3D(z = coexist.BOUNCE.SUBSAMPLE, theta = my.theta, phi = my.phi, 
		contour = T, zlim = my.zlim, zlab = my.zlab)

########################################################################
### Step 6: Vector field decomposition #################################
# First get all the vector fields that you will need with VecDecomAll
VDAll <- VecDecomAll(surface = coexist.BOUNCE, x.rhs = dx, 
				y.rhs = dy, x.bound = bounds.x, y.bound = bounds.y)

# now set some values that will make plotting easier
x.skel = VDAll[,,1]		# deterministic skeleton of x equation
y.skel = VDAll[,,2]		# deterministic skeleton of y equation
x.grad = VDAll[,,3]		# gradient field of x  
y.grad = VDAll[,,4]		# gradient field of y 
x.redr = VDAll[,,5]		# remainder of x
y.redr = VDAll[,,6]		# remainder of y

# Now plot the vector fields
# deterministic skeleton, should look familiar from step 1
VecDecomPlot(x.field = x.skel, y.field = y.skel, dens = c(25, 25), 
		x.bound = bounds.x, y.bound = bounds.y,  arrow.type = "equal", 
		xlim = c(0, 15), ylim = c(0, 15), 
		tail.length = 0.25, head.length = 0.025)

# gradient field: first plot with arrow size proportional to gradient
VecDecomPlot(x.field = x.grad, y.field = y.grad, dens = c(25, 25), 
		x.bound = bounds.x, y.bound = bounds.y, 
		arrow.type = "proportional", tail.length = 0.25, 
		head.length = 0.025, xlim = c(0,15), ylim = c(0,15))

# gradient: plot with all arrows the same size, independent of gradient
VecDecomPlot(x.field = x.grad, y.field = y.grad, dens = c(25, 25), 
		x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", 
		tail.length = 0.25, head.length = 0.025, xlim = c(0,15),
		ylim = c(0,15))

#remainder field: plot with arrow size proportional to remainder
VecDecomPlot(x.field = x.redr, y.field = y.redr, dens = c(25, 25), 
	x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", 
	tail.length = 0.35, head.length = 0.025)

#remainder field: plot with arrow size independent of remainder
VecDecomPlot(x.field = x.redr, y.field = y.redr, dens = c(25, 25), 
	x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", 
	tail.length = 0.35, head.length = 0.025, xlim = c(0,15), 
	ylim = c(0,15))

# HARMPI Tutorial by Sasha Tchekhovskoy

##How to set up HARMPI: choose the problem, compile, and run the code:

* To install the code, you do:

		git clone git@github.com:atchekho/harmpi.git
		cd harmpi
		make clean
		make

	If you would like to use ICC instead of GCC, then edit makefile and replace "USEICC = 0" with "USEICC = 1".

	Before running one of the problems with HARMPI, you choose the problem you want to run in the file decs.h.

	You do so by modifying decs.h. There are 8 problems to choose from:

		//which problem
		#define MONOPOLE_PROBLEM_1D 1
		#define MONOPOLE_PROBLEM_2D 2
		#define BZ_MONOPOLE_2D 3
		#define TORUS_PROBLEM 4
		#define BONDI_PROBLEM_1D 5
		#define BONDI_PROBLEM_2D 6
		#define SNDWAVE_TEST 7
		#define ENTWAVE_TEST 8

	The default choice is

		#define WHICHPROBLEM TORUS_PROBLEM

	This is a magnetized torus accretion problem in 2D. For starters, however, it might be easier to start with a 1D problem without magnetic field, e.g., with `BONDI_PROBLEM_1D`.

* Compile:

		make clean
		make

* Choosing desired resolution

	You can choose the desired resolution per MPI *tile* by changing the values of `N1`, `N2`, and `N3` for the chosen problem in [decs.h](decs.h).

* To run the code

	You choose the number of tiles along each of the 3 directions at run time. E.g., if `N1 = N2 = 64` and `N3 = 1`, and you choose to run with two cores in r- and theta-directions and one core in the phi-direction, you get a resolution of 128x128x1 cells by running

		mpirun -n 4 ./harm 2 2 1

* To run the code in serial (on a single core), do the following:

		./harm 1 1 1

	For the example above, this would give a total resolution of 64x64x1 cells.

	As the code runs, it produces sequential dump files in the 'dumps' subdirectory, i.e., dumps/dump###. It also produces dumps/gdump file that contains the information about the metric and the grid.

## How to read in the output into Python

You can use any visualization program you are used to (Python, Gnuplot, IDL, Matlab, etc)
* The scripts for reading in the data in Python are located in [harm_script.py](harm_script.py)
* For a great Python tutorial, check out <http://www.learnpython.org/>

## To use these Python scripts, do:

* Start Interactive Python Shell:

		ipython --pylab

	You should do this from the location of the harm executable.

* Initialize the scripts:

		%run -i harm_script.py

* Read in grid information:

		rg("gdump")

* Read in dump ﬁle:

		rd("dump000")

	where you can replace "000" with any other number of a dump file (provided it was already written out).

* Once you read in a dump file, you can print its current time via

		print t

## Understanding the grid

HARMPI (like HARM2D) is capable of using non-uniform grids. The examples below will automatically take care of converting from the grid coordinates, which are x1, x2, x3, to the physical coordinates, which are r, theta, phi.

## How to visualize the output in Python
* Python Visualization Examples

	For a 1D test: scroll down to the bottom of harm_script.py and replace

		if False:
			#1D plot example

	with

		if True:
			#1D plot example

	For a 2D test: scroll down to the bottom of harm_script.py and replace

		if False:
			#2D plot example

	with

		if True:
			#2D plot example

* Plot the 1D radial dependence of density along the equatorial plane

		plt.loglog(r[:,ny/2,0], rho[:,ny/2,0])
		plt.xlabel("r")
		plt.ylabel("rho")

* Plot 2D (R-z) contours of logarithm of density in color (with color bar):

		plco(np.log10(rho),xy=1,xmax=100,ymax=50,cb=True)

	Here `xmax = 100` and `ymax = 50` set the plotting range. You can omit these arguments, and a default range will be used.

* Compute magnetic ﬂux function:

		psi = psicalc()

* Overplot magnetic ﬁeld lines in black:

		plc(psi,xy=1,colors="black")

* Compute Lorentz factor:

		alpha = (-guu[0,0])**(-0.5) #lapse
		gamma = alpha*uu[0]

* Compute energy output of a black hole:

		aux() #compute OmegaF and T^\mu_\nu
		dEr = -gdet*Tud[1,0]*_dx2*_dx3
		Er = dEr.sum(2).sum(1) #sum in phi (axis=2) and theta (axis=1)

* Plot total radial energy ﬂux vs. radius:

		x=r[:,0,0] #x is now the 1D radial coordinate
		plt.plot(x,Er) #make a 1D plot of energy ﬂux vs. x
		plt.xscale("log") #logarithmic x-axis
		plt.xlim(rhor,20) #limit x-range to 20r_g

* Plot the angular frequency of field line rotation vs. theta on the black hole horizon:

		OmegaF = omegaf2 #the value of omegaf2 is computed inside aux()
		OmegaH = a / (2*rhor)
		rhor = 1+(1-a*a)**0.5
		ihor = ti[r>=rhor][0] #compute the index of the cell at the location of the horizon
		plt.plot(h[ihor,:,0],OmegaF[ihor,:,0]/OmegaH)
		plt.ylabel("OmegaF/OmegaH")
		plt.xlabel("theta")
		plt.ylim(0,1)
		plt.xlim(0,np.pi)

* Show ergosphere with thick red line:

		rergo = 1+np.sqrt(1-(a*np.cos(h))**2)
		plc(r-rergo,levels=(0,),xy=1,colors="red",linewidths=2)

* Show black hole with black circle:

		#generate ellipse object
		#import the Ellipse function
		from matplotlib.patches import Ellipse
		#create Ellipse "actor"
		el = Ellipse(
			(0,0), #origin of ellipse
			2*rhor, 2*rhor, #major axes of ellipse
			facecolor='k', #color of ellipse
			alpha=1) #opacity of ellipse
		ax=plt.gca() #get current axes
		art=ax.add_artist(el) #show object in axes
		art.set_zorder(20) #bring ellipse to front
		plt.draw() #make sure it's drawn

* Save ﬁgure ﬁle

		plt.savefig("figure.pdf")

* Clear the current window

		plt.clf()

* Open a new window

		plt.figure()

	

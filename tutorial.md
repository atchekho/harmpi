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

	This is a magnetized torus accretion problem in 2D. For starters,
    however, it might be easier to start with a 1D problem without
    magnetic field, e.g., with `BONDI_PROBLEM_1D`.

* Brief descriptions of the problems

	* Monopole problems (1D and 2D versions; `MONOPOLE_PROBLEM_1D` and
      `MONOPOLE_PROBLEM_2D`)

		Use these problems for studying *acceleration* and black hole *power
		output* in an unconfined magnetosphere: radial magnetic field
		lines, low-density atmosphere, outer radius `Rout = 1e3`
		gravitational radii. You can increase the size of the domain to
		study acceleration over longer range in distance. To increase
		`Rout` from the default value of `1e3` to e.g. `1e4`, go to
		[init.c:74](init.c#L74) and change `init_monopole(1e3);` to
		`init_monopole(1e4);`. Note that in order to keep the effective
		resolution the same, you'd need to increase the resolution by a
		factor of `log(1e4)/log(1e3) = 4/3`. There are two versions of the
		problem:

		* `MONOPOLE_PROBLEM_1D`: 1D version of the problem,
		with the grid in the midplane (`z = 0`). In this setup, the
		magnetic poloidal field lines are forced to be pointing radially outward
		and are pinned (e.g., cannot move around).

		* `MONOPOLE_PROBLEM_2D`: 2D version of the problem. Here, the
		poloidal magnetic field lines are free to move around in the R-z
		plane.

	* BZ monopole problem (`BZ_MONOPOLE_2D`)

		Same as `MONOPOLE_PROBLEM_2D` but with a smaller outer
		radius, `Rout = 1e2` gravitational radii. Due to the limited
		radial range, this problem would not allow you to study the
		acceleration of the outflow, for which >= 3 orders of magnitudes
		in radius is ideal. However, this problem is great for studying
		the power output of the black hole.
    
	* Torus problem (`TORUS_PROBLEM`)

		Consists of an equilibrium hydrodynamic torus on
		an orbit around the black hole. The torus has an inner radius of
		`rin = 6` gravitational radii, pressure maximum at `rmax =
		13`. This leads to torus' outer extent of about `50` gravitational
		radii. Without any magnetic field, the torus would sit on an orbit
		forever. To get some action, we insert a *poloidal* (in R-z plane)
		magnetic field loop into the torus. This magnetic field is
		unstable to the magnetorotational instability (MRI), which
		transports the angular momentum outward and causes the gas and
		magnetic field to move inward, toward the black hole. This
		culminates in black hole receiving large-scale magnetic flux,
		which we originally inserted into the disk, and this flux powers a
		twin relativistic jets. You can make a movie of this by following
		instructions below.

	* Bondi problem (1D and 2D; `BONDI_PROBLEM_1D`, `BONDI_PROBLEM_2D`)

		These two problems start with a uniform density and temperature
		at `r >= 10` gravitational radii. Inside of this "hole", I put
		very low density and pressure. As you run the simulation, the
		"hole" gets quickly eaten by the black hole. What results is a
		spherically-symmetric accretion from a uniform-density ambient
		medium.

	* Sound wave problem (`SNDWAVE_TEST`)

	   This is a sound wave test


    * Entropy wave problem (`ENTWAVE_TEST`)

	   This is an entropy wave test

* Compile:

		make clean
		make

* Resolving possible compilation problems

	* Sometimes `-lm` switch needs to be added (in
    [makefile:74](makefile#L74))

* Choosing desired resolution

	You can choose the desired resolution per *tile*, or per MPI process, by changing the values of `N1`, `N2`, and `N3` for the chosen problem in (decs.h)[decs.h]. If there are more than one tile, this will be smaller than the total resolution.

* Choosing the number of tiles and running the code

	You choose the number of tiles along each of the 3 directions at run time. E.g., if `N1 = N2 = 64` and `N3 = 1`, and you choose to run with 4 cores in the r-direction, 2 cores in the theta-direction and one core in the phi-direction, you get a total resolution of 256x128x1 cells by running

		mpirun -n 8 ./harm 4 2 1

	Here, `-n 8` option tells `mpirun` to use `8` (`= 4*2*1`) total tiles (one
    tile per core), `./harm` specifies the executable ("`./`" is to
    indicate the current directory) and the arguments `4 2 1` specify
    how these tiles are distributed among the 3 directions. Note: OpenMP
    support is presently in development, but it is not yet ready for
    prime time.

* To run the code in serial (on a single core), do the following:

		./harm 1 1 1

	For the example above, this would give a total resolution of 64x64x1 cells. Note: it is a good idea to clean out the `dumps/` directory before you run the code, or else the code would attempt to restart from `dumps/rdump0`.


	As the code runs, it produces sequential dump files in the 'dumps' subdirectory, i.e., dumps/dump###. It also produces dumps/gdump file that contains the information about the metric and the grid.

* Resolving potential run-time problems

	* If you experience `Signal code: Integer divide-by-zero (7)`
    error, you have two options
		* upgrade to 
		[OpenMPI v. 2.0](https://www.open-mpi.org/software/ompi/v2.0/ "OpenMPI
		v. 2.0"). Combination of latest versions of 
		[GCC v. 6.1](https://gcc.gnu.org/gcc-6/ "GCC v. 6.1") and [OpenMPI v. 2.0](https://www.open-mpi.org/software/ompi/v2.0/ "OpenMPI
		v. 2.0") has been
		verified to work. (thanks Matthias Raives for the tip)
		* follow the below solution for `File locking failed` problem

	* If you experience `File locking failed`, try setting `#define
      DO_PARALLEL_WRITE (0)` in [decs.h:369](decs.h#L369). This will force each
      MPI process to write its own file instead combining the output
      into a single file. In this case, the outputs of different
      processes will be combined upon reading in the files in
      python. Note: once you do this, each dump file will no longer be
      a single file but multiple files, one file per core. So, if you
      are running on 4 cores, instead of `dump000` you would get 4
      files, `dump000_0000`, `dump000_0001`, `dump000_0002`,
      `dump000_0003`. You can read in these files the usual way,
      `rd("dump000")`: the python macros will take care of merging the
      files together automatically.


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

* Make a movie:

	You can make a simple movie by running

		mkmov_simple(starti = 0,endi = 100)

	Here `starti` is the first frame number and `endi` is the last frame number. Note: the color limits are tuned for the torus problem (`TORUS_PROBLEM`), and you'd have to adjust them for other problems.

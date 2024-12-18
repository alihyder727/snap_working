# Simulating Non-hydrostatic Atmospheres on Planets (SNAP)
<!-- Jenkins Status Badge in Markdown (with view), unprotected, flat style -->
<!-- In general, need to be on Princeton VPN, logged into Princeton CAS, with ViewStatus access to Jenkins instance to click on unprotected Build Status Badge, but server is configured to whitelist GitHub -->
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
<!--[![Build Status](https://travis-ci.com/luminoctum/athena19-dev.svg?token=AfxC7sH2UkyrrtpsBrob&branch=dev)](https://travis-ci.com/luminoctum/athena19-dev) -->
[![Build Status](https://github.com/luminoctum/athena19-dev/actions/workflows/autotest.yml/badge.svg)](https://github.com/luminoctum/athena19-dev/actions/workflows/autotest.yml)

<!--[![Public GitHub  issues](https://img.shields.io/github/issues/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/issues)
[![Public GitHub pull requests](https://img.shields.io/github/issues-pr/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/pulls) -->

Athena++ Atmospheric Dynamics Code

## Coding Style Guide 
### Functions
1. If a function returns more than one elementary variables, pass pointers as arguments.
1. If a function returns complex variables, pass references.
1. The returned variables are normally placed in the first few arguments in a function
   - Exception #1: In MPI calls, the send argument is always placed in the first.
1. Function names should begin with a lowercase verb indicating the action and dependent
   words using initial capitals.
   - Exception #1: conversion function can be named with either capitalized names or
     uncaptialized names to preserve symmetry. An example is *deg2rad*. It is allowed to
     use either *2* or *To* in the name of the conversion function.
   - Exception #2: if a function name is a single verb, the verb should be
     capitalized. An example is *Vectorize*.
   - Exception #3: C function names can be lowercase verb phrases concatenated by *underscores*.

### Classes and structures
1. Class/structure names should be contatenated noun phrases capitalizating the first
   letter.
1. Class/structure names should be in singular form.
1. Class/structure names should not have *underscores* in the name.

### Variables
1. Variable names should begin with a lower case letter
1. Variable names can contain *underscores*.
1. Private or protected variables must contain an *underscore* at the end of the name.
1. Variable names for pointers should start with the letter *p*.
    - Exception #1: previous or next pointer in a linked list can be named as *prev* or
      *next.

### One dimensional arrays
1. One dimensional arrays should be a C++ vector container.
1. An one dimensional arrays is passed to a function by the pointer to its element and
   the size of the array.

### Multidimensional arrays
1. Spatially dependent multidimensional arrays are allocated as the AthenaArray
1. Other multidimensional arrays are allocated by template functions.

### Choice of units
1. The units of a physical variable should use SI by default. Using units different from SI 
  should be specified by appending *underscore units* at the end of the variable name.
   - Example #1: pres\_bar means *pressure in bars*.

### Physical constants
1. Physical constants are declared as *const Real const* under a Class or under a
   Namespace.

### Pointers
1. Use *nullptr* for null pointer.

### Indentation
1. Use 2 spaces for indentation.

### Documentation
1. Block documentation
1. Function documentation
1. Variable documentation
1. Special commands
    - \todo
    - \bug
    - \warning
    - \deprecated

## General build and compilation:

- In main directory, snap_working/, first patch the code:
	- ./patch.py

- Configure the code for the specific problem and serial/parallel run. This example is for the moist Jupiter runs with PH3 and GeH4 scalars. Setting nscalars=1 will configure the problem for CO, exclusively.
	- ./configure.py --prob=convection --flux=lmars --nvapor=2 --nscalars=2 --nghost=3 -netcdf

- Compile using "make". You can also use "make -jx", where x should be changed to however many cores you are trying to compile over.

- The executable will be created in the bin/ folder. The main executable needed will have the same name as the --prob variable used to set up the configuration. In this case, it would be convection.exe. Make a run directory where the output will be stored. Ensure that the input file ".inp" is located within this run directory.	- mkdir example_run
	- cd example_run/
	- cp ../bin/convection.exe .
	- cp ../moist_jupiter/jupiter2d.inp/.

- Execute the code (assuming serial run)
	- ./convection.ex -i jupiter2d_moist.inp > log.jupiter_moist &

## Output

- The run outputs .nc files, for which you need netCDF installed. To view these, softwares like ncview can be used. There will be many .nc files so first we must combine them using the functionality provided. We will name the final file "test" for this example:
	- ../combine.py -o test
	- ncdump -h jupiter-test-main.nc

- The final command will provide a netcdf output with self-descriptive units and long names. The final .nc file can be viewed as a movie using the ncview functionality of the netCDF library. However, custom scripts can be written as well. If using Python, you must ensure that the netCDF4 library is installed.

NOTE: The softlinks used for the main code will be appropriately linked upon the execution of ./patch.py.

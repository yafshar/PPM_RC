PPM\_RC program
===============
v1.0, Feb 2016

Introduction
------------

Purpose : 2D-3D distributed parallel region competitioni [1] segmentation algorithm 

Images can be so large that they do not fit the main memory of a single computer. 
We address this issue by developing a distributed parallel algorithm for segmentation 
of large fluorescence microscopy images. 
The algorithm is based on the versatile Region Competition [1] method, which has 
previously proven useful in microscopy image analysis. 
The present distributed implementation decomposes the input image into smaller sub-images 
that are distributed across multiple computers. Using network communication, the computers 
orchestrate to collectively solve the global segmentation problem.

1. J. Cardinale, G. Paul, and I. F. Sbalzarini. Discrete region competition for unknown 
   numbers of connected regions. IEEE Trans. Image Process., 21(8):3531–3545, 2012


Initial release of PPM\_RC source program.

Getting the code
------------
You can download the latest version from [here](https://github.com/yafshar/PPM_RC). 
The very latest version is always available via 'github' by invoking one of the following:
````
## For the traditional ssh-based Git interaction:
$ git clone git@github.com:yafshar/PPM_RC.git

## For HTTP-based Git interaction
$ git clone https://github.com/yafshar/PPM_RC.git
````

Package contents
------------

The package you checked out by GIT access should contain on the toplevel the
following files and directories:

~~~~~~~~~~~~~~~~~~~~~~~~~~~
-    34K COPYING       The License
-   7.5K COPYING.LESSER
-    38B Makefile.am
-   4.2K Makefile.in   The Makefile template used by configure
-   653B NOTICE
-    14K README        This file
-    28K README.md   
-    46K autogen.sh
-    42K config.guess
-    35K config.sub
-   524K configure     The build configure script
-    11K configure.ac  Thie build configure script template
-    14K install-sh
-   405B settings.in
-   782B m4            Contains build scripts
-   1.1K src           Contains the source code
-   204B utils         Contains the PPM\_RC debug utility
~~~~~~~~~~~~~~~~~~~~~~~~~~~

After you compile PPM\_RC there will be further directories for the binaries
and include files.

Requirements for building PPM\_RC
----------------------------------
- Fortran compiler with Fortran 2003 support
- C compiler
- C++ compiler

- METIS 5: You may download the latest release of METIS 5 from
  http://glaros.dtc.umn.edu/gkhome/metis/metis/download.

- An MPI distribution (recommended): Either get OpenMPI, mpich2 or any other MPI-3
  compliant MPI library. If you are compiling PPM\_RC on a cluster, most likely your
  sysadmin will have already an MPI installed on the system.

- PPML & PPM core: 
  The PPM\_RC is based on PPML and PPM core 1.2.2 latest development version 
  (object-oriented Fortran 2003) and can only be linked against this version. 
  Please first check out the latest PPM development version (object-oriented 
  Fortran 2003) by anonymous Git access: 
  For the PPML language and compiler:  `git clone http://ppm.mpi-cbg.de/git/cg.git` 
  For the PPM library core:            `git clone http://ppm.mpi-cbg.de/git/ppm.git`.
  For more information and Installation Manual Please refer to: 
  http://mosaic.mpi-cbg.de/?q=downloads/ppm\_lib
  Compile PPM core 1.2.2 latest development version before attempting to 
  compiling this package.
  Make sure that all requirements are compiled with the same compiler that you
  will be using to build PPM core.

- TIFF Library:  
  For handling huge images, or very large collections of images, breaking the 
  4 gigabytes boundary, you need a TIFF Library which supports BigTIFF.

- BOOST C++ Library: 
  For more information Please refer to:
  http://www.boost.org/


Building PPM\_RC 
-----------------

PPM\_RC is built in 3 simple steps:

* Step 1: Confguring PPM\_RC

Run the `configure` script to allow the build system to determine the correct
options to compile PPM\_RC.

It is very important to give `configure` the correct settings to make sure PPM\_RC
is compiled correctly. To find out which settings are supported type

`./configure --help`

This is what will be returned:

`configure` configures PPM\_RC to adapt to many kinds of systems.

Usage: `./configure [OPTION]... [VAR=VALUE]...`

To assign environment variables (e.g., `CC, CFLAGS...`), specify them as
`VAR=VALUE`.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

~~~~~~~~~~
Configuration:
  -h, --help              display this help and exit
  -V, --version           display version information and exit
~~~~~~~~~~

By default, `make` will create an executable file in `$ppm_rc_root/run`.  

For better control, use the options below.
Optional Features:
~~~~~~~~~~
  --enable-mpi            use MPI (default is no), If the MPI implementation
                          of your choice provides compile wrappers that are in
                          PATH, I can set them myself, choose: guess (I will
                          choose the first implementation I can find). Else,
                          set CC, CXX and FC to the appropriate compiler
                          wrappers (safest)
  --enable-debug          enable debug data generation (default is no)
  --enable-linux          compile for linux (default is no)
~~~~~~~~~~

Optional Packages:
~~~~~~~~~~
  --with-ppm=path         set the path to the ppm core library - THIS FLAG IS MANDATORY
  --with-metis=path       user defined path to METIS library
  --with-boost=path       Specify the root directory for boost library
  --with-tiff=path        user defined path to TIFF library
~~~~~~~~~~

Some influential environment variables:
~~~~~~~~~~
  FC          Fortran compiler command
  FCFLAGS     Fortran compiler flags
  LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
              nonstandard directory <lib dir>
  LIBS        libraries to pass to the linker, e.g. -l<library>
  CC          C compiler command
  CFLAGS      C compiler flags
  CPPFLAGS    (Objective) C/C++ preprocessor flags, e.g. -I<include dir> if
              you have headers in a nonstandard directory <include dir>
  CXX         C++ compiler command
  CXXFLAGS    C++ compiler flags
  CPP         C preprocessor
  CLIBS       libraries to pass to the linker for C compiler, e.g. -l<library>
~~~~~~~~~~

Use these variables to override the choices made by `configure' or to help
it to find libraries and programs with nonstandard names/locations.

Report bugs to the package provider.

Following options are especially important:

- `--enable-mpi`: If you will be running PPM_RC on a parallel environment
  (a cluster) using MPI. 
  If your system is properly configured then this should be enough
  information for PPM_RC build system to find the MPI libraries and compiler
  wrappers needed. If this goes wrong, you may ommit this option and set
  compiler wrapper and libraries in `FC` and `LDFLAGS` respectively.
- `--enable-linux`: Set this if you're compiling/running on a Linux system
- `--prefix`: If you like to install PPM_RC and the target directory is not the
  system's standard directory (`/usr/`) then you have to define this directory
  here. You must provide the full path. 
  It is not necessary to install PPM_RC.
  Building it and leaving it in the compilation directory is sufficient. If you
  provide a directory here it must already exist - it will not be created by the
  build system.
- `FC` etc.: If you wish to not use MPI or you have to specify exactly which
  compiler executable should be used, then you can use this flag to set your
  compiler.
- `LDFLAGS`: If metis was not installed in one of the system's standard library
  directories (e.g. `/usr/lib`) you must specify the directory to the libmetis.a
  file here.

Here two examples on how you could run the configure command

`.configure` on Linux cluster using OpenMPI (and intel compilers, wrapped)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
./configure CPP=cpp --enable-mpi --enable-linux --with-ppm=/home/usr/ppmcore 
--with-metis=/home/usr/metis --with-boost=/home/usr/boost_1_58_2 --with-tiff=/home/usr/tiff4
~~~~~~~~~~~~~~~~~~~~~~~~~~~

`./configure` on Mac OS X workstation with HomeBrew compilers
~~~~~~~~~~~~~~~~~~~~~~~~~~~
./configure CPP=cpp-5 --enable-mpi --with-ppm=/home/usr/ppmcore 
--with-metis=/home/usr/metis --with-boost=/home/usr/boost_1_58_2 --with-tiff=/home/usr/tiff4
~~~~~~~~~~~~~~~~~~~~~~~~~~~

`./configure` on a computer with OpenMPI installed in a non-standard location
~~~~~~~~~~~~~~~~~~~~~~~~~~~
./configure --enable-mpi FC=/opt/openmpi/1.10.2/bin/mpif90 CC=/opt/openmpi/1.10.2/bin/mpicc CXX=/opt/openmpi/1.10.2/bin/mpic++ 
--with-ppm=/home/usr/ppmcore --with-metis=/home/usr/metis --with-boost=/home/usr/boost_1_58_2 --with-tiff=/home/usr/tiff4
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Step 2: Compiling PPM_RC

If the configure process finished successfully you should see on your screen a
message that the Makefile has been generated (and you can now find this
Makefile in this directory).

Now you can simply run make to compile PPM_RC:

`make`

If you encounter problems in the compilation process (compile errors) please,
first check if you have set everything correctly in your environment. If the
error persists, please send us a bug-report detailing the previous steps you
have performed. Also, please include the `config.log` file and the output of
`export`. Finally, if yu are using MPI, please include which MPI library you are
using.


* Step 3: Installing PPM_RC (optional)

If you wish to install PPM_RC you can now use the `make install` command to do
so:

`make install`

If the target directory is part of the system, you will most probably get a
message that you have insufficient rights. If you have a root account you can 
use in this case the sudo command to override this security setting.

`sudo make install`

Your PPM_RC is installed.

Enjoy the PPM_RC experience!


Contributors
------------
PPM_RC package maintainer: Yaser Afshar <afshar@mpi-cbg.de>


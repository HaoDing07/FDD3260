# Unit Testing for MIMICA v5 with pFUnit

pFUnit is a unit testing framework for fortran. <br>
https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git

For more background information on unit testing (and software testing in general) <br>
https://coderefinery.github.io/2021-03-17-testing-hackathon/

## Installation
### Requirements
see also https://github.com/Goddard-Fortran-Ecosystem/pFUnit#prerequisites

Most importaintly ...
* cmake > 3.12
* ifort > 19.0.3

To install pFUnit 4.2 on NSC's Tetralith the following packages can be used (May/2021):

``` bash 
module load CMake/3.16.5
module load buildtool-easybuild/4.3.3-nscf4a9479 
module load intel/2020b
```


``` bash 
export FC=/software/sse/manual/mpprun/4.1.3/nsc-wrappers/ifort
```

### Installing pFUnit
```
# Clone the pFUnit repository
git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git

cd pFUnit

# Configure with Cmake
mkdir build
cd build
cmake -DMPI=YES -DOPENMP=YES -DCMAKE_INSTALL_PREFIX=./installed ../

# Build and install to local directory 'installed'
make -j tests
make install
```
Define the PFUNIT_DIR path (replace {XXX} with your pFUnit version)
```
export PFUNIT_DIR=<YOURWORKINGDIR>/build/installed/PFUNIT-{XXX}
```
#### Troubleshooting

In case of failures during the installation that give an error message like this:
```
You need to run this command from the toplevel of the working tree.
CMake Error at cmake/build_submodule.cmake:29 (message):
git submodule update --init failed with 1, please checkout submodules
```
(the important part is ``You need to run this command from the toplevel of the working tree.``)
you need to replace ``CMAKE_CURRENT_SOURCE_DIR`` with ``PROJECT_SOURCE_DIR`` in line 25 in ``cmake/build_submodule.cmake``. This can also occur in e.g. ``extern/fArgParse/extern/gFTL-shared/cmake/build_submodule.cmake``, line 29.

### Testing the install

The installation of pFunit can be testes by evaluating the pFUnit_demo cases

```
# Clone the pFUnit_demos repository
git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit_demos.git
```

Note: one has to add/replace the lines in the Makefile of the respective demo cases <br>
Replace {XXX} with your pFUnit version

```
include $(PFUNIT_DIR)/PFUNIT-{XXX}/include/PFUNIT.mk 
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
```


Try out the Trivial, Basic, and MPI examples
```
cd pFUnit_demos

cd Trivial
./build_with_cmake_and_run.x
cd ../Basic
./build_with_cmake_and_run.x
cd ../MPI
./build_with_cmake_and_run.x
```
#### Troubleshooting

In case of failures during the installation that give an error message like this:
```
CMake Error at CMakeLists.txt:7 (find_package):
  By not providing "FindPFUNIT.cmake" in CMAKE_MODULE_PATH this project has
  asked CMake to find a package configuration file provided by "PFUNIT", but
  CMake did not find one.

  Could not find a package configuration file provided by "PFUNIT" with any
  of the following names:

    PFUNITConfig.cmake
    pfunit-config.cmake

  Add the installation prefix of "PFUNIT" to CMAKE_PREFIX_PATH or set
  "PFUNIT_DIR" to a directory containing one of the above files.  If "PFUNIT"
  provides a separate development package or SDK, be sure it has been
  installed.


-- Configuring incomplete, errors occurred!
```
add the pFUnit cmake config to ``$CMAKE_PREFIX_PATH`` like this: 
```
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:<YOURWORKINGDIR>/build/installed/PFUNIT-{XXX}/cmake
```
## Usage

Once pFUnit is installed it can be used to test individual parts of MIMICA. The available tests for MIMICA will be collected and stored in mimicav5/pFUnitTests

If  _UNITTEST_ flag in ***start*** is set to _TRUE_ the testing routines will automaticall be built when MIMICA is compiled. Once MIMICA is compiled, the executables are available in ```./pFUnitTests/``` .

Have a look at the pFUnit demos  https://github.com/Goddard-Fortran-Ecosystem/pFUnit_demos.git to lean how to write unit tests in pFUnit

Principle currently steps are: 
- write test file *.pf
- add new test to ./pFUnitTests/Makefile
- compile MIMICA
- run executeables in ./pFUnitTests/
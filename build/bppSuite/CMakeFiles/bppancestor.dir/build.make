# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /powerapps/share/cmake/bin/cmake

# The command to remove a file.
RM = /powerapps/share/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /groups/itay_mayrose/halabikeren/biopp/bppsuite

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /groups/itay_mayrose/halabikeren/biopp/bppsuite/build

# Include any dependencies generated for this target.
include bppSuite/CMakeFiles/bppancestor.dir/depend.make

# Include the progress variables for this target.
include bppSuite/CMakeFiles/bppancestor.dir/progress.make

# Include the compile flags for this target's objects.
include bppSuite/CMakeFiles/bppancestor.dir/flags.make

bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o: bppSuite/CMakeFiles/bppancestor.dir/flags.make
bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o: ../bppSuite/bppAncestor.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o"
	cd /groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite && /powerapps/share/mpi/openmpi-1.10.4.c7/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bppancestor.dir/bppAncestor.cpp.o -c /groups/itay_mayrose/halabikeren/biopp/bppsuite/bppSuite/bppAncestor.cpp

bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bppancestor.dir/bppAncestor.cpp.i"
	cd /groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite && /powerapps/share/mpi/openmpi-1.10.4.c7/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /groups/itay_mayrose/halabikeren/biopp/bppsuite/bppSuite/bppAncestor.cpp > CMakeFiles/bppancestor.dir/bppAncestor.cpp.i

bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bppancestor.dir/bppAncestor.cpp.s"
	cd /groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite && /powerapps/share/mpi/openmpi-1.10.4.c7/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /groups/itay_mayrose/halabikeren/biopp/bppsuite/bppSuite/bppAncestor.cpp -o CMakeFiles/bppancestor.dir/bppAncestor.cpp.s

bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.requires:

.PHONY : bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.requires

bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.provides: bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.requires
	$(MAKE) -f bppSuite/CMakeFiles/bppancestor.dir/build.make bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.provides.build
.PHONY : bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.provides

bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.provides.build: bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o


# Object files for target bppancestor
bppancestor_OBJECTS = \
"CMakeFiles/bppancestor.dir/bppAncestor.cpp.o"

# External object files for target bppancestor
bppancestor_EXTERNAL_OBJECTS =

bppSuite/bppancestor: bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o
bppSuite/bppancestor: bppSuite/CMakeFiles/bppancestor.dir/build.make
bppSuite/bppancestor: /groups/itay_mayrose/halabikeren/biopp/lib64/libbpp-phyl.so.12.0.0
bppSuite/bppancestor: /groups/itay_mayrose/halabikeren/biopp/lib64/libbpp-popgen.so.8.0.0
bppSuite/bppancestor: /groups/itay_mayrose/halabikeren/biopp/lib64/libbpp-seq.so.12.0.0
bppSuite/bppancestor: /groups/itay_mayrose/halabikeren/biopp/lib64/libbpp-core.so.4.1.0
bppSuite/bppancestor: bppSuite/CMakeFiles/bppancestor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bppancestor"
	cd /groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bppancestor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
bppSuite/CMakeFiles/bppancestor.dir/build: bppSuite/bppancestor

.PHONY : bppSuite/CMakeFiles/bppancestor.dir/build

bppSuite/CMakeFiles/bppancestor.dir/requires: bppSuite/CMakeFiles/bppancestor.dir/bppAncestor.cpp.o.requires

.PHONY : bppSuite/CMakeFiles/bppancestor.dir/requires

bppSuite/CMakeFiles/bppancestor.dir/clean:
	cd /groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite && $(CMAKE_COMMAND) -P CMakeFiles/bppancestor.dir/cmake_clean.cmake
.PHONY : bppSuite/CMakeFiles/bppancestor.dir/clean

bppSuite/CMakeFiles/bppancestor.dir/depend:
	cd /groups/itay_mayrose/halabikeren/biopp/bppsuite/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /groups/itay_mayrose/halabikeren/biopp/bppsuite /groups/itay_mayrose/halabikeren/biopp/bppsuite/bppSuite /groups/itay_mayrose/halabikeren/biopp/bppsuite/build /groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite /groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/CMakeFiles/bppancestor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : bppSuite/CMakeFiles/bppancestor.dir/depend


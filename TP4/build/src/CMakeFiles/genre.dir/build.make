# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/local.isima.fr/imbenayad/shared/TPGEO/TP4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build

# Include any dependencies generated for this target.
include src/CMakeFiles/genre.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/genre.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/genre.dir/flags.make

src/CMakeFiles/genre.dir/genre.cpp.o: src/CMakeFiles/genre.dir/flags.make
src/CMakeFiles/genre.dir/genre.cpp.o: ../src/genre.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/genre.dir/genre.cpp.o"
	cd /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/genre.dir/genre.cpp.o -c /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/src/genre.cpp

src/CMakeFiles/genre.dir/genre.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/genre.dir/genre.cpp.i"
	cd /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/src/genre.cpp > CMakeFiles/genre.dir/genre.cpp.i

src/CMakeFiles/genre.dir/genre.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/genre.dir/genre.cpp.s"
	cd /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/src/genre.cpp -o CMakeFiles/genre.dir/genre.cpp.s

# Object files for target genre
genre_OBJECTS = \
"CMakeFiles/genre.dir/genre.cpp.o"

# External object files for target genre
genre_EXTERNAL_OBJECTS =

genre: src/CMakeFiles/genre.dir/genre.cpp.o
genre: src/CMakeFiles/genre.dir/build.make
genre: /usr/lib/x86_64-linux-gnu/libgmpxx.so
genre: /usr/lib/x86_64-linux-gnu/libmpfr.so
genre: /usr/lib/x86_64-linux-gnu/libgmp.so
genre: src/CMakeFiles/genre.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../genre"
	cd /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/genre.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/genre.dir/build: genre

.PHONY : src/CMakeFiles/genre.dir/build

src/CMakeFiles/genre.dir/clean:
	cd /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/src && $(CMAKE_COMMAND) -P CMakeFiles/genre.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/genre.dir/clean

src/CMakeFiles/genre.dir/depend:
	cd /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/local.isima.fr/imbenayad/shared/TPGEO/TP4 /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/src /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/src /home/local.isima.fr/imbenayad/shared/TPGEO/TP4/build/src/CMakeFiles/genre.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/genre.dir/depend


# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser"

# Include any dependencies generated for this target.
include Visualiser/CMakeFiles/Visualise.dir/depend.make

# Include the progress variables for this target.
include Visualiser/CMakeFiles/Visualise.dir/progress.make

# Include the compile flags for this target's objects.
include Visualiser/CMakeFiles/Visualise.dir/flags.make

Visualiser/CMakeFiles/Visualise.dir/main.cpp.o: Visualiser/CMakeFiles/Visualise.dir/flags.make
Visualiser/CMakeFiles/Visualise.dir/main.cpp.o: Visualiser/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Visualiser/CMakeFiles/Visualise.dir/main.cpp.o"
	cd "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Visualise.dir/main.cpp.o -c "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser/main.cpp"

Visualiser/CMakeFiles/Visualise.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Visualise.dir/main.cpp.i"
	cd "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser/main.cpp" > CMakeFiles/Visualise.dir/main.cpp.i

Visualiser/CMakeFiles/Visualise.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Visualise.dir/main.cpp.s"
	cd "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser" && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser/main.cpp" -o CMakeFiles/Visualise.dir/main.cpp.s

# Object files for target Visualise
Visualise_OBJECTS = \
"CMakeFiles/Visualise.dir/main.cpp.o"

# External object files for target Visualise
Visualise_EXTERNAL_OBJECTS =

Visualiser/Visualise: Visualiser/CMakeFiles/Visualise.dir/main.cpp.o
Visualiser/Visualise: Visualiser/CMakeFiles/Visualise.dir/build.make
Visualiser/Visualise: Audio-Processing/libAudio.a
Visualiser/Visualise: Visualiser/CMakeFiles/Visualise.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Visualise"
	cd "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Visualise.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Visualiser/CMakeFiles/Visualise.dir/build: Visualiser/Visualise

.PHONY : Visualiser/CMakeFiles/Visualise.dir/build

Visualiser/CMakeFiles/Visualise.dir/clean:
	cd "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser" && $(CMAKE_COMMAND) -P CMakeFiles/Visualise.dir/cmake_clean.cmake
.PHONY : Visualiser/CMakeFiles/Visualise.dir/clean

Visualiser/CMakeFiles/Visualise.dir/depend:
	cd "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser" "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser" "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser" "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser" "/mnt/c/Users/Rufus Vijayaratnam/Dev/Audio Visualiser/Visualiser/CMakeFiles/Visualise.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : Visualiser/CMakeFiles/Visualise.dir/depend


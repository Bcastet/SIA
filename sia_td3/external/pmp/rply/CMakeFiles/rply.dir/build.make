# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_SOURCE_DIR = /home/xenesis/SIA/sia_td3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xenesis/SIA/sia_td3

# Include any dependencies generated for this target.
include external/pmp/rply/CMakeFiles/rply.dir/depend.make

# Include the progress variables for this target.
include external/pmp/rply/CMakeFiles/rply.dir/progress.make

# Include the compile flags for this target's objects.
include external/pmp/rply/CMakeFiles/rply.dir/flags.make

external/pmp/rply/CMakeFiles/rply.dir/rply.c.o: external/pmp/rply/CMakeFiles/rply.dir/flags.make
external/pmp/rply/CMakeFiles/rply.dir/rply.c.o: external/pmp/rply/rply.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xenesis/SIA/sia_td3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object external/pmp/rply/CMakeFiles/rply.dir/rply.c.o"
	cd /home/xenesis/SIA/sia_td3/external/pmp/rply && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/rply.dir/rply.c.o   -c /home/xenesis/SIA/sia_td3/external/pmp/rply/rply.c

external/pmp/rply/CMakeFiles/rply.dir/rply.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/rply.dir/rply.c.i"
	cd /home/xenesis/SIA/sia_td3/external/pmp/rply && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/xenesis/SIA/sia_td3/external/pmp/rply/rply.c > CMakeFiles/rply.dir/rply.c.i

external/pmp/rply/CMakeFiles/rply.dir/rply.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/rply.dir/rply.c.s"
	cd /home/xenesis/SIA/sia_td3/external/pmp/rply && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/xenesis/SIA/sia_td3/external/pmp/rply/rply.c -o CMakeFiles/rply.dir/rply.c.s

rply: external/pmp/rply/CMakeFiles/rply.dir/rply.c.o
rply: external/pmp/rply/CMakeFiles/rply.dir/build.make

.PHONY : rply

# Rule to build all files generated by this target.
external/pmp/rply/CMakeFiles/rply.dir/build: rply

.PHONY : external/pmp/rply/CMakeFiles/rply.dir/build

external/pmp/rply/CMakeFiles/rply.dir/clean:
	cd /home/xenesis/SIA/sia_td3/external/pmp/rply && $(CMAKE_COMMAND) -P CMakeFiles/rply.dir/cmake_clean.cmake
.PHONY : external/pmp/rply/CMakeFiles/rply.dir/clean

external/pmp/rply/CMakeFiles/rply.dir/depend:
	cd /home/xenesis/SIA/sia_td3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xenesis/SIA/sia_td3 /home/xenesis/SIA/sia_td3/external/pmp/rply /home/xenesis/SIA/sia_td3 /home/xenesis/SIA/sia_td3/external/pmp/rply /home/xenesis/SIA/sia_td3/external/pmp/rply/CMakeFiles/rply.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/pmp/rply/CMakeFiles/rply.dir/depend

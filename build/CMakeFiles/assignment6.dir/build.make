# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.22.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.22.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/build

# Include any dependencies generated for this target.
include CMakeFiles/assignment6.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/assignment6.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/assignment6.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/assignment6.dir/flags.make

CMakeFiles/assignment6.dir/src/main.cpp.o: CMakeFiles/assignment6.dir/flags.make
CMakeFiles/assignment6.dir/src/main.cpp.o: ../src/main.cpp
CMakeFiles/assignment6.dir/src/main.cpp.o: CMakeFiles/assignment6.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/assignment6.dir/src/main.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/assignment6.dir/src/main.cpp.o -MF CMakeFiles/assignment6.dir/src/main.cpp.o.d -o CMakeFiles/assignment6.dir/src/main.cpp.o -c /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/src/main.cpp

CMakeFiles/assignment6.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/assignment6.dir/src/main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/src/main.cpp > CMakeFiles/assignment6.dir/src/main.cpp.i

CMakeFiles/assignment6.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/assignment6.dir/src/main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/src/main.cpp -o CMakeFiles/assignment6.dir/src/main.cpp.s

# Object files for target assignment6
assignment6_OBJECTS = \
"CMakeFiles/assignment6.dir/src/main.cpp.o"

# External object files for target assignment6
assignment6_EXTERNAL_OBJECTS =

assignment6: CMakeFiles/assignment6.dir/src/main.cpp.o
assignment6: CMakeFiles/assignment6.dir/build.make
assignment6: /Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd
assignment6: libimguizmo.a
assignment6: imgui/libimgui.a
assignment6: glfw/src/libglfw3.a
assignment6: glad/libglad.a
assignment6: CMakeFiles/assignment6.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable assignment6"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/assignment6.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/assignment6.dir/build: assignment6
.PHONY : CMakeFiles/assignment6.dir/build

CMakeFiles/assignment6.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/assignment6.dir/cmake_clean.cmake
.PHONY : CMakeFiles/assignment6.dir/clean

CMakeFiles/assignment6.dir/depend:
	cd /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6 /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6 /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/build /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/build /Users/guangzhu/Desktop/ETHz/S2/Shape_Modelling/gp22-Mere19/assignment6/build/CMakeFiles/assignment6.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/assignment6.dir/depend


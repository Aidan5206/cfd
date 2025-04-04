# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/aidanhsu/Documents/cfd

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/aidanhsu/Documents/cfd

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "No interactive CMake dialog available..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake to regenerate build system..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Install the project..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Install the project..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing only the local directory..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing only the local directory..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing the project stripped..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Installing the project stripped..."
	/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/aidanhsu/Documents/cfd/CMakeFiles /Users/aidanhsu/Documents/cfd//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/aidanhsu/Documents/cfd/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named cfd

# Build rule for target.
cfd: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 cfd
.PHONY : cfd

# fast build rule for target.
cfd/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cfd.dir/build.make CMakeFiles/cfd.dir/build
.PHONY : cfd/fast

#=============================================================================
# Target rules for targets named minimal

# Build rule for target.
minimal: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 minimal
.PHONY : minimal

# fast build rule for target.
minimal/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/minimal.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/minimal.dir/build
.PHONY : minimal/fast

#=============================================================================
# Target rules for targets named basic

# Build rule for target.
basic: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 basic
.PHONY : basic

# fast build rule for target.
basic/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/basic.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/basic.dir/build
.PHONY : basic/fast

#=============================================================================
# Target rules for targets named modern

# Build rule for target.
modern: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 modern
.PHONY : modern

# fast build rule for target.
modern/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/modern.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/modern.dir/build
.PHONY : modern/fast

#=============================================================================
# Target rules for targets named animation

# Build rule for target.
animation: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 animation
.PHONY : animation

# fast build rule for target.
animation/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/animation.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/animation.dir/build
.PHONY : animation/fast

#=============================================================================
# Target rules for targets named nonblock

# Build rule for target.
nonblock: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 nonblock
.PHONY : nonblock

# fast build rule for target.
nonblock/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/nonblock.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/nonblock.dir/build
.PHONY : nonblock/fast

#=============================================================================
# Target rules for targets named xkcd

# Build rule for target.
xkcd: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 xkcd
.PHONY : xkcd

# fast build rule for target.
xkcd/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/xkcd.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/xkcd.dir/build
.PHONY : xkcd/fast

#=============================================================================
# Target rules for targets named bar

# Build rule for target.
bar: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 bar
.PHONY : bar

# fast build rule for target.
bar/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/bar.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/bar.dir/build
.PHONY : bar/fast

#=============================================================================
# Target rules for targets named fill_inbetween

# Build rule for target.
fill_inbetween: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 fill_inbetween
.PHONY : fill_inbetween

# fast build rule for target.
fill_inbetween/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/fill_inbetween.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/fill_inbetween.dir/build
.PHONY : fill_inbetween/fast

#=============================================================================
# Target rules for targets named fill

# Build rule for target.
fill: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 fill
.PHONY : fill

# fast build rule for target.
fill/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/fill.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/fill.dir/build
.PHONY : fill/fast

#=============================================================================
# Target rules for targets named update

# Build rule for target.
update: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 update
.PHONY : update

# fast build rule for target.
update/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/update.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/update.dir/build
.PHONY : update/fast

#=============================================================================
# Target rules for targets named subplot2grid

# Build rule for target.
subplot2grid: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 subplot2grid
.PHONY : subplot2grid

# fast build rule for target.
subplot2grid/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/subplot2grid.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/subplot2grid.dir/build
.PHONY : subplot2grid/fast

#=============================================================================
# Target rules for targets named lines3d

# Build rule for target.
lines3d: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 lines3d
.PHONY : lines3d

# fast build rule for target.
lines3d/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/lines3d.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/lines3d.dir/build
.PHONY : lines3d/fast

#=============================================================================
# Target rules for targets named surface

# Build rule for target.
surface: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 surface
.PHONY : surface

# fast build rule for target.
surface/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/surface.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/surface.dir/build
.PHONY : surface/fast

#=============================================================================
# Target rules for targets named colorbar

# Build rule for target.
colorbar: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 colorbar
.PHONY : colorbar

# fast build rule for target.
colorbar/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/colorbar.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/colorbar.dir/build
.PHONY : colorbar/fast

#=============================================================================
# Target rules for targets named contour

# Build rule for target.
contour: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 contour
.PHONY : contour

# fast build rule for target.
contour/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/contour.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/contour.dir/build
.PHONY : contour/fast

#=============================================================================
# Target rules for targets named spy

# Build rule for target.
spy: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 spy
.PHONY : spy

# fast build rule for target.
spy/fast:
	$(MAKE) $(MAKESILENT) -f _deps/matplotlibcpp-build/CMakeFiles/spy.dir/build.make _deps/matplotlibcpp-build/CMakeFiles/spy.dir/build
.PHONY : spy/fast

Grid.o: Grid.cpp.o
.PHONY : Grid.o

# target to build an object file
Grid.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cfd.dir/build.make CMakeFiles/cfd.dir/Grid.cpp.o
.PHONY : Grid.cpp.o

Grid.i: Grid.cpp.i
.PHONY : Grid.i

# target to preprocess a source file
Grid.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cfd.dir/build.make CMakeFiles/cfd.dir/Grid.cpp.i
.PHONY : Grid.cpp.i

Grid.s: Grid.cpp.s
.PHONY : Grid.s

# target to generate assembly for a file
Grid.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cfd.dir/build.make CMakeFiles/cfd.dir/Grid.cpp.s
.PHONY : Grid.cpp.s

main.o: main.cpp.o
.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cfd.dir/build.make CMakeFiles/cfd.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i
.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cfd.dir/build.make CMakeFiles/cfd.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s
.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cfd.dir/build.make CMakeFiles/cfd.dir/main.cpp.s
.PHONY : main.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... list_install_components"
	@echo "... rebuild_cache"
	@echo "... animation"
	@echo "... bar"
	@echo "... basic"
	@echo "... cfd"
	@echo "... colorbar"
	@echo "... contour"
	@echo "... fill"
	@echo "... fill_inbetween"
	@echo "... lines3d"
	@echo "... minimal"
	@echo "... modern"
	@echo "... nonblock"
	@echo "... spy"
	@echo "... subplot2grid"
	@echo "... surface"
	@echo "... update"
	@echo "... xkcd"
	@echo "... Grid.o"
	@echo "... Grid.i"
	@echo "... Grid.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system


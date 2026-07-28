#-------------------------------------------------------------------------------
# Makefile for the SuiteSparse CMake support layer
#-------------------------------------------------------------------------------

# SPDX-FileCopyrightText: 2026 Sergiu Deitsch
# SPDX-License-Identifier: Apache-2.0

#-------------------------------------------------------------------------------

export CMAKE_OPTIONS ?=
export JOBS ?= 8

BUILD_DIR ?= build

.PHONY: all configure library local global debug demos tests test install \
    uninstall clean purge distclean

default: library

all: library

configure:
	cmake -S . -B $(BUILD_DIR) $(CMAKE_OPTIONS)

library: configure
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)

local:
	cmake -S . -B $(BUILD_DIR) $(CMAKE_OPTIONS) \
		-DCMAKE_INSTALL_PREFIX=$(CURDIR)
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)

global:
	cmake -S . -B $(BUILD_DIR) $(CMAKE_OPTIONS)
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)

debug:
	cmake -S . -B $(BUILD_DIR) $(CMAKE_OPTIONS) \
		-DCMAKE_BUILD_TYPE=Debug
	cmake --build $(BUILD_DIR) --config Debug --parallel $(JOBS)

demos:
	cmake -S . -B $(BUILD_DIR) $(CMAKE_OPTIONS) -DWITH_DEMOS=ON
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)

tests:
	cmake -S . -B $(BUILD_DIR) $(CMAKE_OPTIONS) -DBUILD_TESTING=ON
	cmake --build $(BUILD_DIR) --config Release --parallel $(JOBS)
	ctest --test-dir $(BUILD_DIR) --output-on-failure

test: tests

install:
	cmake --install $(BUILD_DIR)

uninstall:
	-$(RM) $(BUILD_DIR)/install_manifest.txt

clean:
	cmake -E rm -rf $(BUILD_DIR)

purge: clean

distclean: clean


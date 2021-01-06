WARTHOG_SRC = $(wildcard domains/*.cpp) $(wildcard util/*.cpp) \
			  $(wildcard search/*.cpp) $(wildcard experimental/*.cpp) \
			  $(wildcard heuristics/*.cpp) $(wildcard jps/*.cpp) \
			  $(wildcard contraction/*.cpp) $(wildcard label/*.cpp) \
			  $(wildcard memory/*.cpp) $(wildcard mapf/*.cpp) $(wildcard sys/*.cpp) \
			  $(wildcard sipp/*.cpp) $(wildcard cpd/*.cpp)

WARTHOG_EXE = $(wildcard programs/*.cpp)
WARTHOG_TEST = $(wildcard test/*.cpp)
WARTHOG_OBJ = $(WARTHOG_SRC:%.cpp=build/%.o)

D_WARTHOG_INCLUDES = -I../ -I../domains -I../util -I../search -I../experimental \
					 -I../heuristics -I../jps -I../contraction -I../label \
					 -I../memory -I../mapf -I../sys -I../sipp -I../cpd -I../third_party

D_INCLUDES = $(D_WARTHOG_INCLUDES) -I/usr/include -I/usr/local/include
D_LIBS = -L./lib -L/usr/local/lib

warthog = ./lib/libwarthog.a

CFLAGS = -std=c++11 -pedantic -pthread -Wall -Werror -Wno-conversion \
	-Wno-unused-result -Wno-unused-but-set-variable -fopenmp
# PROFILE_CFLAGS = $(DEV_CFLAGS) -pg -DNDEBUG

FLAVOURS = fast dev debug
PROGRAMS = $(WARTHOG_EXE:programs/%.cpp=bin/%)
PROGRAMS += $(WARTHOG_TEST:.cpp=)

ifeq ("$(findstring Darwin, "$(shell uname -s)")", "Darwin")
  CC = g++
  CFLAGS += -DOS_MAC
else
  CC = g++
  ifeq ("$(findstring Linux, "$(shell uname -s)")", "Linux")
    D_LIBS += -lrt
  endif
endif

.EXPORT_ALL_VARIABLES:

.PHONY: all help default clean cleaner $(FLAVOURS)

.SILENT: all clean cleaner

default:
	+$(MAKE) fast main

##@ Options

# These need to be explicitly split
all:														## Build all flavours
	for dir in $(FLAVOURS); do \
		$(MAKE) $$dir; \
	done

help:														## Print this message
	@echo "Compile warthog in different flavours, the executables will be built"
	@echo "in their respective directory."
	@echo ""
	@echo "Usage: make [cmd | flavour [target]]"
	@awk 'BEGIN {FS = ":.*##";} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)
	@echo "$(PROGRAMS)"

##@ Cleaning

clean:													## Clean build artifacts, can be run as a flavour
	for dir in $(FLAVOURS); do \
		$(MAKE) -C $$dir clean; \
	done

cleaner:												## Remove sub-directories
	for dir in $(FLAVOURS); do \
		$(RM) -rf $$dir; \
	done

##@ Flavours

# We need to call make explicitly on the subdirs as they are the ones deciding
# what needs to run.
fast: CFLAGS += -O3 -DNDEBUG -Wno-unused-variable
fast: fast/Makefile							## Compile with opti flags
# Little hack: the first argument is read as the subdir
	+$(MAKE) -C $(MAKECMDGOALS)

dev: CFLAGS += -ggdb -O0 -gdwarf-4 -no-pie -DSINGLE_THREADED -DQUICK_LOG=1
dev: dev/Makefile								## Compile with debug symbols
	+$(MAKE) -C $(MAKECMDGOALS)

debug: CFLAGS += -ggdb -O3 -gdwarf-4 -no-pie -DQUICK_LOG=3
debug: debug/Makefile						## Opti flags with debug symbols
	+$(MAKE) -C $(MAKECMDGOALS)

# Generate the Makefile of the appropriate flavour
%/Makefile: make.file
	@$(shell mkdir -p $(@D))
	@$(shell cp $< $@)

# Gobble extra targets
%:
	@true

# Used for help
##@ Targets

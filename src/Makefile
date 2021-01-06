WARTHOG_SRC = $(wildcard domains/*.cpp) $(wildcard util/*.cpp) \
			  $(wildcard search/*.cpp) $(wildcard experimental/*.cpp) \
			  $(wildcard heuristics/*.cpp) $(wildcard jps/*.cpp) \
			  $(wildcard contraction/*.cpp) $(wildcard label/*.cpp) \
			  $(wildcard memory/*.cpp) $(wildcard mapf/*.cpp) $(wildcard sys/*.cpp) \
			  $(wildcard sipp/*.cpp) $(wildcard cpd/*.cpp)


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

.PHONY: all info dev fast default

default:
	+$(MAKE) fast main

# These need to be explicitly split
all:
	+$(MAKE) fast
	+$(MAKE) dev

fast: CFLAGS += -O3 -DNDEBUG -Wno-unused-variable
fast: fast/Makefile
# Little hack: the first argument is read as the subdir
	+$(MAKE) -C $(MAKECMDGOALS)

dev: CFLAGS += -ggdb -O0 -gdwarf-4 -no-pie -DSINGLE_THREADED -DQUICK_LOG=1
dev: dev/Makefile
	+$(MAKE) -C $(MAKECMDGOALS)

fast/Makefile: make.file
	@$(shell mkdir -p $(@D))
	@$(shell cp $< $@)

dev/Makefile: make.file
	@$(shell mkdir -p $(@D))
	@$(shell cp $< $@)

info:
	@echo "Compile warthog in different flavours, the executables will be built"
	@echo "in their respective directory."
	@echo ""
	@echo "Usage:   make [fast | dev] <target>"
	@echo "Options:"
	@echo "  fast          build with opti flags"
	@echo "  dev           build with debug info"

# Gobble extra targets
%:
	@true

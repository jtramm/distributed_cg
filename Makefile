#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPENMP      = yes
MPI         = yes
OPTIMIZE    = yes
PROFILE     = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = conj

source = \
main.c \
conj.c \
utils.c

obj = $(source:.c=.o)

#===============================================================================
# Sets Flags
#===============================================================================

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
endif

# Standard Flags
CFLAGS := -std=gnu99

# MPI Compiler
ifeq ($(MPI),yes)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g -fno-omit-frame-pointer
  #CFLAGS += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  #CFLAGS += -pg -O0 -fno-omit-frame-pointer
  CFLAGS += -pg -fno-omit-frame-pointer
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

# OpenMP
ifeq ($(OPENMP),yes)
  CFLAGS += -fopenmp
endif

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) conj_header.h Makefile
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c Makefile conj_header.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) *\.lst conj.dSYM

edit:
	vim -p $(source) conj_header.h

run:
	./$(program)

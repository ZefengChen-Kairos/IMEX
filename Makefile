# =============================================================================
#  Semi-implicit IMEX finite-volume MHD solver with constrained transport
# =============================================================================
#  make            build (optimised) -> build/imex_mhd
#  make debug      build with -O0 -g -fsanitize=address,undefined
#  make run        build and run in ./
#  make clean      remove build artefacts
#
#  Eigen and Boost.Math are header-only. If they live somewhere unusual:
#      make EIGEN_INC=/path/to/eigen3 BOOST_INC=/path/to/boost/parent
# =============================================================================

TARGET   := imex_mhd
BUILDDIR := build
SRCDIR   := src
INCDIR   := include

CXX      ?= g++
CXXSTD   := -std=c++17
WARN     := -Wall -Wextra -Wno-sign-compare
OPT      := -O2

# --- header-only dependency discovery -----------------------------------------
# Checked in order; first hit wins. Override on the command line if needed.
EIGEN_CANDIDATES := /usr/include/eigen3 \
                    /usr/local/include/eigen3 \
                    /opt/homebrew/include/eigen3 \
                    $(HOME)/include/eigen3
BOOST_CANDIDATES := /usr/include \
                    /usr/local/include \
                    /opt/homebrew/include

EIGEN_INC ?= $(firstword $(foreach d,$(EIGEN_CANDIDATES),$(if $(wildcard $(d)/Eigen/Dense),$(d))))
BOOST_INC ?= $(firstword $(foreach d,$(BOOST_CANDIDATES),$(if $(wildcard $(d)/boost/math/quadrature/trapezoidal.hpp),$(d))))

ifeq ($(strip $(EIGEN_INC)),)
  $(error Eigen not found. Install it (apt install libeigen3-dev / brew install eigen) or pass EIGEN_INC=/path/to/eigen3)
endif
ifeq ($(strip $(BOOST_INC)),)
  $(error Boost.Math not found. Install it (apt install libboost-math-dev / brew install boost) or pass BOOST_INC=/path/to/include)
endif

CPPFLAGS := -I$(INCDIR) -I$(EIGEN_INC) -I$(BOOST_INC) -MMD -MP
CXXFLAGS = $(CXXSTD) $(WARN) $(OPT)
LDFLAGS  :=
LDLIBS   :=

# GCC 8 and earlier need this for <filesystem>
ifeq ($(shell $(CXX) -dumpversion | cut -d. -f1),8)
  LDLIBS += -lstdc++fs
endif

SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

.PHONY: all debug run clean help

all: $(BUILDDIR)/$(TARGET)

$(BUILDDIR)/$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	@echo "built $@"

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

debug: OPT := -O0 -g -fsanitize=address,undefined
debug: LDFLAGS += -fsanitize=address,undefined
debug: clean all


run: all
	./$(BUILDDIR)/$(TARGET)

clean:
	$(RM) -r $(BUILDDIR)

help:
	@echo "targets: all (default), debug, run, clean"
	@echo "EIGEN_INC = $(EIGEN_INC)"
	@echo "BOOST_INC = $(BOOST_INC)"

-include $(DEPS)

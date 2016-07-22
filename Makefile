###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# MOOSE_DIR    := /Users/kerblm/projects/moose
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
ALL_MODULES := yes
include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# Serpent dep app
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := SerpentCouplerII
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk
# ADDITIONAL_CXXFLAGS += -fopenmp -lgd -lm -omp

# Serpent (optional)
# SERPENT_DIR          ?= $(CURDIR)/contrib/Serpent
# ifneq ($(wildcard $(SERPENT_DIR)/src/main.c),)
#  APPLICATION_DIR    := $(SERPENT_DIR)
#  APPLICATION_NAME   := serpent
#  DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
#  include            $(FRAMEWORK_DIR)/app.mk
#  ADDITIONAL_CXXFLAGS += -fopenmp -lgd -lm
# endif

#include            $(CURDIR)/serpent.mk

###############################################################################
# Additional special case targets should be added here

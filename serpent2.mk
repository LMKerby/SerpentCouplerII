# Serpent coupling
# Created by Leslie Kerby, 5/2016

SERPENT_DIR := $(CURDIR)/contrib/Serpent
SERPENT_SRC_DIR := $(CURDIR)/contrib/Serpent/src
SERPENT_INC_DIR := $(CURDIR)/contrib/Serpent/include
serpent_LIB := $(SERPENT_DIR)/libserpent-$(METHOD).la
serpent_APP := $(SERPENT_DIR)/serpent-$(METHOD)-dyn

#serpent_fsrcfiles    := $(shell find $(SERPENT_SRC_DIR) -name "*.f90" -not -name "main*" -not -name "serpent_mesh_generator*" -not -name "ampx2xml*")
serpent_csrcfiles    := $(shell find $(SERPENT_SRC_DIR) -name *.c)
serpent_mainsrcfiles := $(shell find $(SERPENT_SRC_DIR) -name main*.c)

serpent_objects     := $(patsubst %.c, %.$(obj-suffix), $(serpent_csrcfiles))
#serpent_objects     += $(patsubst %.c, %.$(obj-suffix), $(serpent_csrcfiles))
serpent_app_objects := $(patsubst %.c, %.$(obj-suffix), $(SERPENT_SRC_DIR)/main.c)

serpent_INCLUDE := $(foreach i, $(SERPENT_INC_DIR), -I$(i))
libmesh_INCLUDE := $(serpent_INCLUDE) $(libmesh_INCLUDE)

serpent_FLAGS := -fopenmp -lgd -lm

ifneq ($(wildcard $(CURDIR)/contrib/triangle/triangle.c),)
  TRIANGLE_DIR ?= $(CURDIR)/contrib/triangle
  TRIANGLE_OBJ := $(TRIANGLE_DIR)/triangle.o
  yak_FFLAGS += -Duse_Triangle
$(TRIANGLE_OBJ):
	@echo "Making Triangle ..."
	cd $(TRIANGLE_DIR) && make trilibrary
endif


xsgenclean:
	-rm -fr  $(XSGEN_APP) $(XSGEN_LIB) $(XSGEN_DIR)/.libs  $(XSGEN_DIR)/*dylib $(XSGEN_DIR)/src/.libs $(XSGEN_DIR)/*.lo* $(XSGEN_SRC_DIR)/*.lo*

serpent: $(serpent_LIB)

$(serpent_LIB): $(serpent_objects)
	@echo "Linking library "$@"..."
	@echo "serpent_INCLUDE = "$(serpent_INCLUDE)
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(serpent_objects) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(serpent_FLAGS) -rpath $(SERPENT_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(serpent_LIB) $(SERPENT_DIR)

serpent: $(serpent_APP)

$(serpent_APP): $(serpent_LIB) $(serpent_app_objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
  $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(serpent_app_objects) $(serpent_LIB) $(libmesh_LIBS) $(libmesh_LDFLAGS)

cleanserp:
	$(RM) $(serpent_LIB) $(serpent_objects) $(serpent_APP)

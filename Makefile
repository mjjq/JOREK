# This makefile automatically searches and compiles fortran programs
# Known limitations:
# - All dependencies must be in modules
# - Module name must be equal to file name
# - Dependencies do not work correctly on make < 3.81

# Include jorek-specific things and settings
include Makefile.inc
default: jorek2_main
include defaults.mk

.PHONY: .mod/version.h clean cleanall cleandep
cleanall: clean cleandep
clean:
	@echo ">> Deleting Object Files <<"
	-@rm -r $(OBJDIR)
	@echo ">> Deleting Module Files <<"
	-@rm -r $(MODDIR)
	-@find . -name '*.mod' -or -name '*.dep' -or -name '*.o' -delete

cleandep:
	@echo ">> Deleting Dependency Files <<"
	-@rm -r $(DEPDIR)


# Most of the sources are actually fortran 2003/8 but use the .f90 suffix
DIRS =  .				\
	timing				\
        models/$(MODEL) 		\
        datatypes			\
	models				\
	communication			\
	elements			\
	grids				\
	grids/grid_utils		\
	matrix				\
	solvers 			\
	plots				\
	diagnostics			\
	diagnostics/new_diag		\
	refinement			\
	postproc			\
	tools                           \
	util                            \
	vacuum

sources:=$(shell find $(DIRS) -maxdepth 1 -iname '*.f90')
depends:=$(patsubst %.f90,$(DEPDIR)/%.d,$(notdir $(sources)))
# check if there are any duplicate files
duplicates:=$(shell echo $(notdir $(sources)) | tr ' ' "\n" | sort | uniq -d)
ifneq ($(strip $(duplicates)),)
  $(error Duplicate filenames found: $(duplicates))
endif

# For each source file add an explicit rule with the template
$(foreach dir,$(sort $(dir $(sources))),$(eval $(call F90_O_TEMPLATE,$(dir))))
# For each source file add a rule to create dependency files
$(foreach dir,$(sort $(dir $(sources))),$(eval $(call F90_D_TEMPLATE,$(dir))))
# Cancel all rules that make a .f90 file, .d file, Makefile or Makefile.inc
# to avoid unnecessary checks
%.f90: ;
Makefile: ;
Makefile.inc: ;
%.d: ;
-include $(depends)



# A list of programs to compile by filename
PROGRAM_SOURCES = jorek2_main.f90                       \
		  diagnostics/jorek2_poincare.f90       \
		  diagnostics/rst_bin2hdf5.f90          \
		  diagnostics/rst_hdf52bin.f90          \
		  diagnostics/jorek2_four.f90           \
		  postproc/jorek2_postproc.f90          \
		  diagnostics/jorek2_povray.f90         \
		  diagnostics/jorek2_connection.f90     \
		  diagnostics/jorek2_connection2.f90    \
		  diagnostics/jorek2_strikes_ordered.f90\
		  diagnostics/enbiggen.f90              \
		  diagnostics/jordel.f90                \
		  diagnostics/jorpol.f90                \
		  diagnostics/jorek2vtk.f90             \
		  diagnostics/jorek2_diagno.f90         \
		  diagnostics/jorek_to_helena.f90       \
		  diagnostics/jorek2_target2vtk.f90     \
		  diagnostics/jorek2_powers.f90         \
		  diagnostics/new_diag_demo.f90 	\
		  diagnostics/jorek_extract_data.f90    \
		  diagnostics/jorek2_connection_stan.f90\
		  diagnostics/jorek2_fieldlines_vtk.f90 \
		  diagnostics/jorek2vtk_3d.f90          \
		  diagnostics/jorek2vtk_GausVortTerms.f90          \
		  diagnostics/jorek2_stan.f90           \
		  util/import_eqdsk.f90                 \
		  diagnostics/jorek2_import_perturbation.f90 \

$(foreach prog,$(PROGRAM_SOURCES) $(EXAMPLES),$(eval $(call PROGRAM_TEMPLATE,$(prog))))

# The all target compiles all of the programs above
all: $(notdir $(PROGRAM_SOURCES))

.mod/version.h:
	@echo "Generate version.h"
	@rm -f $@
	@`pwd`/util/version.sh 2>/dev/null > $@
	@echo "#define compile_command '$(FC)'" >> $@
	@echo "#define compile_flags '$(FFLAGS)'" >> $@
	@echo "#define compile_includes '$(INCLUDES)'" >> $@
	@echo "#define compile_defines '$(DEFINES)'" >> $@
	@echo "#define compile_libs '$(LIBS)'" >> $@
	-@echo "#define compile_dir '`pwd`'" >> $@
	-@echo "#define compile_time '`date \"+%F %T\"`'" >> $@
	-@echo "#define compile_user '`whoami`'" >> $@
	-@echo "#define compile_machine '`hostname`'" >> $@
	@echo "#define compile_modules '$(LOADEDMODULES)'" >> $@

# Special cases
jorek2_main.o: INCLUDES += -Itiming
#r3_info.o: INCLUDES += -Itiming # not sure if this one is neded, as it will probably look in .
# Is this used by anyone?
#-include forcheck.mk

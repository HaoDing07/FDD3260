
include $(PFUNIT_DIR)/include/PFUNIT.mk
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

all: test_shared_thermo test_activation 

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I$(BUILDDIR)
FFLAGS += -L$(BUILDDIR)

test_shared_thermo_TESTS := test_shared_thermo.pf 
test_shared_thermo_OTHER_LIBRARIES := -L$(BUILDDIR) -lsut
$(eval $(call make_pfunit_test,test_shared_thermo))

test_activation_TESTS := test_activation.pf 
test_activation_OTHER_LIBRARIES := -L$(BUILDDIR) -lsut
$(eval $(call make_pfunit_test,test_activation))

clean:
	$(RM) *.o *.mod *.a  *.inc
	$(RM) test_test_shared_thermo.F90 test_test_shared_thermo test_activation.F90 test_activation


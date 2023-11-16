
exe:
	$(CXX) -c $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) *.cpp
	$(LD) *.o -o $@ $(PLUMED_LOAD)

exe-c:
	$(CC) -c $(CPPFLAGS) $(ADDCPPFLAGS) *.c
	$(LD) *.o -o exe $(PLUMED_LOAD)

test-c11:
	$(CC) -c $(CPPFLAGS) $(ADDCPPFLAGS) __test_c11.c
	rm -f __test_c11.o
	@echo "SUCCESS=YES"

exe-fortran:
	$(FC) $(PLUMED_FORTRAN_FLAGS) -c $(PLUMED_FORTRAN) *.f90
	$(FC) $(PLUMED_FORTRAN_FLAGS) *.o -o exe $(PLUMED_LOAD)

exe-fortran08:
	$(FC) $(PLUMED_FORTRAN_FLAGS) -c $(PLUMED_FORTRAN08) *.f90
	$(FC) $(PLUMED_FORTRAN_FLAGS) *.o -o exe $(PLUMED_LOAD)

test-fortran08:
	$(FC) $(PLUMED_FORTRAN_FLAGS) -c __test_fortran08.f90
	rm -f __test_fortran08*
	@echo "SUCCESS=YES"

test-fortran-arg-mismatch:
	$(FC) $(PLUMED_FORTRAN_FLAGS) -c __test_fortran_arg_mismatch.f90
	@echo "SUCCESS=YES"

print-fortran:
	@echo FC=$(FC)

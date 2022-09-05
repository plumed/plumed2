
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
	$(FC) -c $(PLUMED_FORTRAN) *.f90
	$(FC) *.o -o exe $(PLUMED_LOAD)

exe-fortran08:
	$(FC) -c $(PLUMED_FORTRAN08) *.f90
	$(FC) *.o -o exe $(PLUMED_LOAD)

test-fortran08:
	$(FC) -c __test_fortran08.f90
	rm -f __test_fortran08*
	@echo "SUCCESS=YES"

print-fortran:
	@echo FC=$(FC)

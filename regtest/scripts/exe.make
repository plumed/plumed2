
exe:
	$(CXX) -c $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) *.cpp
	$(LD) *.o -o $@ $(PLUMED_LOAD)

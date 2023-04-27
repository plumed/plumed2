#include <iostream>

#include "plumed/wrapper/Plumed.h"
#include "plumed/tools/Communicator.h"


int main(int, char**){
    //GIVEN an installation without MPI
    if (PLMD::Communicator::PlumedHasMPI) {;
        return 0; //no problem: passing
    }
    //AND GIVEN an initialized Plumed interface 
    PLMD::Plumed plumed;
    //intitalize the interpreter
    unsigned int natoms=10;

    std::vector<double> positions(3*natoms,0.0);
    for(unsigned i=0;i<3*natoms;i++) positions[i]=i;
    std::vector<double> masses(natoms,1.0);
    std::vector<double> forces(3*natoms,0.0);
    std::vector<double> box(9,0.0);
    std::vector<double> virial(9,0.0);
    plumed.cmd("setNatoms",&natoms);
    plumed.cmd("setLogFile","test.log");
    plumed.cmd("init");
    plumed.cmd("readInputLine","d: DISTANCE ATOMS=1,2");
    plumed.cmd("readInputLine","d1: DISTANCE ATOMS={1 2}");
    //WHEN the user give an input file with WALKERS_MPI keyword in it
    
    const std::string mokedLine=
    "METAD ARG=d,d1 SIGMA=0.1,0.2 HEIGHT=0.1 PACE=2 RESTART=YES WALKERS_MPI";
    //THEN plumed should gracefully exit with a clear error message
    try{
        plumed.cmd("readInputLine",mokedLine.c_str());
    } catch(PLMD::Plumed::ExceptionError &e){//correct throw, we are happy
        std::string exceptionText{e.what()};
        if (exceptionText.find("requires MPI compilation") != std::string::npos) {
            //throws with the wanted message and we are happy
           return 0;
        }
        std::cout << "Exception thrown, wrong message: "
                  << e.what() << '\n';
        return 1;
    }
    std::cout << "Exception not thrown\n";
    return 1;
}
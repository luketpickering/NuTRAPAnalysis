#include <iostream>

#include "EventLoop.hxx"

int main(int argc, char const * argv[]){

  RunLoop(argc,argv,[](Long64_t ent, TransversityVarsB const &tv){
    std::cout << "[" << ent << "]" << tv.NFinalStateParticles << std::endl;
    return true;
  });

}

#ifndef __EVENTLOOP_HXX_SEEN__
#define __EVENTLOOP_HXX_SEEN__

#include <functional>

#include "TransversityVariableObjects.hxx"

int RunLoop(int argc,
  char const * argv[],
  std::function<bool(Long64_t LoopNum, TransversityVarsB const &tv)> EvLoopFunc);

#endif

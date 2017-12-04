#include "zeitgeist.hpp"    

/// Constructs a Zeitgeist 
Zeitgeist::Zeitgeist(){
      tges = 0.0;
}

/// Deletes the Zeitgeist
  Zeitgeist::~Zeitgeist(){
      
  }

//starts the timer
void Zeitgeist::Tic(){
    tstart = clock();
}

void Zeitgeist::Tac(){
    tges = tges + clock() - tstart;
}

//ends the timer
real_t Zeitgeist::Toc() const {
    return (real_t) (tges)/CLOCKS_PER_SEC;
}

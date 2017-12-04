//------------------------------------------------------------------------------
#ifndef __ZEITGEIST_HPP
#define __ZEITGEIST_HPP

//------------------------------------------------------------------------------

#include <iostream>
#include "typedef.hpp"
#include <time.h>

//------------------------------------------------------------------------------
class Zeitgeist {
public:
  /// Constructs a Zeitgeist 
  Zeitgeist();
  
  /// Deletes the Zeitgeist
  ~Zeitgeist();

  //starts the timer
  void Tic();
  
  //blablabladföhi
  void Tac();

  //ends the timer
  real_t Toc() const;
  
  
    //tend = time(0); 
    //std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< std::endl;
    
  

private:
  double tstart;
  double tges;
};
//------------------------------------------------------------------------------
#endif // __GRID_HPP

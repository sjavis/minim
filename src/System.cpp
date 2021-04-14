#include "System.h"


System::System(int ndof_) : ndof(ndof_) {
  state = new double[ndof_]();
}


System::System(int ndof_, double* state_)
 : ndof(ndof_), state(state_)
{}


System::~System() {
  //delete this->state;
}

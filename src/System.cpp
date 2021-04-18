#include <vector>
#include "System.h"


System::System(int ndof_)
: ndof(ndof_), state(ndof)
{}


System::System(int ndof_, std::vector<double> state_)
: ndof(ndof_), state(state_)
{}

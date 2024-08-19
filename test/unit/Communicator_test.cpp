#include "test_main.cpp"
#include "communicators/CommUnstructured.h"

#include "Potential.h"
#include "State.h"
#include "utils/vec.h"

using namespace minim;

// Define potential using energy elements to test CommUnstructured
class ParallelPot: public NewPotential<ParallelPot> {
  public:
    ParallelPot() { _parallelDef = true; };

    void init(const vector<double>& coords) {
      elements = {};
      elements.push_back({0, {0,1}});
      elements.push_back({0, {2,3}});
      elements.push_back({0, {4,5}});
      constraints = {};
      constraints.push_back({{1}});
      constraints.push_back({{4}});
      constraints.push_back({{1,4}});
    }

    void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) {
      switch (el.type) {
        case 0: {
          if (e) *e += coords[el.idof[0]] + coords[el.idof[1]];
          if (g) {
            (*g)[el.idof[0]] += 1;
            (*g)[el.idof[1]] += 1;
          }
        } break;
      }
    }
};


TEST(CommUnstructured, TestConstraintDistribution) {
  vector<double> coords{0, 1, 2, 3, 4, 5};
  ParallelPot pot;
  pot.init(coords);
  CommUnstructured comm;
  comm.setup(pot, 6, {0,1});

  EXPECT_EQ(pot.constraints.size(), 2);
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {1}));
  EXPECT_TRUE(ArraysMatch(pot.constraints[1].idof, (mpi.rank==0)?vector<int>{1,4}:vector<int>{3,1}));
}

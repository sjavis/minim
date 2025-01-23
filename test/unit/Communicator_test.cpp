// NPROCS 4
#include "test_main.cpp"
#include "communicators/CommUnstructured.h"
#include "communicators/CommGrid.h"

#include "Potential.h"
#include "State.h"
#include "utils/vec.h"
#include "utils/mpi.h"

using namespace minim;

// Define potential using energy elements to test CommUnstructured
class UnstructuredPot: public NewPotential<UnstructuredPot> {
  public:
    int potentialType() const override { return Potential::UNSTRUCTURED; };

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
  if (mpi.rank >= 2) return;
  vector<double> coords{0, 1, 2, 3, 4, 5};
  UnstructuredPot pot;
  pot.init(coords);
  CommUnstructured comm;
  comm.setup(pot, 6, {0,1});

  EXPECT_EQ(pot.constraints.size(), 2);
  EXPECT_TRUE(ArraysMatch(pot.constraints[0].idof, {1}));
  EXPECT_TRUE(ArraysMatch(pot.constraints[1].idof, (mpi.rank==0)?vector<int>{1,4}:vector<int>{3,1}));
}


// Define a potential to test CommGrid
class GridPot: public NewPotential<GridPot> {
  public:
    int potentialType() const override { return Potential::GRID; };
    GridPot(vector<int> gridSize) { this->gridSize = gridSize; }
};


TEST(CommGrid, TestCommArray) {
  if (mpi.rank >= 2) return;
  GridPot pot({2,6});

  // Default commArray
  CommGrid comm1;
  comm1.setup(pot, 12, {0,1});
  EXPECT_TRUE(ArraysMatch(comm1.commArray, {1, 2}));

  // commArray defined in the Communicator
  CommGrid comm2;
  comm2.commArray = {2, 1};
  comm2.setup(pot, 12, {0,1});
  EXPECT_TRUE(ArraysMatch(comm2.commArray, {2, 1}));

  // commArray defined in the Potential
  pot.setCommArray({2, 1});
  CommGrid comm3;
  comm3.setup(pot, 12, {0,1});
  EXPECT_TRUE(ArraysMatch(comm3.commArray, {2, 1}));

  // Invalid commArray
  EXPECT_ANY_THROW({
    pot.setCommArray({1, 3});
    CommGrid().setup(pot, 12, {0,1});
  });
  EXPECT_ANY_THROW({
    pot.setCommArray({2});
    CommGrid().setup(pot, 12, {0,1});
  });
}


TEST(CommGrid, TestAssign) {
  if (mpi.rank >= 2) return;
  // Initialise communicator
  CommGrid comm;
  comm.commArray = {1, 2};
  GridPot pot({2,6});
  comm.setup(pot, 12, {0,1});

  // Test assignment of global data
  vector<double> globalData{10,11,12,13,14,15, 20,21,22,23,24,25};
  vector<double> procData;
  if (comm.rank() == 0) {
    procData = {15,10,11,12,13, 25,20,21,22,23};
  } else if (comm.rank() == 1) {
    procData = {12,13,14,15,10, 22,23,24,25,20};
  }
  EXPECT_TRUE(ArraysMatch(comm.assignProc(globalData), procData));
  EXPECT_TRUE(ArraysMatch(comm.assignBlock(globalData), procData));

  // Test assignment of block data (no halo)
  vector<double> blockData = {10,11,12, 20,21,22};
  procData = {0,10,11,12,0, 0,20,21,22,0};
  EXPECT_TRUE(ArraysMatch(comm.assignBlock(blockData), procData));
}


TEST(CommGrid, TestGet) {
  if (mpi.rank >= 2) return;
  // Initialise communicator
  CommGrid comm;
  comm.commArray = {1, 2};
  GridPot pot({2,6});
  comm.setup(pot, 12, {0,1});

  // Test getLocalIdx first
  EXPECT_EQ(comm.getLocalIdx(6, 0), (comm.rank()==0) ? 6 : -1);
  EXPECT_EQ(comm.getLocalIdx(10, 1), (comm.rank()==0) ? -1 : 7);
  EXPECT_EQ(comm.getLocalIdx(6), (comm.rank()==0) ? 6 : -1);
  EXPECT_EQ(comm.getLocalIdx(10), (comm.rank()==0) ? -1 : 7);

  // Get points from a distributed global index array
  vector<double> globalData{0,1,2,3,4,5, 6,7,8,9,10,11};
  vector<double> procData = comm.assignProc(globalData);
  EXPECT_EQ(comm.get(procData, 6), 6);
  EXPECT_EQ(comm.get(procData, 10), 10);
}


TEST(CommGrid, TestCommunicate) {
  // Initialise communicator
  CommGrid comm;
  comm.commArray = {2, 2};
  GridPot pot({4,4});
  comm.setup(pot, 12, {});

  vector<double> data, result;
  if (comm.rank() == 0) {
    data   = { 0,  0,  0,  0,
               0, 11, 12,  0,
               0, 13, 14,  0,
               0,  0,  0,  0};
    result = {44, 33, 34, 43,
              22, 11, 12, 21,
              24, 13, 14, 23,
              42, 31, 32, 41};
  } else if (comm.rank() == 1) {
    data   = { 0,  0,  0,  0,
               0, 21, 22,  0,
               0, 23, 24,  0,
               0,  0,  0,  0};
    result = {34, 43, 44, 33,
              12, 21, 22, 11,
              14, 23, 24, 13,
              32, 41, 42, 31};
  } else if (comm.rank() == 2) {
    data   = { 0,  0,  0,  0,
               0, 31, 32,  0,
               0, 33, 34,  0,
               0,  0,  0,  0};
    result = {24, 13, 14, 23,
              42, 31, 32, 41,
              44, 33, 34, 43,
              22, 11, 12, 21};
  } else if (comm.rank() == 3) {
    data   = { 0,  0,  0,  0,
               0, 41, 42,  0,
               0, 43, 44,  0,
               0,  0,  0,  0};
    result = {14, 23, 24, 13,
              32, 41, 42, 31,
              34, 43, 44, 33,
              12, 21, 22, 11};
  }
  comm.communicate(data);
  EXPECT_TRUE(ArraysMatch(data, result));
}


TEST(CommGrid, TestCommunicateAccumulate) {
  // Initialise communicator
  CommGrid comm;
  comm.commArray = {2, 2};
  GridPot pot({4,6});
  comm.setup(pot, 24, {});

  vector<double> data = {3,  1,  1,  1, 3,
                         2, 10, 20, 30, 2,
                         2, 40, 50, 60, 2,
                         3,  1,  1,  1, 3};
  vector<double> result = {3,  1,  1,  1, 3,
                           2, 16, 21, 36, 2,
                           2, 46, 51, 66, 2,
                           3,  1,  1,  1, 3};
  comm.communicateAccumulate(data);
  EXPECT_TRUE(ArraysMatch(data, result));
}


TEST(CommGrid, TestGather) {
  if (mpi.rank >= 2) return;
  // Initialise communicator
  CommGrid comm;
  comm.commArray = {1, 2};
  GridPot pot({2,4});
  comm.setup(pot, 8, {0,1});

  auto data = 10*(comm.rank()+1) + vector<double>{0,1,2,0, 0,3,4,0};
  EXPECT_TRUE(ArraysMatch(comm.gather(data), {11,12,21,22, 13,14,23,24}));
}


TEST(CommGrid, TestDotProduct) {
  CommGrid comm;
  comm.commArray = {2, 2};
  GridPot pot({4,6});
  comm.setup(pot, 24, {});

  vector<double> a = {1,  1,  1,  1, 1,
                      1, 10, 10, 10, 1,
                      1, 20, 20, 20, 1,
                      1,  1,  1,  1, 1};
  vector<double> b = {1,  1,  1,  1, 1,
                      1, 20, 20, 20, 1,
                      1, 10, 10, 10, 1,
                      1,  1,  1,  1, 1};
  EXPECT_EQ(comm.dotProduct(a, b), 4*1200);
}

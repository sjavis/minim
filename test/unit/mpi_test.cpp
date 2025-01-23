#include "test_main.cpp"
#include "utils/mpi.h"

using namespace minim;

TEST(MpiTest, sizeRank) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  EXPECT_EQ(mpi.size, 2);
  EXPECT_EQ(mpi.rank, rank);
}

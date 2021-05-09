#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

void mpiInit();
void mpiInit(int *argc, char ***argv);


class Communicator {
  public:
    int size;
    int rank;

    Communicator();
    ~Communicator();
    void init(int *argc, char ***argv);

  private:
    bool _init = false;
};

namespace minim {
  extern Communicator mpi;
}

#endif

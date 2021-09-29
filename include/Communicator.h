#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>

namespace minim {
  class Args;
  class Priv;

  class Communicator {
    private:
     typedef std::vector<double> Vector;
     Priv *priv;

    public:
      int ndof;
      int nproc;
      int nblock;

      Communicator(int ndof, Args &args);
      ~Communicator();

      Vector assignBlock(const Vector &in);

      void communicate(Vector &vector);
      double get(const Vector &vector, int i);
      double dotProduct(const Vector &a, const Vector &b);

      Vector gather(const Vector &block, int root=-1);
      Vector scatter(const Vector &data, int root=-1);

      void bcast(int &value, int root=0);
      void bcast(double &value, int root=0);
      void bcast(Vector &value, int root=0);
  };

}

#endif

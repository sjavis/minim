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

      Vector assignBlock(const Vector &in) const;

      void communicate(Vector &vector) const;
      double get(const Vector &vector, int i) const;
      double dotProduct(const Vector &a, const Vector &b) const;

      Vector gather(const Vector &block, int root=-1) const;
      Vector scatter(const Vector &data, int root=-1) const;

      void bcast(int &value, int root=0) const;
      void bcast(double &value, int root=0) const;
      void bcast(Vector &value, int root=0) const;
  };

}

#endif

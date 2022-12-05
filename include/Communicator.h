#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include <memory>
#include "Potential.h"

namespace minim {

  class Communicator {
    typedef std::vector<double> Vector;

    public:
      size_t ndof;   //!< Total number of degrees of freedom
      size_t nproc;  //!< Total number of degrees of freedom on processor (including halo)
      size_t nblock; //!< Number of degrees of freedom assigned to processor (excluding halo)
      int iblock;    //!< The starting index for this processor
      bool usesThisProc = true;

      Communicator(Potential& pot, size_t ndof, std::vector<int> ranks);
      Communicator(const Communicator& comm);
      Communicator& operator=(const Communicator& comm);
      ~Communicator();

      int rank() const;
      int size() const;

      Vector assignBlock(const Vector& in) const;
      Vector assignProc(const Vector& in) const;

      void communicate(Vector& vector) const;
      double get(const Vector& vector, int i) const;
      double dotProduct(const Vector& a, const Vector& b) const;

      Vector gather(const Vector& block, int root=-1) const;
      Vector scatter(const Vector& data, int root=-1) const;

      void bcast(int& value, int root=0) const;
      void bcast(double& value, int root=0) const;
      void bcast(Vector& value, int root=0) const;

    private:
      class Priv;
      std::unique_ptr<Priv> p;
  };

}

#endif

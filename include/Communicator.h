#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>
#include <memory>
#include "Potential.h"

namespace minim {
  using std::vector;

  class Communicator {
    public:
      size_t ndof;   //!< Total number of degrees of freedom
      size_t nproc;  //!< Total number of degrees of freedom on processor (including halo)
      size_t nblock; //!< Number of degrees of freedom assigned to processor (excluding halo)
      int iblock;    //!< The starting index for this processor

      bool usesThisProc = true;
      vector<int> ranks = vector<int>();

      Communicator();
      Communicator(const Communicator& comm);
      Communicator& operator=(const Communicator& comm);
      ~Communicator();
      void setup(Potential& pot, size_t ndof, vector<int> ranks);

      int rank() const;
      int size() const;

      vector<double> assignBlock(const vector<double>& in) const;
      vector<double> assignProc(const vector<double>& in) const;

      void communicate(vector<double>& vector) const;
      double get(const vector<double>& vector, int i) const;

      double sum(double a) const;
      double sum(const vector<double>& a) const;
      double dotProduct(const vector<double>& a, const vector<double>& b) const;

      vector<double> gather(const vector<double>& block, int root=-1) const;
      vector<double> scatter(const vector<double>& data, int root=-1) const;

      void bcast(int& value, int root=0) const;
      void bcast(double& value, int root=0) const;
      void bcast(vector<double>& value, int root=0) const;

    private:
      class Priv;
      std::unique_ptr<Priv> p;
  };

}

#endif

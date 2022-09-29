#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <memory>

namespace minim {
  class State;

  class Potential {
    public:
      class Args {
        public:
          struct Element {
            int id;
            int type;
            std::vector<int> idof;
            std::vector<double> parameters;
          };

          const int ndof;
          std::vector<Element> elements;
          std::vector<Element> elements_halo;

          Args(int ndof) : ndof(ndof) {};
          Args(int ndof, std::vector<std::vector<int>> idofs);
          Args(int ndof, std::vector<std::vector<int>> idofs, std::vector<int> types,
               std::vector<std::vector<double>> parameters);
          ~Args() {};
      };

    private:
      typedef std::vector<double> Vector;
      typedef double (*EFunc)(const Vector&, const Args&);
      typedef Vector (*GFunc)(const Vector&, const Args&);
      EFunc _energy;
      GFunc _gradient;

    public:
      Potential(EFunc energy, GFunc gradient) : _energy(energy), _gradient(gradient) {};
      ~Potential() {};

      virtual double energy(const Vector &coords, const Args &args) const;
      virtual Vector gradient(const Vector &coords, const Args &args) const;

      virtual std::unique_ptr<Args> newArgs(int ndof);
      State newState(Vector coords);
      State newState(int ndof);

    protected:
      Potential() {};
  };

}

#endif

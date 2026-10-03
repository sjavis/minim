#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>
#include <memory>
#include <functional>

namespace minim {
  class State;
  class Communicator;

  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;


  /// Abstract base class for system potentials.
  ///
  /// The potential computes the energy and gradient of a set of coordinates.
  /// It can be constructed from user-supplied functions, or derived from to
  /// create more complex potentials.
  class Potential {
    typedef double (*EFunc)(const vector<double>&);
    typedef vector<double> (*GFunc)(const vector<double>&);
    typedef void (*EGFunc)(const vector<double>&, double*, vector<double>*);
    EFunc _energy;
    GFunc _gradient;
    EGFunc _energyGradient;

    public:
      /// Convergence criterion used when creating a new state.
      double convergence = 1e-6;

      // Basic user functions

      /// Construct a potential from separate energy and gradient functions.
      Potential(EFunc energy, GFunc gradient);
      /// Construct a potential from a combined energy and gradient function.
      Potential(EGFunc energyGradient);

      /// Compute the energy of the given coordinates.
      virtual double energy(const vector<double>& coords) const;
      /// Compute the energy gradient of the given coordinates.
      virtual vector<double> gradient(const vector<double>& coords) const;
      /// Compute the energy and gradient of the given coordinates.
      ///
      /// @param coords The coordinates to compute the energy and gradient of.
      /// @param comm The communicator, used for distributed calculations.
      /// @param e Output: the total energy.
      /// @param g Output: the gradient of the total energy.
      virtual void energyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const;
      // Note: In any derived classes, either energyGradient or (energy and gradient) MUST be overridden

      /// Create a new state of the given number of degrees of freedom.
      ///
      /// @param ndof The number of degrees of freedom.
      /// @param ranks The MPI ranks to distribute the state over
      ///   (default: all ranks).
      State newState(int ndof, const vector<int>& ranks={});
      /// Create a new state from initial coordinates.
      ///
      /// @param coords The initial coordinates.
      /// @param ranks The MPI ranks to distribute the state over
      ///   (default: all ranks).
      State newState(const vector<double>& coords, const vector<int>& ranks={});

      // Constraints

      /// A constraint on a set of degrees of freedom.
      ///
      /// Constrains the given degrees of freedom to move only along a
      /// normal direction, optionally with a correction function.
      struct Constraint {
        using NormalFn = std::function<vector<double>(const vector<int>&, const vector<double>&)>;
        using CorrectionFn = std::function<void(const vector<int>&, vector<double>&)>;
        /// The degrees of freedom the constraint applies to.
        vector<int> idof;
        /// Function returning the constrained components of the normal vector.
        NormalFn normalFn;
        /// Optional function to correct the coordinates to satisfy the constraint.
        CorrectionFn correction = nullptr;
        vector<double> normal(const vector<double>& normalVec) const { return normalFn(idof, normalVec); }
      };
      /// The constraints applied to the state.
      vector<Constraint> constraints;

      /// Constrain the given degrees of freedom to be fixed.
      Potential& setConstraints(vector<int> iFix);
      /// Constrain the given degrees of freedom to a linear constraint.
      Potential& setConstraints(vector2d<int> idofs, vector<double> normal);
      /// Constrain the given degrees of freedom with a custom normal and correction function.
      Potential& setConstraints(vector2d<int> idofs, Constraint::NormalFn normal, Constraint::CorrectionFn correction=nullptr);

      /// Apply the constraints to a proposed step.
      void applyConstraints(const vector<double>& coords, const Communicator& comm, vector<double>& step) const;
      /// Correct the coordinates to satisfy the constraints.
      void correctConstraints(vector<double>& coords) const;
      virtual void specialisedConstraints(const vector<double>& coords, const Communicator& comm, vector<double>& data) const {};

      /// Whether the given degree of freedom is fixed.
      bool isFixed(int index) const;
      /// Whether each of the given degrees of freedom is fixed.
      vector<char> isFixed(const vector<int>& indicies) const;

      // Internal

      /// Initialise the potential with the coordinates of a new state.
      virtual void init(const vector<double>& coords) {};
      /// Initialise the potential with local (distributed) data.
      ///
      /// Take care using this, if the potential is cloned any distributed
      /// parameters will be copied as they are.
      virtual void initLocal(const vector<double>& coords, const Communicator& comm) {};
      void energyGradientWrapper(const vector<double>& coords, double* e, vector<double>* g, const Communicator& comm);
      bool isSerial() const;

      // Copy / destruct
      virtual ~Potential() = default;
      /// Create a deep copy of the potential.
      virtual std::unique_ptr<Potential> clone() const {
        return std::make_unique<Potential>(*this);
      }


      // Members and functions specific to different types of potential

      /// The type of parallelisation used by the potential.
      enum{
        SERIAL = 0,
        UNSTRUCTURED = 1,
        GRID = 2,
      };
      virtual int potentialType() const { return SERIAL; };
      virtual std::unique_ptr<Communicator> newComm() const;

      // UNSTRUCTURED: Energy elements for parallelisation

      /// Whether the potential's elements are distributed over processors.
      bool distributed = false;

      /// A single energy element, used for parallelisation of unstructured potentials.
      struct Element {
        /// The element type, used to select the energy calculation.
        int type;
        /// The degrees of freedom in the element.
        vector<int> idof;
        /// Any parameters of the element.
        vector<double> parameters;
      };
      /// The energy elements used for parallelisation.
      vector<Element> elements;
      vector<Element> elements_halo;
      /// Set the elements for an unstructured potential.
      Potential& setElements(vector<Element> elements);
      /// Set the elements for an unstructured potential from degrees of freedom.
      Potential& setElements(vector2d<int> idofs);
      /// Set the elements for an unstructured potential from degrees of freedom, types, and parameters.
      Potential& setElements(vector2d<int> idofs, vector<int> types, vector2d<double> parameters);

      /// Compute the energy and gradient of a single element.
      virtual void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {};
      /// Compute the energy and gradient of the local block.
      virtual void blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const {};

      // GRID

      /// Number of degrees of freedom per grid node.
      int dofPerNode = 1;
      /// Width of the halo region for grid potentials.
      int haloWidth = 1;
      /// Size of the grid along each dimension.
      vector<int> gridSize;
      /// The number of processors along each grid dimension.
      vector<int> commArray;
      /// Set the number of processors along each grid dimension.
      Potential& setCommArray(vector<int> commArray);

    protected:
      Potential() : _energy(nullptr), _gradient(nullptr), _energyGradient(nullptr) {};
  };


  /// An intermediate class used to return the derived type for methods that return a Potential.
  template<typename Derived>
  class NewPotential : public Potential {
    public:
      std::unique_ptr<Potential> clone() const override {
        return std::make_unique<Derived>(static_cast<const Derived&>(*this));
      }

      /// Set the elements for an unstructured potential.
      Derived& setElements(vector<Element> elements) {
        return static_cast<Derived&>(Potential::setElements(elements));
      }
      /// Set the elements for an unstructured potential from degrees of freedom.
      Derived& setElements(vector2d<int> idofs) {
        return static_cast<Derived&>(Potential::setElements(idofs));
      }
      /// Set the elements for an unstructured potential from degrees of freedom, types, and parameters.
      Derived& setElements(vector2d<int> idofs, vector<int> types, vector2d<double> parameters) {
        return static_cast<Derived&>(Potential::setElements(idofs, types, parameters));
      }

      /// Constrain the given degrees of freedom to be fixed.
      Derived& setConstraints(vector<int> iFix) {
        return static_cast<Derived&>(Potential::setConstraints(iFix));
      }
      /// Constrain the given degrees of freedom to a linear constraint.
      Derived& setConstraints(vector2d<int> idofs, vector<double> normal) {
        return static_cast<Derived&>(Potential::setConstraints(idofs, normal));
      }
      /// Constrain the given degrees of freedom with a custom normal and correction function.
      Derived& setConstraints(vector2d<int> idofs, Constraint::NormalFn normal, Constraint::CorrectionFn correction=nullptr) {
        return static_cast<Derived&>(Potential::setConstraints(idofs, normal, correction));
      }
  };

}

#endif

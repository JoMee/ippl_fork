#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H 

#include <mpi.h>
template <typename T, unsigned Dim>
struct SimulationParameters {

    SimulationParameters(unsigned nt_, Vector_t<int, Dim>& nr_, std::string& solver_, double lbt_,
                         Vector_t<T, Dim> rmin_ = 0.0, Vector_t<T, Dim> rmax_ = 10.0,
                         Vector_t<T, Dim> origin_ = 0.0)
        : nt(nt_)
        , nr(nr_)
        , solver(solver_)
        , lbt(lbt_)
        , rmin(rmin_)
        , rmax(rmax_)
        , origin(origin_) {

      this->dr = this->rmax - this->rmin;

      this->hr = this->dr / this->nr;

      this->dt = std::min(0.05, 0.5 * ( *std::min_element(this->hr.begin(), this->hr.end()) ) );

      this->time = 0.0;

      this->it = 0;

      for (unsigned i = 0; i < Dim; i++) {
          this->domain[i] = ippl::Index(this->nr[i]);
      }
      
      this->decomp.fill(true);

      this->mesh = Mesh_t<Dim>(this->domain, this->hr, this->origin);
      this->fl = FieldLayout_t<Dim>(MPI_COMM_WORLD, this->domain, this->decomp, true);
    }

    int nt;
    Vector_t<int, Dim> nr;
    std::string solver;
    double lbt;

    double dt;
    double time;
    int it;
    Vector_t<T, Dim> hr;
    Vector_t<double, Dim> dr;
    ippl::NDIndex<Dim> domain;
    std::array<bool, Dim> decomp;
    Mesh_t<Dim> mesh;
    FieldLayout_t<Dim> fl;

    Vector_t<T, Dim> rmin;
    Vector_t<T, Dim> rmax;
    Vector_t<T, Dim> origin;
};

#endif

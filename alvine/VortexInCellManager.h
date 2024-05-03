#ifndef IPPL_VORTEX_IN_CELL_MANAGER_H
#define IPPL_VORTEX_IN_CELL_MANAGER_H

#include <memory>

#include "AlvineManager.h"
#include "FieldContainer.hpp"
#include "FieldSolver.hpp"
#include "LoadBalancer.hpp"
#include "ParticleFieldStrategy.hpp"
#include "Manager/BaseManager.h"
#include "Particle/ParticleBase.h"
#include "ParticleContainer.hpp"
#include "ParticleDistributions.h"
#include "Random/Distribution.h"
#include "Random/InverseTransformSampling.h"
#include "Random/NormalDistribution.h"
#include "Random/Randu.h"
#include "VortexDistributions.h"

template <typename T, unsigned Dim, typename ParticleDistribution, typename VortexDistribution>
class VortexInCellManager : public AlvineManager<T, Dim> {
public:
    using ParticleContainer_t = ParticleContainer<T, Dim>;
    using FieldContainer_t    = FieldContainer<T, Dim>;
    //using LoadBalancer_t      = LoadBalancer<T, Dim>;
    bool remove_particles;

    VortexInCellManager(unsigned nt_, Vector_t<int, Dim>& nr_, std::string& solver_, double lbt_,
        Vector_t<double, Dim> rmin_ = 0.0,
        Vector_t<double, Dim> rmax_ = 10.0,
        Vector_t<double, Dim> origin_ = 0.0,
        bool remove_particles = false)
        : AlvineManager<T, Dim>(nt_, nr_, solver_, lbt_) {
            this->rmin_m = rmin_;
            this->rmax_m = rmax_;
            this->origin_m = origin_;
            this->remove_particles = remove_particles;
        }

    ~VortexInCellManager() {}

    void pre_run() override {

      Inform csvout(NULL, "particles.csv", Inform::OVERWRITE);
      csvout.precision(16);
      csvout.setf(std::ios::scientific, std::ios::floatfield);

      if constexpr (Dim == 2) {
          csvout << "time,index,pos_x,pos_y,vorticity" << endl;
      } else {
          csvout << "time,index,pos_x,pos_y,pos_z" << endl;
      }

      for (unsigned i = 0; i < Dim; i++) {
          this->domain_m[i] = ippl::Index(this->nr_m[i]);
      }

      Vector_t<double, Dim> dr = this->rmax_m - this->rmin_m;

      this->hr_m = dr / this->nr_m;

      // Courant condition
      this->dt_m = std::min(0.05, 0.5 * ( *std::min_element(this->hr_m.begin(), this->hr_m.end()) ) );

      this->it_m = 0;
      this->time_m = 0.0;

      int density = 1; // particles per cell
      set_number_of_particles(density);

      this->decomp_m.fill(true);
      this->isAllPeriodic_m = true;

      this->setFieldContainer(std::make_shared<FieldContainer_t>(
            this->hr_m, this->rmin_m, this->rmax_m, this->decomp_m, this->domain_m, this->origin_m,
            this->isAllPeriodic_m));

      
      if constexpr (Dim == 2) {
          std::shared_ptr<FieldContainer<T, 2>> fc = std::dynamic_pointer_cast<FieldContainer<T, 2>>(this->fcontainer_m);
          this->setParticleContainer(std::make_shared<TwoDimParticleContainer<T>>(fc->getMesh(), fc->getFL()));
          this->setFieldSolver( std::make_shared<TwoDimFFTSolverStrategy<T>>() );
          this->setParticleFieldStrategy( std::make_shared<TwoDimParticleFieldStrategy<T>>() );

      } else if constexpr (Dim == 3) {
          this->setParticleFieldStrategy( std::make_shared<ThreeDimParticleFieldStrategy<T>>() );
      }

      this->fsolver_m->initSolver(this->fcontainer_m);

      //this->setLoadBalancer( std::make_shared<LoadBalancer_t>( this->lbt_m, this->fcontainer_m, this->pcontainer_m, this->fsolver_m) );

      ParticleDistribution particle_dist(this->rmin_m, this->rmax_m, this->np_m);
      VortexDistribution vortex_dist(this->rmin_m, this->rmax_m, this->origin_m);


      this->pcontainer_m->initializeParticles(particle_dist, vortex_dist);
      std::cout << "Particles initialized" << std::endl;

      // if (this->remove_particles) {
      //     int removed = pc->removeInvalidParticles();
      //     this->np_m -= removed;
      // }

      this->par2grid();

      this->fsolver_m->solve(this->fcontainer_m);
      this->updateFields();

      this->grid2par();

      std::shared_ptr<ParticleContainer<T, Dim>> pc = std::dynamic_pointer_cast<ParticleContainer<T, Dim>>(this->pcontainer_m);

      pc->R_old = pc->R;
      pc->R = pc->R_old + pc->P * this->dt_m;
      pc->update();
    }

    void set_number_of_particles(int density) {
        int particles = 1;
        for (unsigned i = 0; i < Dim; i++) {
            particles *= this->nr_m[i];
        }
        this->np_m = particles*density;
        std::cout << "Number of particles: " << this->np_m << std::endl;
    }

    void advance() override {
      LeapFrogStep();     
    }

    void LeapFrogStep() {

      static IpplTimings::TimerRef PTimer           = IpplTimings::getTimer("pushVelocity");
      static IpplTimings::TimerRef RTimer           = IpplTimings::getTimer("pushPosition");
      static IpplTimings::TimerRef updateTimer      = IpplTimings::getTimer("update");
      static IpplTimings::TimerRef SolveTimer       = IpplTimings::getTimer("solve");
      
      std::shared_ptr<ParticleContainer<T, Dim>> pc = std::dynamic_pointer_cast<ParticleContainer<T, Dim>>(this->pcontainer_m);

      // scatter the vorticity to the underlying grid
      this->par2grid();

      // claculate stream function
      IpplTimings::startTimer(SolveTimer);
      this->fsolver_m->solve(this->fcontainer_m);
      IpplTimings::stopTimer(SolveTimer);

      // calculate velocity from stream function
      IpplTimings::startTimer(PTimer);
      this->updateFields();
      IpplTimings::stopTimer(PTimer);

      // gather velocity field
      this->grid2par();

      //drift
      IpplTimings::startTimer(RTimer);
      typename ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>::particle_position_type R_old_temp = pc->R_old;

      pc->R_old = pc->R;
      pc->R = R_old_temp + 2 * pc->P * this->dt_m;
      IpplTimings::stopTimer(RTimer);


      IpplTimings::startTimer(updateTimer);
      pc->update();
      IpplTimings::stopTimer(updateTimer);

    }

    void dump() override {
      std::shared_ptr<TwoDimParticleContainer<T>> pc = std::dynamic_pointer_cast<TwoDimParticleContainer<T>>(this->pcontainer_m);

      Inform csvout(NULL, "particles.csv", Inform::APPEND);

      for (unsigned i = 0; i < this->np_m; i++) {
        csvout << this->it_m << "," << i; 
        for (unsigned d = 0; d < Dim; d++) {
          csvout << "," << pc->R(i)[d];
        }
        csvout << "," << pc->omega(i) << endl;
      }
       
    }
};
#endif

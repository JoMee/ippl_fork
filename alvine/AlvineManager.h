#ifndef IPPL_ALVINE_MANAGER_H
#define IPPL_ALVINE_MANAGER_H

#include <memory>

#include "FieldContainer.hpp"
#include "FieldSolver.hpp"
#include "LoadBalancer.hpp"
#include "Manager/BaseManager.h"
#include "Manager/FieldSolverBase.h"
#include "Manager/PicManager.h"
#include "ParticleContainer.hpp"
#include "ParticleFieldStrategy.hpp"
#include "Random/Distribution.h"
#include "Random/InverseTransformSampling.h"
#include "Random/NormalDistribution.h"
#include "Random/Randn.h"

using view_type = typename ippl::detail::ViewType<ippl::Vector<double, Dim>, 1>::view_type;

template <typename T, unsigned Dim>
class AlvineManager
    : public ippl::PicManager<T, Dim, ParticleContainerBase, FieldContainerBase,
                              LoadBalanceStrategy, FieldSolverStrategy<FieldContainerBase>> {
public:
    using particle_field_strategy_type =
        typename std::shared_ptr<ParticleFieldStrategy<FieldContainerBase, ParticleContainerBase>>;

protected:
    unsigned nt_m;
    unsigned it_m;
    unsigned np_m;
    Vector_t<int, Dim> nr_m;
    std::array<bool, Dim> decomp_m;
    bool isAllPeriodic_m;
    ippl::NDIndex<Dim> domain_m;
    std::string solver_m;
    double lbt_m;
    particle_field_strategy_type particle_field_strategy_m;

public:
    AlvineManager(unsigned nt_, Vector_t<int, Dim>& nr_, std::string& solver_, double lbt_)
        : ippl::PicManager<T, Dim, ParticleContainerBase, FieldContainerBase, LoadBalanceStrategy,
                           FieldSolverStrategy<FieldContainerBase>>()
        , nt_m(nt_)
        , nr_m(nr_)
        , solver_m(solver_)
        , lbt_m(lbt_)
        , particle_field_strategy_m(nullptr) {}

    ~AlvineManager() {}

protected:
    double time_m;
    double dt_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;
    Vector_t<double, Dim> origin_m;
    Vector_t<double, Dim> hr_m;

public:
    void setParticleFieldStrategy(particle_field_strategy_type particle_field_strategy) {
        particle_field_strategy_m = particle_field_strategy;
    }

    double getTime() { return time_m; }

    void setTime(double time_) { time_m = time_; }

    int getNt() const { return nt_m; }

    void setNt(int nt_) { nt_m = nt_; }

    virtual void dump(){/* default does nothing */};

    void pre_step() override {
        Inform m("Pre-step");
        m << "Done" << endl;
    }

    void post_step() override {
        this->time_m += this->dt_m;
        this->it_m++;

        this->dump();
    }

    void grid2par() override {
        if (particle_field_strategy_m) {
            particle_field_strategy_m->grid2par(this->fcontainer_m, this->pcontainer_m);
        } else {
            throw std::runtime_error("Particle-Field strategy not defined");
        }
    }

    void par2grid() override {
        if (particle_field_strategy_m) {
            particle_field_strategy_m->par2grid(this->fcontainer_m, this->pcontainer_m);
        } else {
            throw std::runtime_error("Particle-Field strategy not defined");
        }
    }

    void updateFields() {
        if (particle_field_strategy_m) {
            particle_field_strategy_m->updateFields(this->fcontainer_m);
        } else {
            throw std::runtime_error("Particle-Field strategy not defined");
        }
    }
                                
    void scatterCIC() {
      this->fcontainer_m->getOmegaField() = 0.0;
      if constexpr (Dim == 2) {
          scatter(this->pcontainer_m->omega, this->fcontainer_m->getOmegaField(), this->pcontainer_m->R);
      } else if constexpr (Dim == 3) {
        //TODO: for some reason the scatter method doesn't work in three dimensions but gather does. 
      }
    }

    void initDumpEnergy() {
        Inform energyout(NULL, "energy.csv", Inform::OVERWRITE);
        energyout.precision(16);
        energyout.setf(std::ios::scientific, std::ios::floatfield);
        energyout << "time,energy" << endl;
    }

    void dumpEnergy() {
        static IpplTimings::TimerRef ETimer = IpplTimings::getTimer("energy");
        
        std::shared_ptr<ParticleContainer<T, Dim>> pc = std::dynamic_pointer_cast<ParticleContainer<T, Dim>>(this->pcontainer_m);
        double energy = 0.0;

        IpplTimings::startTimer(ETimer);

        Kokkos::parallel_reduce(
            "Compute energy", this->np_m,
            KOKKOS_LAMBDA(const int i, double& local_sum) {
                for (unsigned d = 0; d < Dim; d++) {
                    local_sum += pc->P(i)[d] * pc->P(i)[d];
                }
            },
            energy);

        IpplTimings::stopTimer(ETimer);

        Inform energyout(NULL, "energy.csv", Inform::APPEND);
        energyout << this->it_m << "," << energy << endl;
    }
};
#endif

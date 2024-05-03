#ifndef IPPL_PARTICLE_CONTAINER_H
#define IPPL_PARTICLE_CONTAINER_H

#include <memory>

#include "Manager/BaseManager.h"
#include "ParticleDistributions.h"
#include "VortexDistributions.h"


class ParticleContainerBase {
public:
    virtual ~ParticleContainerBase() = default;
    
    virtual int removeInvalidParticles() =0;

    virtual void initializeParticles(BaseParticleDistribution& vortex_dist, BaseVortexDistribution& particle_dist) =0;
};

using view_type     = typename ippl::detail::ViewType<ippl::Vector<double, Dim>, 1>::view_type;
using host_type     = typename ippl::ParticleAttrib<T>::HostMirror;
using vector_type   = Vector_t<double, Dim>;

// Define the ParticlesContainer class
template <typename T, unsigned Dim>
class ParticleContainer : public ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>,
                          public ParticleContainerBase {
    using Base       = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;
    using omega_type = std::conditional<Dim == 2, ippl::ParticleAttrib<T>,
                                        typename Base::particle_position_type>::type;
    using valid_type = ippl::ParticleAttrib<bool>;

public:
    typename Base::particle_position_type R_old;
    typename Base::particle_position_type P;
    omega_type omega;
    valid_type invalid;

private:
    PLayout_t<T, Dim> pl_m;

public:
    ParticleContainer(Mesh_t<Dim>& mesh, FieldLayout_t<Dim>& FL)
        : pl_m(FL, mesh) {
        this->initialize(pl_m);
        registerAttributes();
        setupBCs();
    }

    virtual ~ParticleContainer() = default;

    std::shared_ptr<PLayout_t<T, Dim>> getPL() { return pl_m; }
    void setPL(std::shared_ptr<PLayout_t<T, Dim>>& pl) { pl_m = pl; }

    void setupBCs() { setBCAllPeriodic(); }

private:
    void setBCAllPeriodic() { this->setParticleBC(ippl::BC::PERIODIC); }
    void registerAttributes() {
        this->addAttribute(P);
        this->addAttribute(R_old);
        this->addAttribute(omega);
        this->addAttribute(invalid);
    }
};

template <typename T>
class TwoDimParticleContainer : public ParticleContainer<T, 2> {
public:
    TwoDimParticleContainer(Mesh_t<2>& mesh, FieldLayout_t<2>& FL)
        : ParticleContainer<T, 2>(mesh, FL) {}

    ~TwoDimParticleContainer() {}

    void initializeParticles(BaseParticleDistribution& particle_dist,  BaseVortexDistribution& vortex_dist) override {
        size_t np = particle_dist.np_m;

        // Create np_m particles in container
        this->create(np);  // TODO: local number of particles? from kokkos?

        // Assign positions
        view_type R = this->R.getView();  // Position vector
        particle_dist.setR(R);
        Kokkos::parallel_for(np, particle_dist);
        this->R.print();

        // Assign vorticity
        host_type omega_host = this->omega.getHostMirror();  // Vorticity values
        vortex_dist.setR(R);
        vortex_dist.setOmega(omega_host);
        Kokkos::parallel_for(np, vortex_dist);
        Kokkos::deep_copy(this->omega.getView(), vortex_dist.omega);

        Kokkos::fence();
        ippl::Comm->barrier();
    }

    int removeInvalidParticles() override {

        int totalP = this->omega.getView().size();
        std::cout << "Total particles: " << totalP << std::endl;

        int total_invalid = 0;
            Kokkos::parallel_for(
                "Mark vorticity null as invalid", totalP, KOKKOS_LAMBDA(const size_t i) {
                    this->invalid.getView()(i) = false;
                    if (this->omega.getView()(i) == 0) {
                        this->invalid.getView()(i) = true;
                    }
                });

            for (int i = 0; i < totalP; i++) {
                this->invalid.getView()(i) ? total_invalid++ : total_invalid;
            }

            if (total_invalid and (total_invalid < int(totalP))) {
                std::cout << "Removing " << total_invalid << " particles" << std::endl;
                const auto invalid = this->invalid.getView();
                this->destroy(invalid, total_invalid);
            } else {
                std::cout << "No particles removed" << std::endl;
            }
        return total_invalid;
    }
};

#endif

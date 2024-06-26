#ifndef IPPL_LOAD_BALANCER_H
#define IPPL_LOAD_BALANCER_H

#include "ParticleContainer.hpp"
#include <memory>

class LoadBalanceStrategy {

    virtual ~LoadBalanceStrategy() = default;
};

template <typename T, unsigned Dim>
class LoadBalancer{
    using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;
    using FieldSolver_t = FieldSolverStrategy<FieldContainer<T, 2>>;
    using vorticity_field_type = std::conditional<Dim == 2, Field<T, Dim>, VField_t<T, Dim>>::type;
    
    private:
        double loadbalancethreshold_m;
        vorticity_field_type* omega_m;
        std::shared_ptr<ParticleContainer<T, Dim>> pc_m;
        std::shared_ptr<FieldSolver_t> fs_m;
        unsigned int loadbalancefreq_m;
        ORB<T, Dim> orb;
    public:
        LoadBalancer(double lbs, std::shared_ptr<FieldContainer<T,Dim>> &fc, std::shared_ptr<ParticleContainer<T, Dim>> &pc, std::shared_ptr<FieldSolver_t> &fs)
           :loadbalancethreshold_m(lbs), omega_m(&fc->getOmegaField()), pc_m(pc), fs_m(fs) {}

        ~LoadBalancer() {  }

        double getLoadBalanceThreshold() const { return loadbalancethreshold_m; }
        void setLoadBalanceThreshold(double threshold) { loadbalancethreshold_m = threshold; }

        Field_t<Dim>* getOmega() const { return omega_m; }
        void setOmega(Field_t<Dim>* omega) { omega_m = omega; }

        std::shared_ptr<ParticleContainer<T, Dim>> getParticleContainer() const { return pc_m; }
        void setParticleContainer(std::shared_ptr<ParticleContainer<T, Dim>> pc) { pc_m = pc; }

        std::shared_ptr<FieldSolver_t> getFieldSolver() const { return fs_m; }
        void setFieldSolver(std::shared_ptr<FieldSolver_t> fs) { fs_m = fs; }

        void updateLayout(ippl::FieldLayout<Dim>* fl, ippl::UniformCartesian<T, Dim>* mesh, bool& isFirstRepartition) {
            // Update local fields

            static IpplTimings::TimerRef tupdateLayout = IpplTimings::getTimer("updateLayout");
            IpplTimings::startTimer(tupdateLayout);
            (*omega_m).updateLayout(*fl);

            // Update layout with new FieldLayout
            PLayout_t<T, Dim>* layout = &pc_m->getLayout();
            (*layout).updateLayout(*fl, *mesh);
            IpplTimings::stopTimer(tupdateLayout);
            static IpplTimings::TimerRef tupdatePLayout = IpplTimings::getTimer("updatePB");
            IpplTimings::startTimer(tupdatePLayout);
            if (!isFirstRepartition) {
                pc_m->update();
            }
            IpplTimings::stopTimer(tupdatePLayout);
        }

        void initializeORB(ippl::FieldLayout<Dim>* fl, ippl::UniformCartesian<T, Dim>* mesh) {
            orb.initialize(*fl, *mesh, *omega_m);
        }

        void repartition(ippl::FieldLayout<Dim>* fl, ippl::UniformCartesian<T, Dim>* mesh, bool& isFirstRepartition) {
            // Repartition the domains

            using Base = ippl::ParticleBase<ippl::ParticleSpatialLayout<T, Dim>>;
            typename Base::particle_position_type *R;
            R = &pc_m->R;
            bool res = orb.binaryRepartition(*R, *fl, isFirstRepartition);
            if (res != true) {
                std::cout << "Could not repartition!" << std::endl;
                return;
            }
            // Update
            this->updateLayout(fl, mesh, isFirstRepartition);
            if constexpr (Dim == 2 || Dim == 3) {
                if (fs_m->getStype() == "FFT") {
                    std::get<FFTSolver_t<T, Dim>>(fs_m->getSolver()).setRhs(*omega_m);
                }
            }
        }

        bool balance(size_type totalP, const unsigned int nstep) {
            if (ippl::Comm->size() < 2) {
                return false;
            }
            if (std::strcmp(TestName, "UniformPlasmaTest") == 0) {
                return (nstep % loadbalancefreq_m == 0);
            } else {
                int local = 0;
                std::vector<int> res(ippl::Comm->size());
                double equalPart = (double)totalP / ippl::Comm->size();
                double dev       = std::abs((double)pc_m->getLocalNum() - equalPart) / totalP;
                if (dev > loadbalancethreshold_m) {
                   local = 1;
                }
                MPI_Allgather(&local, 1, MPI_INT, res.data(), 1, MPI_INT,
                          ippl::Comm->getCommunicator());

               for (unsigned int i = 0; i < res.size(); i++) {
                   if (res[i] == 1) {
                        return true;
                   }
               }
               return false;
            }
        }
};

#endif

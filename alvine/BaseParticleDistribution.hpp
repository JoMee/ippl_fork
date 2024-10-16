#ifndef IPPL_PARTICLE_DISTRIBUTION_H
#define IPPL_PARTICLE_DISTRIBUTION_H
#include "Types/Vector.h"

#include "BaseDistributionFunction.hpp"

template <typename T, unsigned Dim>
class IParticleDistribution {
public:
    virtual ~IParticleDistribution()                                                      = default;
    virtual void generateDistribution()                                                   = 0;
    virtual int getNumParticles() const                                                   = 0;
    virtual void applyFilter(std::function<bool(const ippl::Vector<T, Dim>&)> filterFunc) = 0;
    virtual void merge(const IParticleDistribution<T, Dim>& other)                        = 0;

    virtual const typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type&
    getParticles() const = 0;
};

template <typename T, unsigned Dim>
class PlacementStrategy {
public:
    virtual ~PlacementStrategy() = default;
    virtual void placeParticles(
        typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type& container,
        ippl::Vector<T, Dim> rmin, ippl::Vector<T, Dim> rmax) const = 0;
};

template <typename T, unsigned Dim>
class ParticleDistributionBase : public IParticleDistribution<T, Dim> {
protected:
    ippl::Vector<T, Dim> rmin, rmax;
    typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type particle_container;
    std::unique_ptr<PlacementStrategy<T, Dim>> placementStrategy;

public:
    ParticleDistributionBase(ippl::Vector<T, Dim> rmin_, ippl::Vector<T, Dim> rmax_,
                             PlacementStrategy<T, Dim>* strategy)
        : rmin(rmin_)
        , rmax(rmax_)
        , placementStrategy(strategy) {}

    virtual void generateDistribution() override {
        if (placementStrategy) {
            placementStrategy->placeParticles(particle_container, rmin, rmax);
        } else {
            throw std::logic_error("Placement strategy not initialized.");
        }
    }

    void applyFilter(std::function<bool(const ippl::Vector<T, Dim>&)> filterFunc) override {
        Kokkos::View<int*> indices("indices", particle_container.extent(0));
        int numInside = 0;

        Kokkos::parallel_scan(
            "FilterParticles", particle_container.extent(0),
            KOKKOS_LAMBDA(const int& i, int& update, const bool final) {
                if (filterFunc(particle_container(i))) {
                    if (final) {
                        indices(update) = i;
                    }
                    update += 1;
                }
            },
            numInside);

        Kokkos::fence();

        typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type new_particles(
            "FilteredParticles", numInside);

        Kokkos::parallel_for(
            "CopyFilteredParticles", numInside,
            KOKKOS_LAMBDA(const int& i) { new_particles(i) = particle_container(indices(i)); });

        Kokkos::fence();
        particle_container = new_particles;
    }

    virtual void merge(const IParticleDistribution<T, Dim>& other) override {
        auto& otherParticles = other.getParticles();
        size_t originalSize  = particle_container.extent(0);
        size_t newSize       = originalSize + otherParticles.extent(0);

        typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type newContainer(
            "new_particles", newSize);

        Kokkos::parallel_for(
            "CopyOriginalParticles", originalSize,
            KOKKOS_LAMBDA(const int i) { newContainer(i) = particle_container(i); });

        Kokkos::parallel_for(
            "AppendNewParticles", otherParticles.extent(0),
            KOKKOS_LAMBDA(const int i) { newContainer(originalSize + i) = otherParticles(i); });

        Kokkos::fence();

        particle_container = newContainer;
    }

    ParticleDistributionBase& operator+=(const ParticleDistributionBase& other) {
        this->merge(other);
        return *this;
    }

    void setPlacementStrategy(PlacementStrategy<T, Dim>* strategy) {
        placementStrategy.reset(strategy);
    }

    int getNumParticles() const override { return particle_container.extent(0); }

    const typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type& getParticles()
        const override {
        return particle_container;
    }
};

template <typename T, unsigned Dim>
class GridPlacement : public PlacementStrategy<T, Dim> {
    ippl::Vector<int, Dim> num_points;
    ippl::Vector<T, Dim> rmin;
    ippl::Vector<T, Dim> rmax;

public:
    GridPlacement(ippl::Vector<int, Dim> num_points_)
        : num_points(num_points_) {}

    void placeParticles(
        typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type& container,
        ippl::Vector<T, Dim> rmin, ippl::Vector<T, Dim> rmax) const override {
        ippl::Vector<T, Dim> dr = (rmax - rmin) / num_points;

        size_t total_num =
            std::reduce(num_points.begin(), num_points.end(), 1, std::multiplies<int>());
        Kokkos::resize(container, total_num);

        if constexpr (Dim == 2) {
            Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0, 0}, {num_points(0), num_points(1)});

            Kokkos::parallel_for(
                "2DGridInit", policy, KOKKOS_LAMBDA(const int i, const int j) {
                    ippl::Vector<int, 2> loc(i, j);
                    container(i * num_points(0) + j) = rmin + dr * loc;
                });
        } else if constexpr (Dim == 3) {
            Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy(
                {0, 0, 0}, {num_points(0), num_points(1), num_points(2)});

            Kokkos::parallel_for(
                "3DLoopInit", policy, KOKKOS_LAMBDA(const int i, const int j, const int k) {
                    ippl::Vector<int, 3> loc(i, j, k);
                    this->particle_container(i * num_points(0) + j * num_points(1) * num_points(0)
                                             + k) = dr * loc;
                });
        }
        Kokkos::fence();
    }
};

template <typename T, unsigned Dim>
class RandomPlacement : public PlacementStrategy<T, Dim> {
    ippl::Vector<int, Dim> num_points;
    ippl::Vector<T, Dim> rmin;
    ippl::Vector<T, Dim> rmax;

    int seed = 42;
    GeneratorPool rand_pool = GeneratorPool((size_type)(seed + 100 * ippl::Comm->rank()));

public:
    RandomPlacement(ippl::Vector<int, Dim> num_points_)
        : num_points(num_points_) {}

    void placeParticles(
        typename ippl::detail::ViewType<ippl::Vector<T, Dim>, 1>::view_type& container,
        ippl::Vector<T, Dim> rmin, ippl::Vector<T, Dim> rmax) const override {

        size_t total_num =
            std::reduce(num_points.begin(), num_points.end(), 1, std::multiplies<int>());
        Kokkos::resize(container, total_num);

        if constexpr (Dim == 2) {
            Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0, 0}, {num_points(0), num_points(1)});

            Kokkos::parallel_for(
                "2DRandInit", policy, KOKKOS_LAMBDA(const int i, const int j) {
                    typename GeneratorPool::generator_type rand_gen = rand_pool.get_state();
                    ippl::Vector<int, 2> loc;

                    for (unsigned d = 0; d < 2; ++d) {
                        loc[d] = rand_gen.drand(rmin[d], rmax[d]);
                    }
                    
                    container(i * num_points(0) + j) = loc;
                    rand_pool.free_state(rand_gen);
                });

        } else if constexpr (Dim == 3) {
            Kokkos::MDRangePolicy<Kokkos::Rank<3>> policy(
                {0, 0, 0}, {num_points(0), num_points(1), num_points(2)});

            Kokkos::parallel_for(
                "3DRandInit", policy, KOKKOS_LAMBDA(const int i, const int j, const int k) {
                    typename GeneratorPool::generator_type rand_gen = rand_pool.get_state();
                    ippl::Vector<int, 3> loc;

                    for (unsigned d = 0; d < 3; ++d) {
                        loc[d] = rand_gen.drand(rmin[d], rmax[d]);
                    }

                    this->particle_container(i * num_points(0) + j * num_points(1) * num_points(0)
                                             + k) = loc;
                    rand_pool.free_state(rand_gen);
                });
        }
        Kokkos::fence();
    }
};

template <typename T, unsigned Dim>
class GridDistribution : public ParticleDistributionBase<T, Dim> {
public:
    GridDistribution(ippl::Vector<int, Dim> num_points, ippl::Vector<T, Dim> rmin_,
                     ippl::Vector<T, Dim> rmax_)
        : ParticleDistributionBase<T, Dim>(rmin_, rmax_, new GridPlacement<T, Dim>(num_points)) {
        this->generateDistribution();
    }
};

template <typename T, unsigned Dim>
class FilteredDistribution : public ParticleDistributionBase<T, Dim> {
private:
    const AbstractDistributionFunction<ippl::Vector<T, Dim>, T>& distFunction;

public:
    FilteredDistribution(const AbstractDistributionFunction<ippl::Vector<T, Dim>, T>& distFunction_,
                         ippl::Vector<T, Dim> rmin, ippl::Vector<T, Dim> rmax,
                         PlacementStrategy<T, Dim>* strategy, double threshold = 1e-8)
        : ParticleDistributionBase<T, Dim>(rmin, rmax, strategy)
        , distFunction(distFunction_) {
        this->generateDistribution();

        this->applyFilter([this, threshold](const ippl::Vector<T, Dim>& point) -> bool {
            return (this->distFunction.evaluate(point) > threshold);
        });
    }
};

#endif

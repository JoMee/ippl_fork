#ifndef IPPL_PLACEMENT_CONFIGURATIONS_H
#define IPPL_PLACEMENT_CONFIGURATIONS_H

#include "BaseDistributionFunction.hpp"
#include "BaseParticleDistribution.hpp"

class InitialPlacement {
    // Wrapper but is just to illustrate how to combine and add the distributions

public:
    using vector_type             = ippl::Vector<T, Dim>;
    using particle_container_type = std::shared_ptr<ParticleContainer<T, Dim>>;
    using field_continer_type     = std::shared_ptr<FieldContainer<T, Dim>>;

    using view_type          = ippl::detail::ViewType<vector_type, 1>::view_type;
    using dist_function_type = CompositeDistributionFunction<ippl::Vector<T, Dim>, T>;

    particle_container_type pcontainer_m;
    field_continer_type fcontainer_m;
    int np_m;
    vector_type center;

    InitialPlacement(particle_container_type& pc_, field_continer_type fc_) {
        this->pcontainer_m = pc_;
        this->fcontainer_m = fc_;
        this->np_m         = 0;

        this->center = 0.5 * (this->fcontainer_m->getRMax() - this->fcontainer_m->getRMin());
    }

    ~InitialPlacement() {}

    void initalizeParticles(vector_type nr_) {
        // Section 1: Place particles
        ParticleDistributionBase<T, Dim> filteredDist(this->fcontainer_m->getRMin(),
                                                      this->fcontainer_m->getRMax(),
                                                      new GridPlacement<T, Dim>(nr_));

        filteredDist.generateDistribution();

        this->np_m = filteredDist.getNumParticles();
        std::cout << "Number of particles: " << this->np_m << std::endl;

        // Section 2: Put them in the particle cointainer
        view_type particle_view = filteredDist.getParticles();

        this->pcontainer_m->create(this->np_m);
        Kokkos::parallel_for(
            "AddParticles", this->np_m, KOKKOS_LAMBDA(const int& i) {
                this->pcontainer_m->R(i)     = particle_view(i);
                this->pcontainer_m->omega(i) = 1;
            });

        Kokkos::fence();
    }
};

class EquidistantBand : public InitialPlacement {
public:
    EquidistantBand(particle_container_type& pc_, field_continer_type fc_)
        : InitialPlacement(pc_, fc_) {}

    void initalizeParticles(vector_type nr_) {
        // Section 1: Create composite distribution
        Band<T, Dim> band(1);

        ShiftTransformation<T, Dim> shift_(-this->center);
        band.applyTransformation(shift_);

        // Section 2: Place particles
        FilteredDistribution<T, Dim> filteredDist(band, this->fcontainer_m->getRMin(),
                                                  this->fcontainer_m->getRMax(),
                                                  new GridPlacement<T, Dim>(nr_));

        this->np_m = filteredDist.getNumParticles();
        std::cout << "Number of particles: " << this->np_m << std::endl;

        // Section 3: Put them in the particle cointainer
        view_type particle_view = filteredDist.getParticles();

        this->pcontainer_m->create(this->np_m);
        Kokkos::parallel_for(
            "AddParticles", this->np_m, KOKKOS_LAMBDA(const int& i) {
                this->pcontainer_m->R(i)     = particle_view(i);
                this->pcontainer_m->omega(i) = band.evaluate(this->pcontainer_m->R(i));
            });

        Kokkos::fence();
    }
};

class EquidistantCircle : public InitialPlacement {
public:
    EquidistantCircle(particle_container_type& pc_, field_continer_type fc_)
        : InitialPlacement(pc_, fc_) {}

    void initalizeParticles(vector_type nr_) {
        // Section 1: Create composite distribution
        Circle<T, Dim> circle(1);

        ShiftTransformation<T, Dim> shift_(-this->center);
        circle.applyTransformation(shift_);

        // Section 2: Place particles
        FilteredDistribution<T, Dim> filteredDist(circle, this->fcontainer_m->getRMin(),
                                                  this->fcontainer_m->getRMax(),
                                                  new GridPlacement<T, Dim>(nr_));

        this->np_m = filteredDist.getNumParticles();
        std::cout << "Number of particles: " << this->np_m << std::endl;

        // Section 3: Put them in the particle cointainer
        view_type particle_view = filteredDist.getParticles();

        this->pcontainer_m->create(this->np_m);
        Kokkos::parallel_for(
            "AddParticles", this->np_m, KOKKOS_LAMBDA(const int& i) {
                this->pcontainer_m->R(i)     = particle_view(i);
                this->pcontainer_m->omega(i) = circle.evaluate(this->pcontainer_m->R(i));
            });

        Kokkos::fence();
    }
};

class TwoCircles : public InitialPlacement {
public:
    TwoCircles(particle_container_type& pc_, field_continer_type fc_)
        : InitialPlacement(pc_, fc_) {}

    void initalizeParticles(vector_type nr_) {
        // Section 1: Create composite distribution
        Circle<T, Dim> circle1(1);
        Circle<T, Dim> circle2(1);

        ShiftTransformation<T, Dim> shift_(-this->center);
        circle1.applyTransformation(shift_);
        circle2.applyTransformation(shift_);

        vector_type y_disp(0);
        y_disp[1] = 1.5;
        ShiftTransformation<T, Dim> y_shift_p(y_disp);
        ShiftTransformation<T, Dim> y_shift_m(-y_disp);
        circle1.applyTransformation(y_shift_p);
        circle2.applyTransformation(y_shift_m);


        circle1 += circle2;


        // Section 2: Place particles
        FilteredDistribution<T, Dim> filteredDist(circle1, this->fcontainer_m->getRMin(),
                                                  this->fcontainer_m->getRMax(),
                                                  new GridPlacement<T, Dim>(nr_));

        this->np_m = filteredDist.getNumParticles();
        std::cout << "Number of particles: " << this->np_m << std::endl;

        // Section 3: Put them in the particle cointainer
        view_type particle_view = filteredDist.getParticles();

        this->pcontainer_m->create(this->np_m);
        Kokkos::parallel_for(
            "AddParticles", this->np_m, KOKKOS_LAMBDA(const int& i) {
                this->pcontainer_m->R(i)     = particle_view(i);
                this->pcontainer_m->omega(i) = circle1.evaluate(this->pcontainer_m->R(i));
            });

        Kokkos::fence();
    }
};

#endif

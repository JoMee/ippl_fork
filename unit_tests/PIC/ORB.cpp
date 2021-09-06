//
// Unit tests ORB for class OrthogonalRecursiveBisection
//   Test volume and charge conservation in PIC operations.
//
// Copyright (c) 2021, Michael Ligotino, ETH, Zurich; 
// Paul Scherrer Institut, Villigen; Switzerland
// All rights reserved
//
// This file is part of IPPL.
//
// IPPL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with IPPL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Ippl.h"

#include <cmath>
#include "gtest/gtest.h"

#include <random>

class ORBTest : public ::testing::Test {

public:
    static constexpr size_t dim = 3;
    typedef ippl::Field<double, dim> field_type;
    typedef ippl::FieldLayout<dim> flayout_type;
    typedef ippl::UniformCartesian<double, dim> mesh_type;
    typedef ippl::ParticleSpatialLayout<double, dim> playout_type;
    typedef ippl::OrthogonalRecursiveBisection<double, dim, mesh_type> ORB;

    template<class PLayout>
    struct Bunch : public ippl::ParticleBase<PLayout>
    {
        Bunch(PLayout& playout)
        : ippl::ParticleBase<PLayout>(playout)
        {
            this->addAttribute(Q);
        }
        
        ~Bunch(){ }
        
        typedef ippl::ParticleAttrib<double> charge_container_type;
        charge_container_type Q;
        
        void update() {
            PLayout& layout = this->getLayout();
            layout.update(*this);
        }

        void updateLayout(flayout_type fl, mesh_type mesh) {
            PLayout& layout = this->getLayout();
            layout.updateLayout(fl, mesh);
        }
    };


    typedef Bunch<playout_type> bunch_type;


    ORBTest()
    : nParticles(std::pow(256,3))
    , nPoints(512)
    {
        setup();
    }

    void setup() {
        ippl::Index I(nPoints);
        ippl::NDIndex<dim> owned(I, I, I);

        ippl::e_dim_tag allParallel[dim];    // Specifies SERIAL, PARALLEL dims
        for (unsigned int d = 0; d < dim; d++)
            allParallel[d] = ippl::PARALLEL;

        layout_m = flayout_type (owned, allParallel);

        double dx = 1.0 / double(nPoints);
        ippl::Vector<double, dim> hx = {dx, dx, dx};
        ippl::Vector<double, dim> origin = {0, 0, 0};

        mesh_m = mesh_type(owned, hx, origin);

        field = std::make_unique<field_type>(mesh_m, layout_m);

        pl_m = playout_type(layout_m, mesh_m);
        
        bunch = std::make_unique<bunch_type>(pl_m);
        
        int nRanks = Ippl::Comm->size();
        if (nParticles % nRanks > 0) {
            if (Ippl::Comm->rank() == 0) {
                std::cerr << nParticles << " not a multiple of " << nRanks << std::endl;
            }
            exit(1);
        }

        size_t nloc = nParticles / nRanks;
        bunch->create(nloc);
        
        std::mt19937_64 eng;
        eng.seed(42);
        eng.discard( nloc * Ippl::Comm->rank());
        std::uniform_real_distribution<double> unif(hx[0]/2, 1-(hx[0]/2));

        typename bunch_type::particle_position_type::HostMirror R_host = bunch->R.getHostMirror();
        for(size_t i = 0; i < nloc; ++i) {
            for (int d = 0; d<3; d++) {
                R_host(i)[d] =  unif(eng);
            }
        }

        Kokkos::deep_copy(bunch->R.getView(), R_host);

        orb.initialize(layout_m, mesh_m);
    }


    void repartition() {
        orb.binaryRepartition(bunch->R, layout_m);
        field->updateLayout(layout_m);
        bunch->updateLayout(layout_m, mesh_m);
    }

    ippl::NDIndex<dim> getDomain() {
        return layout_m.getDomain();
    }

    std::unique_ptr<field_type> field;
    std::unique_ptr<bunch_type> bunch;
    size_t nParticles;
    size_t nPoints;

private:
    flayout_type layout_m;
    mesh_type mesh_m;
    playout_type pl_m;
    ORB orb;
};

TEST_F(ORBTest, Volume) {
 
    ippl::NDIndex<dim> dom = getDomain();

    repartition();

    bunch->update();

    ippl::NDIndex<dim> ndom = getDomain();   

    ASSERT_DOUBLE_EQ(dom[0].length() * dom[1].length() * dom[2].length(), 
                     ndom[0].length() * ndom[1].length() * ndom[2].length());

}


TEST_F(ORBTest, Charge) {

    *field = 0.0;

    double charge = 0.5;

    bunch->Q = charge;

    bunch->update();

    scatter(bunch->Q, *field, bunch->R);

    repartition();

    bunch->update();

    double totalcharge = field->sum();

    ASSERT_DOUBLE_EQ(nParticles * charge, totalcharge);

}

int main(int argc, char *argv[]) {
    Ippl ippl(argc,argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
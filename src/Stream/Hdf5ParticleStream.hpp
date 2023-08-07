namespace ippl {

    namespace hdf5 {

        template <class ParticleContainer>
        void ParticleStream<ParticleContainer>::operator<<(const ParticleContainer& obj) {

            std::cout << "This is a HDF5 stream:" << std::endl;
            size_t nAttrib = obj.getAttributeNum();
            for (size_t i = 0; i < nAttrib; ++i) {
                auto attr = obj.getAttribute(i);
                std::cout << "Name:      " << attr->name() << std::endl;
                std::cout << "Long name: " << attr->long_name() << std::endl;
                std::cout << "Units:     " << attr->unit() << std::endl;
//                 std::cout << "Data:      ";
//                 for (int i = 0; i < attr->size(); ++i) {
//                     std::cout << (*attr)[i] << " ";
//                 }
//                 std::cout << std::endl;
            }
            std::cout << std::endl;

        }

        template <class ParticleContainer>
        void ParticleStream<ParticleContainer>::operator>>(ParticleContainer& /*obj*/) {


        }

    }

}  // namespace ippl

#ifndef SURFACE_COMMUNICATOR_H
#define SURFACE_COMMUNICATOR_H

#include <deal.II/base/tensor.h>
#include <iostream>
#include <fstream>
#include <vector>

namespace KirasFM{
  using namespace dealii;

template <int dim>
class SurfaceCommunicator {
  public:
  
    SurfaceCommunicator();

    SurfaceCommunicator(unsigned int n_domains_);

    SurfaceCommunicator(unsigned int n_domains_, std::string, unsigned int n_face_q_points, unsigned int n_rows);

    SurfaceCommunicator(const SurfaceCommunicator<dim>& copy);

    std::vector<
      std::vector<Tensor<1, dim, std::complex<double>>>
    > value(unsigned int i, unsigned int j);

    void value (
      std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> in,
      unsigned int i, unsigned int j
    );

    void update(SurfaceCommunicator<dim> sc, unsigned int i);

    SurfaceCommunicator<dim>& operator=(const SurfaceCommunicator<dim>& copy);

    void to_file (
      std::string name,
      unsigned int i,
      unsigned int j
    );
    

  private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar& n_domains;
        ar& n_faces;
        ar& value_data;
        v = version;
    }


    unsigned int n_domains;
    unsigned int n_faces;
    unsigned int v;

    std::vector<std::vector<
      std::vector<Tensor<1, dim, std::complex<double>>>
    >> value_data;
};

} // namespace KirasFM
#endif

#include <surface_communicator.h>

namespace KirasFM{
  using namespace dealii;

template <int dim>
SurfaceCommunicator<dim>::SurfaceCommunicator() {

}

template <int dim>
SurfaceCommunicator<dim>::SurfaceCommunicator(unsigned int n_domains_) {
  n_domains = n_domains_;

  // initialize value
  std::vector<
    std::vector<Tensor<1, dim, std::complex<double>>>
  > tmp_value;

  for(unsigned int i = 0; i < n_domains * n_domains; i++){
    value_data.push_back(tmp_value);
  }

}

// Copy Constructor (Rule of 3: 2/3)
template <int dim>
SurfaceCommunicator<dim>::SurfaceCommunicator(const SurfaceCommunicator<dim>& copy) {
  n_domains  = copy.n_domains;
  n_faces    = copy.n_faces;

  value_data = copy.value_data;
}

template <int dim>
std::vector<
  std::vector<Tensor<1, dim, std::complex<double>>>
> SurfaceCommunicator<dim>::value(unsigned int i, unsigned int j) {
  AssertIndexRange(i, n_domains);
  AssertIndexRange(j, n_domains);
  return value_data[(i * n_domains) + j];
}

template <int dim>
void SurfaceCommunicator<dim>::value(
  std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> in,
  unsigned int i, unsigned int j
) {
  AssertIndexRange(i, n_domains);
  AssertIndexRange(j, n_domains);
  value_data[(i * n_domains) + j] = in;
}

template <int dim>
void SurfaceCommunicator<dim>::update(SurfaceCommunicator<dim> sc, unsigned int i) {
  AssertIndexRange(i, n_domains);
  for ( unsigned int j = 0; j < n_domains; j++ ) { 
    value_data[(i * n_domains) + j] = sc.value(i, j);
  }
}

// Assaignment operator (Rule of 3: 3/3)
template <int dim>
SurfaceCommunicator<dim>& SurfaceCommunicator<dim>::operator=(const SurfaceCommunicator<dim>& copy) {
  if (this == &copy) 
    return *this;

  n_domains  = copy.n_domains;
  n_faces    = copy.n_faces;

  value_data = copy.value_data;
  
  return *this;
}       

template <int dim>
void SurfaceCommunicator<dim>::to_file(std::string name, unsigned int i, unsigned int j) {
  AssertIndexRange(i, n_domains);
  AssertIndexRange(j, n_domains);

  unsigned int subfaces = value_data[(i * n_domains) + j].size();
  if ( subfaces != 0 ) {
    unsigned int n_face_q_points = value_data[(i * n_domains) + j][0].size();
    std::ofstream file;
    file.open(name);
    // loop over all subfaces
    for ( unsigned int k = 0; k < subfaces; k++ ) {
      // loop over all quadrature points
      for ( unsigned int q_point = 0; q_point < n_face_q_points; q_point++ ) {
        file << value_data[(i * n_domains) + j][k][q_point];
        if ( q_point != n_face_q_points - 1 )
          file << "; ";
      }
      if ( k != n_faces - 1 )
        file << "\n";
    }
    file.close();

  } // fi

}

template class SurfaceCommunicator<2>;
template class SurfaceCommunicator<3>;

} // namespace: KirasFM

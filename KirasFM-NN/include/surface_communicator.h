/*---------------------------------------------------------------------
 *
 * Copyright (C) 2024 Sebastian Kinnewig
 *
 * The code is licensed under the GNU Lesser General Public License as
 * published by the Free Software Foundation in version 2.1
 * The full text of the license can be found in the file LICENSE.md
 *
 * ---------------------------------------------------------------------
 * Contact:
 *   Sebastian Kinnewig
 *   Leibniz Universität Hannover (LUH)
 *   Institut für Angewandte Mathematik (IfAM)
 *
 * Questions?
 *   E-Mail: kinnewig@ifam.uni-hannover.de
 *
 * Date: February 2023
 *       Update: August 2024
 */

#ifndef SURFACE_COMMUNICATOR_H
#define SURFACE_COMMUNICATOR_H

#include <deal.II/base/tensor.h>

#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <fstream>
#include <iostream>

namespace KirasFM
{
  using namespace dealii;

  template <int dim>
  class SurfaceCommunicator
  {
  public:
    // Constructor
    SurfaceCommunicator();
    SurfaceCommunicator(unsigned int n_domains_,
                        std::string,
                        unsigned int n_face_q_points,
                        unsigned int n_rows);
    SurfaceCommunicator(unsigned int n_domains_);
    SurfaceCommunicator(const SurfaceCommunicator<dim> &copy);

    // return functions
    std::vector<std::vector<Tensor<1, dim, std::complex<double>>>>
    value(unsigned int i, unsigned int j);

    std::vector<std::vector<std::complex<double>>>
    curl(unsigned int i, unsigned int j);

    // update functions
    void
    value(std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> in,
          unsigned int                                                   i,
          unsigned int                                                   j);

    void
    curl(std::vector<std::vector<std::complex<double>>> in,
         unsigned int                                   i,
         unsigned int                                   j);

    void
    update(SurfaceCommunicator<dim> sc, unsigned int i);

    // copy assignment operator
    SurfaceCommunicator<dim> &
    operator=(const SurfaceCommunicator<dim> &copy);

    void
    to_file(std::string name, unsigned int i, unsigned int j);

  private:
    friend class boost::serialization::access;

    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version)
    {
      ar & n_domains;
      ar & n_faces;
      ar & value_data;
      ar & curl_data;
      v = version;
    }

    unsigned int n_domains;
    unsigned int n_faces;
    unsigned int v;

    std::vector<std::vector<std::vector<Tensor<1, dim, std::complex<double>>>>>
      value_data;

    std::vector<std::vector<std::vector<std::complex<double>>>> curl_data;
  };

} // namespace KirasFM
#endif

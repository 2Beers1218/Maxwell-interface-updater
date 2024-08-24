/* ---------------------------------------------------------------------
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
 * Date: November 2021
 *       Update: August 2024
 */

#ifndef SURFACE_COMMUNICATOR_H
#define SURFACE_COMMUNICATOR_H

#include <deal.II/base/tensor.h>

#include <fstream>
#include <iostream>
#include <vector>

namespace KirasFM
{
  using namespace dealii;

  template <int dim>
  class SurfaceCommunicator
  {
  public:
    SurfaceCommunicator();

    SurfaceCommunicator(unsigned int n_domains_);

    SurfaceCommunicator(unsigned int n_domains_,
                        std::string,
                        unsigned int n_face_q_points,
                        unsigned int n_rows);

    SurfaceCommunicator(const SurfaceCommunicator<dim> &copy);

    std::vector<std::vector<Tensor<1, dim, std::complex<double>>>>
    value(unsigned int i, unsigned int j);

    void
    value(std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> in,
          unsigned int                                                   i,
          unsigned int                                                   j);

    void
    update(SurfaceCommunicator<dim> sc, unsigned int i);

    SurfaceCommunicator<dim> &
    operator=(const SurfaceCommunicator<dim> &copy);

    void
    to_file(std::string name, unsigned int i, unsigned int j);


  private:
    unsigned int n_domains;
    unsigned int n_faces;
    unsigned int v;

    std::vector<std::vector<std::vector<Tensor<1, dim, std::complex<double>>>>>
      value_data;
  };

} // namespace KirasFM
#endif

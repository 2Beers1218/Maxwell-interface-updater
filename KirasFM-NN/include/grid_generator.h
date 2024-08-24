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


#ifndef MY_GRID_GENERATOR_H
#define MY_GRID_GENERATOR_H

#include <fstream>

// Distributed grid generator
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

// Grid generator
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

namespace DDM_Grid_Generator
{
  using namespace dealii;

  template <int dim>
  class DDMGridGenerator
  {
  public:
    // constructor
    DDMGridGenerator(const unsigned int domain_id,
                     const unsigned int N_domains,
                     const unsigned int refinements);

    void
    make_simpleblock(Triangulation<dim> &in);

  private:
    // internal functions
    void
    set_boundary_ids(Triangulation<dim> &in);

    // Number of domains & domain_id
    const unsigned int domain_id;
    const unsigned int N_domains;
    const unsigned int refinements;
  };

} // namespace DDM_Grid_Generator

#endif

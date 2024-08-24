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

#include <grid_generator.h>

namespace DDM_Grid_Generator
{
  using namespace dealii;

  const unsigned int robin     = 0;
  const unsigned int dirichlet = 1;

  template <int dim>
  void
  mark_block(const Triangulation<dim> &triangulation,
             const Point<dim>          bottom_left,
             const Point<dim>          upper_right,
             types::material_id        id)
  {
    Assert(dim == 2 || dim == 3, ExcInternalError());

    Point<dim> center;
    for (auto &cell : triangulation.active_cell_iterators())
      {
        center = cell->center();

        if (dim == 2)
          {
            if (center(0) > bottom_left(0) && center(0) < upper_right(0) &&
                center(1) > bottom_left(1) && center(1) < upper_right(1))
              {
                cell->set_material_id(id);
              }
          }

        else if (dim == 3)
          {
            if (center(0) > bottom_left(0) && center(0) < upper_right(0) &&
                center(1) > bottom_left(1) && center(1) < upper_right(1) &&
                center(2) > bottom_left(2) && center(2) < upper_right(2))
              {
                cell->set_material_id(id);
              }
          }
      }
  }


  template <int dim>
  DDMGridGenerator<dim>::DDMGridGenerator(const unsigned int domain_id,
                                          const unsigned int N_domains,
                                          const unsigned int refinements)
    : domain_id(domain_id)
    , N_domains(N_domains)
    , refinements(refinements)
  {}

  template <int dim>
  void
  DDMGridGenerator<dim>::set_boundary_ids(Triangulation<dim> &in)
  {
    if (domain_id == 0 && domain_id == N_domains - 1)
      {
        for (auto &cell : in.cell_iterators())
          {
            cell->set_material_id(0);
            for (unsigned int face = 0;
                 face < GeometryInfo<dim>::faces_per_cell;
                 face++)
              {
                if (cell->face(face)->at_boundary())
                  {
                    switch (cell->face(face)->boundary_id())
                      {
                        case 2:
                          cell->face(face)->set_boundary_id(dirichlet);
                          break;

                        default:
                          cell->face(face)->set_boundary_id(robin); // normaly 1
                          break;
                      } // switch: boundary_id
                  } // fi: at_boundary
              } // rof: face
          } // rof: cell
      }
    else if (domain_id == 0)
      {
        for (auto &cell : in.cell_iterators())
          {
            cell->set_material_id(0);
            for (unsigned int face = 0;
                 face < GeometryInfo<dim>::faces_per_cell;
                 face++)
              {
                if (cell->face(face)->at_boundary())
                  {
                    switch (cell->face(face)->boundary_id())
                      {
                        case 2:
                          cell->face(face)->set_boundary_id(dirichlet);
                          break;

                        case 3:
                          cell->face(face)->set_boundary_id(domain_id + 3);
                          break;

                        default:
                          cell->face(face)->set_boundary_id(robin); // normaly 1
                          break;
                      } // switch: boundary_id
                  } // fi: at_boundary
              } // rof: face
          } // rof: cell
      }
    else if (domain_id == N_domains - 1)
      {
        for (auto &cell : in.cell_iterators())
          {
            cell->set_material_id(0);
            for (unsigned int face = 0;
                 face < GeometryInfo<dim>::faces_per_cell;
                 face++)
              {
                if (cell->face(face)->at_boundary())
                  {
                    switch (cell->face(face)->boundary_id())
                      {
                        case 2:
                          cell->face(face)->set_boundary_id(domain_id + 1);
                          break;

                        default:
                          cell->face(face)->set_boundary_id(robin); // normaly 1
                          break;
                      } // switch: boundary_id
                  } // fi: at_boundary
              } // rof: face
          } // rof: cell
      }
    else
      {
        for (auto &cell : in.cell_iterators())
          {
            cell->set_material_id(0);
            for (unsigned int face = 0;
                 face < GeometryInfo<dim>::faces_per_cell;
                 face++)
              {
                if (cell->face(face)->at_boundary())
                  {
                    switch (cell->face(face)->boundary_id())
                      {
                        case 2:
                          cell->face(face)->set_boundary_id(domain_id + 1);
                          break;

                        case 3:
                          cell->face(face)->set_boundary_id(domain_id + 3);
                          break;

                        default:
                          cell->face(face)->set_boundary_id(robin); // normaly 1
                          break;

                      } // switch: boundary_id
                  } // fi: at_boundary
              } // rof: face
          } // rof: cell
      }
  }

  template <int dim>
  void
  DDMGridGenerator<dim>::make_simpleblock(Triangulation<dim> &in)
  {
    const double absolut_lower_y = 0.0;
    const double absolut_upper_y = 1.0;

    double step    = (absolut_upper_y - absolut_lower_y) / (1.0 * N_domains);
    double lower_y = step * domain_id;
    double upper_y = step * (domain_id + 1);

    // left lower corner of the rectangle
    const Point<dim> left_edge =
      (dim == 2) ? Point<dim>(0.0, lower_y) : Point<dim>(0.0, lower_y, 0.0);

    // right upper corner of the rectangle
    const Point<dim> right_edge =
      (dim == 2) ? Point<dim>(1.0, upper_y) : Point<dim>(1.0, upper_y, 1.0);

    std::vector<unsigned int> n_subdivisions;
    n_subdivisions.push_back(N_domains);
    n_subdivisions.push_back(1);
    if (dim == 3)
      {
        n_subdivisions.push_back(N_domains);
      }

    // create the rectangle
    GridGenerator::subdivided_hyper_rectangle(
      in, n_subdivisions, left_edge, right_edge, true);

    in.refine_global(refinements);

    set_boundary_ids(in);

    // mark the waveguide inside of the grid
    const Point<dim> box_lower_left  = dim == 2 ?
                                         Point<dim>(2.0 / 8.0, 0.0) :
                                         Point<dim>(3.0 / 8.0, 0.0, 3.0 / 8.0);
    const Point<dim> box_upper_right = dim == 2 ?
                                         Point<dim>(6.0 / 8.0, 1.0) :
                                         Point<dim>(5.0 / 8.0, 1.0, 5.0 / 8.0);

    mark_block<dim>(in, box_lower_left, box_upper_right, 1);
  }

  template class DDMGridGenerator<2>;
  // template class DDMGridGenerator<3>;

} // namespace DDM_Grid_Generator

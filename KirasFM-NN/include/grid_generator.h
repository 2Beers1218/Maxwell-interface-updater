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

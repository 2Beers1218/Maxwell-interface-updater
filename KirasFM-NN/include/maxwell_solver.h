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

#ifndef MAXWELL_SOLVER_H 
#define MAXWELL_SOLVER_H 

// === deal.II includes ==
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

// Grid and triangulation
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_refinement.h>

#include <deal.II/distributed/shared_tria.h>

// dof handler
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

// FE - Elements:
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nedelec_sz.h>

// Dirichlet boundary condition
#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools.h>

// Solution transfer
#include <deal.II/numerics/solution_transfer.h>

// Solver
#include <deal.II/lac/trilinos_solver.h>

// Trilinos
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

// === C++ includes ===
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

// === KirasFM includes ===
#include <boundary_values.h>
#include <cross_product.h>
#include <parameter_reader.h>
#include <post_processing.h>
#include <surface_communicator.h>
#include <refinement_communicator.h>

namespace KirasFM {
  using namespace dealii;

// === The Maxwell Solver class ===
/*
 * This is a very simple Maxwell Solver, which does not provide much
 * functionality. It is mainly used to test features and for tracking
 * bugs.
 *
 * We aim to solve the partial differential equation
 *      curl ( curl ( E ) - \omega^2 E = f(x)             on \Omega
 *      trace( E )                     = \trace (E_{inc}} on \Gamma_inc
 *
 *
 */
template <int dim>
class MaxwellProblem {
  public:
    // constructor
    MaxwellProblem(
      ParameterReader    &param, 
      ConditionalOStream pcout,
      TimerOutput        timer,
      const unsigned int domain_id,
      const unsigned int N_domains,
      MPI_Comm           local_mpi_comm
    );

    MaxwellProblem(
      const MaxwellProblem &copy
    );

    // execute:
    void initialize();
    void solve(); 
    void print_results() const;
    void solution_to_file(std::string name);

    // assemble interface
    void update_interface_rhs();
    void assemble_interface_rhs();

    // return:
    Triangulation<dim>& return_triangulation();
    SurfaceCommunicator<dim> return_g_out();
    SurfaceCommunicator<dim> return_g_in();

    RefinementCommunicator<dim> return_refine();

    void setup_system();

    // assemble system
    void assemble_system();

    // interpolate
    void interpolate();
    void interpolate_global();
    void recall_solution();

    // refinement functions:
    void prepare_mark_interface_for_refinement();
    void apply_mark_interface_for_refinement();

    void refine();

    // update:
    void update_g_out( SurfaceCommunicator<dim> g );
    void update_g_in( SurfaceCommunicator<dim> g );

    void update_refine( RefinementCommunicator<dim> r );

  private:
    // === internal functions ===
    void setup_system( DoFHandler<dim>& dof_handler_local );


    // error estimator
    void error_estimator( Vector<float> &error_indicators ) const;

    // output
    void output_results() const;


    // MPI communicator
    MPI_Comm                       mpi_communicator;

    // out put
    ConditionalOStream             pcout;

    //timer
    TimerOutput                    timer;

    // content
    Triangulation<dim>             triangulation;
    DoFHandler<dim>                dof_handler;
    FESystem<dim>                  fe;

    AffineConstraints<double>      constraints; 

    // linear system
    TrilinosWrappers::SparseMatrix system_matrix;
    TrilinosWrappers::MPI::Vector  solution, system_rhs, rhs_backup, solution_backup;

    // locally_owned_dofs:
    IndexSet                       locally_owned_dofs;
    IndexSet                       locally_relevant_dofs;

    // configuration
    ParameterReader                prm;

    const std::complex<double>     imag_i = std::complex<double>(0.0, 1.0); 

    // Interface
    SurfaceCommunicator<dim>       SurfaceOperator, g_out;
    RefinementCommunicator<dim>    RefinementOperator; /* Refinement over intefaces */

    const unsigned int             domain_id;
    const unsigned int             N_domains;

    bool                           first_rhs = true;
    bool                           solved    = false;
    bool                           rebuild   = false;


};

} // namespace: KirasFM
#endif

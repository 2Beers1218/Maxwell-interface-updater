/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 *  Written by Sebastian Kinnewig
 *  November 2021
 *
 */

#include <maxwell_solver.h>

namespace KirasFM {
  using namespace dealii;

  const unsigned int robin = 0;
  const unsigned int dirichlet = 1;


// === constructor ===
// standart constructor:
template <int dim>
MaxwellProblem<dim>::MaxwellProblem(
  ParameterReader    &param,
  ConditionalOStream pcout,
  TimerOutput        timer,
  const unsigned int domain_id,
  const unsigned int N_domains,
  MPI_Comm           local_mpi_comm
):
  mpi_communicator( local_mpi_comm ),

  pcout( pcout ),

  timer( timer ),

  triangulation(
    mpi_communicator,
    typename Triangulation<dim>::MeshSmoothing (
      Triangulation<dim>::smoothing_on_refinement 
      | Triangulation<dim>::smoothing_on_coarsening
    )
  ),

  dof_handler( triangulation ),

  fe( FE_NedelecSZ<dim>( 
      param.get_integer(
        "Mesh & geometry parameters", 
        "Polynomial degree" ) ), 2 ),

  prm(param),

  // Surface Communicator
  g_in(SurfaceCommunicator<dim>(N_domains)),

  domain_id(domain_id),
  N_domains(N_domains),
  first_rhs(true)
{}

// copy constructor
template <int dim>
MaxwellProblem<dim>::MaxwellProblem (
  const MaxwellProblem<dim> &copy
):
  mpi_communicator( copy.mpi_communicator ),

  pcout( copy.pcout ),

  timer( copy.timer ),

  triangulation(
    mpi_communicator,
    typename Triangulation<dim>::MeshSmoothing (
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening
    )
  ),

  dof_handler( triangulation ),

  fe( FE_NedelecSZ<dim>( 
      copy.prm.get_integer(
        "Mesh & geometry parameters", 
        "Polynomial degree" ) ), 2 ),

  prm( copy.prm ),

  // Surface Communicator       
  g_in(SurfaceCommunicator<dim>(copy.N_domains)),

  domain_id(copy.domain_id),

  N_domains(copy.N_domains),

  first_rhs(copy.first_rhs)
{}


template <int dim>
void MaxwellProblem<dim>::setup_system() {
  timer.enter_subsection ( "Setup dof system" );
  pcout << "setup system...";

  // initialize the dof handler
  dof_handler.distribute_dofs(fe);

  // get the locally owned dofs
  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // Initialize the vectors:
  solution.reinit(
    locally_owned_dofs,
    locally_relevant_dofs,
    mpi_communicator
  );

  system_rhs.reinit(
    locally_owned_dofs, 
    mpi_communicator
  );

  rhs_backup.reinit(
    locally_owned_dofs, 
    mpi_communicator
  );

  // Constraits
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);

  // FE_Nedelec boundary condition.
  VectorTools::project_boundary_values_curl_conforming_l2(
    dof_handler,   
    0 /* vector component*/ , 
    DirichletBoundaryValues<dim>(), 
    dirichlet /* boundary id*/,  
    constraints
  ); 

  // FE_Nedelec boundary condition.
  VectorTools::project_boundary_values_curl_conforming_l2(
    dof_handler,   
    dim /* vector component*/ , 
    DirichletBoundaryValues<dim>(), 
    dirichlet /* boundary id*/,  
    constraints
  ); 

  constraints.close();

  // create sparsity pattern:
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(
    dof_handler, 
    dsp,
    constraints, 
    false
  );

  // create the system matrix
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    dof_handler.locally_owned_dofs(),
    mpi_communicator,
    locally_relevant_dofs
  );

  system_matrix.reinit(
    locally_owned_dofs,
    locally_owned_dofs,
    dsp,
    mpi_communicator
  );

  pcout << " done!" << std::endl;
  timer.leave_subsection();
}

/*
 * Assemble the system matrix: M
 * We decompose the system into real and imaginary part.
 *
 *         |  A   B |   | Re(E) |   | Re(f) |
 * M * E = |        | . |       | = |       |
 *         | -B   A |   | Im(E) |   | Im(f) |
 *
 * where E = Re(E) + i * Im(E), where i is the imaginary unit.
 *
 * And A = \int_{cell} curl( \phi_i(x) ) * curl( \mu * phi_j(x) ) dx 
 *       - \omega^2 \int_{cell} \phi_i(x) * \phi_j(x) dx  
 *
 *     B = \int_{face} \Trace( \psi_i(s) ) * \Trace( \psi_j(s) ) ds
 *
 * which corresponds to the maxwell equations, with robin boundary conditions 
 * in the weak form
 */
template <int dim>
void MaxwellProblem<dim>::assemble_system() {
  timer.enter_subsection ("Assemble system matrix");
  pcout << "assemble system matrix...\t";

  system_matrix         = 0;
  system_rhs            = 0;

  const unsigned int curl_dim = (dim == 2) ? 1 : 3;

  // choose the quadrature formulas
  QGauss<dim>     quadrature_formula(fe.degree + 2);
  QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

  // get the number of quadrature points and dofs
  const unsigned int n_q_points      = quadrature_formula.size(),
                     n_face_q_points = face_quadrature_formula.size(),
                     dofs_per_cell   = fe.dofs_per_cell;

  // set update flags
  FEValues<dim> fe_values(
    fe, quadrature_formula, 
    update_values 
    | update_gradients 
    | update_quadrature_points 
    | update_JxW_values
  );
  FEFaceValues<dim> fe_face_values(
    fe, face_quadrature_formula, 
    update_values 
    | update_quadrature_points 
    | update_normal_vectors 
    | update_gradients
    | update_JxW_values
  );

  // Extractors to real and imaginary parts
  const FEValuesExtractors::Vector E_re(0);
  const FEValuesExtractors::Vector E_im(dim);
  std::vector<FEValuesExtractors::Vector> vec(2);
  vec[0] = E_re;
  vec[1] = E_im;

  // create the local left hand side and right hand side
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // global material constants:
  const double omega      = prm.get_wavenumber();

  const std::complex<double> alpha(1.0, 0);

  // loop over all cells
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (cell->is_locally_owned() == false)
      continue;

    // initialize values:
    cell_matrix = 0;
    cell_rhs = 0;
    fe_values.reinit(cell);

    // cell dependent material constants:
    // for more details see: base_dir/mathematica/Material_Parameter.pdf
      // magnetic permability (approximativly = 1)
      const double mu_term  = 1.0;
      // electric permability: since mu_term = 1 => n = sqrt(eps_term) <=> eps_term = n^2
      const double eps_term = std::pow( prm.get_refrective_index( cell->material_id() ), 2 );

    // === Domain ===
    for ( const unsigned int i : fe_values.dof_indices() ) {
      const unsigned int block_index_i = fe.system_to_block_index(i).first;

      // we only want to compute this once
      std::vector<Tensor<1, dim>>      phi_i(n_q_points);
      std::vector<Tensor<1, curl_dim>> curl_phi_i(n_q_points);
      for (unsigned int q_point = 0; q_point < n_q_points; q_point++) {
        phi_i[q_point]        = fe_values[vec[block_index_i]].value(i, q_point);
        curl_phi_i[q_point]   = fe_values[vec[block_index_i]].curl(i, q_point);
      }
   
      for ( unsigned int j = i; j < dofs_per_cell; j++ ) { // we are using the symmetry of the here
        const unsigned int block_index_j = fe.system_to_block_index(j).first;

        /*
         * Assamble the block A, from the system matrix M
         *
         *  | A   0 |
         *  |       |
         *  | 0   A |
         *
         * where A = \int_{cell} curl( \phi_i(x) ) * curl( \mu * phi_j(x) ) dx 
         *         - \omega^2 \int_{cell} \phi_i(x) * \phi_j(x) dx  
         */
        if ( block_index_i != block_index_j ) 
          continue;

        double mass_part = 0;
        double curl_part = 0;

        for ( unsigned int q_point = 0; q_point < n_q_points; q_point++ ) {
          Tensor<1, dim> phi_j           = fe_values[vec[block_index_i]].value(j, q_point);
          Tensor<1, curl_dim> curl_phi_j = fe_values[vec[block_index_i]].curl(j, q_point);

          curl_part += 
            curl_phi_i[q_point]
            * curl_phi_j
            * fe_values.JxW(q_point);
          
          mass_part +=
            phi_i[q_point]
            * phi_j
            * fe_values.JxW(q_point);

        } // rof: q_point

        // Use skew-symmetry to fill matrices:
        cell_matrix(i,j) = ( mu_term * curl_part ) - ( ( omega * omega ) * eps_term * mass_part );
        cell_matrix(j,i) = ( mu_term * curl_part ) - ( ( omega * omega ) * eps_term * mass_part );

      } // rof: dof_j

    } // rof: dof_i

    // === boundary condition ===
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {
      fe_face_values.reinit(cell, face);

      if ( cell->face(face)->at_boundary() == false ) 
        continue;

      // --- Robin boundary conditions ---
      if ( cell->face(face)->boundary_id() == robin ) {
        /*
         * Assamble the block B
         *
         *  |  0   B |
         *  |        |
         *  | -B   0 |
         *
         *  where B = \int_{face} \Trace( \psi_i(s) ) * \Trace( \psi_j(s) ) ds
         */

        // compute the normal
        std::vector<Tensor<1, dim>> normal( n_face_q_points );
        for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++)
          normal[q_point] = fe_face_values.normal_vector(q_point);

        for ( const unsigned int i : fe_face_values.dof_indices() ) {
          const unsigned int block_index_i = fe.system_to_block_index(i).first;

          if ( fe.has_support_on_face(i, face) == false )
            continue;

          // we only want to compute this once
          std::vector<Tensor<1, dim>> phi_i(n_face_q_points);
          for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
            phi_i[q_point] = 
              CrossProduct::trace_tangential( fe_face_values[vec[block_index_i]].value(i, q_point), normal[q_point] );
          }

          for ( const unsigned int j : fe_face_values.dof_indices() ) {
            const unsigned int block_index_j = fe.system_to_block_index(j).first;

            if ( fe.has_support_on_face(j, face) == false )
              continue;

            double robin = 0;

            if ( block_index_i == block_index_j ) 
              continue;

            for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
              Tensor<1, dim> phi_j = CrossProduct::trace_tangential(fe_face_values[vec[block_index_j]].value(j, q_point), normal[q_point]);

              if ( block_index_i == 0 ) {
                robin -= 
                phi_i[q_point]  
                  * phi_j 
                  * fe_face_values.JxW(q_point);
              } else { // if ( block_index_i == 1 ) 
                robin += 
                  phi_i[q_point]  
                  * phi_j 
                  * fe_face_values.JxW(q_point);
              }

            } // rof: q_point

            cell_matrix(i, j) += omega * std::sqrt(eps_term) * robin;

          } // rof: dof_j

        } // rof: dof_i

      } // fi: Robin boundary condition

      // --- interface condition ---
      if ( cell->face(face)->boundary_id() >= 2) {

        // compute the normal
        std::vector<Tensor<1, dim>> normal( n_face_q_points );
        for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
          normal[q_point] = fe_face_values.normal_vector(q_point);
        }

        for ( const unsigned int i : fe_face_values.dof_indices() ) {
          const unsigned int block_index_i = fe.system_to_block_index(i).first;

          // we only want to compute this once
          std::vector<Tensor<1, dim>> phi_i( n_face_q_points );
          for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
            phi_i[q_point] = 
              CrossProduct::trace_tangential(
                fe_face_values[vec[block_index_i]].value(i, q_point), normal[q_point]);

           } // rof: q_point

          for ( const unsigned int j : fe_face_values.dof_indices() ) {
            const unsigned int block_index_j = fe.system_to_block_index(j).first;

            if ( fe.has_support_on_face(i, face) == false 
                 || fe.has_support_on_face(j, face) == false)
              continue;

            double robin           = 0;

            for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
              Tensor<1, dim> phi_j = 
                CrossProduct::trace_tangential(fe_face_values[vec[block_index_j]].value(j, q_point), normal[q_point]);

              // IBC ( = Robin):
              robin += 
                phi_i[q_point]  
                * phi_j 
                * fe_face_values.JxW(q_point);
              
            } // rof: q_point

            if (block_index_i == block_index_j) { // real values: (blocks: (0,0) and (1,1) )
            } else if ( block_index_i == 0 && block_index_j == 1) {
              cell_matrix(i, j) -= omega * robin;
            } else if (block_index_i == 1 && block_index_j == 0) { // if ( block_index_i == 1 ) 
              cell_matrix(i, j) += omega * robin;
            }

          } // rof: dof_j

        } // rof: dof_i

      } // fi: interface_condition

    } // rof: faces
  
    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global(
      cell_matrix, 
      cell_rhs, 
      local_dof_indices, 
      system_matrix, 
      system_rhs
    );

  } // rof: active cell

  // synchronization between all processors
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  // save the initial rhs;
  rhs_backup = system_rhs;
  //rhs_backup.compress(VectorOperation::add);

  pcout << " done!" << std::endl;
  timer.leave_subsection();
}


/* 
 * assemble_rhs sets the right hand side of the equation
 *      curl ( curl ( E ) ) - w^2 E = f(x)
 * where f(x) is defined via the function "CurlRHS" 
 * above.
 */
template <int dim>
void MaxwellProblem<dim>::assemble_system_rhs() {
  timer.enter_subsection ("Assemble system right hand side");
  pcout << "assemble right hand side...";

  // choose the quadrature formulas
  QGauss<dim>     quadrature_formula(fe.degree + 2);
  QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

  // get the number of quadrature points and dofs
  const unsigned int n_q_points         = quadrature_formula.size();
  const unsigned int dofs_per_cell      = fe.dofs_per_cell;

  // set update flags
  FEValues<dim> fe_values(
    fe, quadrature_formula, 
    update_values 
    | update_gradients 
    | update_quadrature_points 
    | update_JxW_values
  );

  // Extractors to real and imaginary parts
  const FEValuesExtractors::Vector E_re(0);
  const FEValuesExtractors::Vector E_im(dim);
  std::vector<FEValuesExtractors::Vector> vec(2);
  vec[0] = E_re;
  vec[1] = E_im;

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // exact solution
  CurlRHS<dim> curl_rhs;

  Vector<double> tmp(2 * dim);
  std::vector<Vector<double>> curl_rhs_values(n_q_points);
  for( unsigned int i = 0; i < n_q_points; i++) {
    curl_rhs_values[i] = tmp;
  }

  Tensor<1, dim> curl_rhs_tensor;

  const std::complex<double> imag_i(0.0, 1.0);

  Vector<double> cell_rhs(dofs_per_cell);

  for ( const auto &cell : dof_handler.active_cell_iterators() ) {

    if (cell->is_locally_owned() == false) 
      continue;

    fe_values.reinit(cell);
    cell->get_dof_indices(local_dof_indices);

    // local right hand side
    Vector<double> cell_rhs(dofs_per_cell);

    for ( const unsigned int i : fe_values.dof_indices() ) {
      const unsigned int block_index_i = fe.system_to_block_index(i).first;
      const unsigned int pos = block_index_i * dim;

      curl_rhs.vector_value_list( 
        fe_values.get_quadrature_points(), 
        curl_rhs_values
      );

      for (unsigned int q_point = 0; q_point < n_q_points; q_point++){

        curl_rhs_tensor[0] = curl_rhs_values[q_point][0 + pos];
        curl_rhs_tensor[1] = curl_rhs_values[q_point][1 + pos];
        if(dim == 3)
          curl_rhs_tensor[2] = curl_rhs_values[q_point][2 + pos];

        cell_rhs[i] += fe_values[vec[block_index_i]].value(i, q_point)
                       * curl_rhs_tensor * fe_values.JxW(q_point);

      } // rof: q_point

    } // rof: dof_i

    for ( const unsigned int i : fe_values.dof_indices() ) {
      system_rhs(local_dof_indices[i]) += cell_rhs[i];
    } // rof: dof_i

  } // rof: cell

  system_rhs.compress(VectorOperation::add);

  pcout << " done!" << std::endl;
  timer.leave_subsection();
}


template<int dim>
void MaxwellProblem<dim>::update_interface_rhs() {
  timer.enter_subsection ("update_interface");
  pcout << "assemble the interface...\t";

  // choose the quadrature formulas
  QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

  const unsigned int curl_dim = (dim == 2) ? 1 : 3;

  // get the number of quadrature points and dofs
  const unsigned int n_face_q_points    = face_quadrature_formula.size(); 
  const unsigned int dofs_per_cell      = fe.dofs_per_cell;

  // set update flags
  FEFaceValues<dim> fe_face_values(
    fe, face_quadrature_formula, 
    update_values 
    | update_quadrature_points 
    | update_normal_vectors 
    | update_gradients
    | update_JxW_values
  );

  // Extractors to real and imaginary parts
  const FEValuesExtractors::Vector E_re(0);
  const FEValuesExtractors::Vector E_im(dim);
  std::vector<FEValuesExtractors::Vector> vec(2);
  vec[0] = E_re;
  vec[1] = E_im;

  // Surface Operators:
  // s_tmp:
  std::vector<Tensor<1, dim, std::complex<double>>> update_value(n_face_q_points);
  std::vector<
    std::vector< std::vector<Tensor<1, dim, std::complex<double>>> >
  > s_value;
  {
    std::vector< std::vector<Tensor<1, dim, std::complex<double>>> > tmp;
    for (unsigned int k = 0; k < N_domains; k++)
      s_value.push_back(tmp); 
  }

  // Material constants:
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  double omega     = prm.get_wavenumber();

  const std::complex<double> beta(1.0 / (3.0 * omega * omega), 0.0);

  for( const auto &cell : dof_handler.active_cell_iterators() ) {
    //if (cell->is_locally_owned() == false) 
    //  continue;

    // cell dependent matirial constant
    //const double mu_term = 1.0;

    cell->get_dof_indices( local_dof_indices );


    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {

      unsigned int face_id = cell->face(face)->boundary_id();
      if (cell->face(face)->at_boundary() == false || face_id < 2)
        continue;

      fe_face_values.reinit(cell, face);

      // compute the normal
      std::vector<Tensor<1, dim>> normal( n_face_q_points );
      for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
        normal[q_point] = fe_face_values.normal_vector(q_point);
      }

      // face values:
      std::vector<Tensor<1,dim>>  E_value_re(n_face_q_points);
      std::vector<Tensor<1,dim>>  E_value_im(n_face_q_points);
      fe_face_values[E_re].get_function_values(solution, E_value_re);
      fe_face_values[E_im].get_function_values(solution, E_value_im);

      // face gradients:
      std::vector<Tensor<2,dim>>  E_grad_re(n_face_q_points);
      std::vector<Tensor<2,dim>>  E_grad_im(n_face_q_points);
      fe_face_values[E_re].get_function_gradients(solution, E_grad_re);
      fe_face_values[E_im].get_function_gradients(solution, E_grad_im);

      // face curls:
      std::vector<Tensor<1, curl_dim>>  E_curl_re(n_face_q_points);
      std::vector<Tensor<1, curl_dim>>  E_curl_im(n_face_q_points);
      fe_face_values[E_re].get_function_curls(solution, E_curl_re);
      fe_face_values[E_im].get_function_curls(solution, E_curl_im);

      for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++){
        update_value[q_point] = imag_i * omega * ( CrossProduct::trace_tangential(E_value_re[q_point], normal[q_point])
                              + (imag_i * CrossProduct::trace_tangential(E_value_im[q_point], normal[q_point])) );

        update_value[q_point] += ( CrossProduct::trace(E_curl_re[q_point], normal[q_point])
                        + (imag_i * CrossProduct::trace(E_curl_im[q_point], normal[q_point])) );
      }

      s_value[face_id - 2].push_back(update_value);

    } // rof: face

  } // rof: cell

  for( unsigned int face_id = 0; face_id < N_domains; face_id++) {
    g_in.value(s_value[face_id], domain_id, face_id);
  }

  pcout << " done!" << std::endl;
  timer.leave_subsection();
}

template <int dim>
void MaxwellProblem<dim>::assemble_interface_rhs() {
  timer.enter_subsection ("assemble_interface_rhs");
  pcout << "assmable the interface right hand side...\t";

  system_rhs = 0;

  // choose the quadrature formulas
  QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

  // get the number of quadrature points and dofs
  const unsigned int n_face_q_points    = face_quadrature_formula.size();
  const unsigned int dofs_per_cell      = fe.dofs_per_cell;

  // set update flags
  FEFaceValues<dim> fe_face_values(
    fe, face_quadrature_formula, 
    update_values 
    | update_gradients
    | update_quadrature_points 
    | update_normal_vectors 
    | update_JxW_values
  );

  // Extractors to real and imaginary parts
  const FEValuesExtractors::Vector E_re(0);
  const FEValuesExtractors::Vector E_im(dim);
  std::vector<FEValuesExtractors::Vector> vec(2);
  vec[0] = E_re;
  vec[1] = E_im;

  // SurfaceCommunicators:
  // g_tmp: (update_value)
  std::vector<Tensor<1, dim, std::complex<double>>>  update_value(n_face_q_points);
  std::vector<
    std::vector< std::vector<Tensor<1, dim, std::complex<double>>> >
  > g_value;
  {
    std::vector< std::vector<Tensor<1, dim, std::complex<double>>> > tmp;
    for (unsigned int k = 0; k < N_domains; k++)
      g_value.push_back(tmp);
  }

  // g_tmp: (update_curl)
  std::vector<std::complex<double>> update_curl(n_face_q_points);
  std::vector<
    std::vector< std::vector<std::complex<double>> >
  > g_curl;
  {
    std::vector< std::vector<std::complex<double>> > tmp;
    for (unsigned int k = 0; k < N_domains; k++)
      g_curl.push_back(tmp);
  }

  Vector<double> cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<unsigned int> face_n(N_domains, 0);

  for(const auto &cell : dof_handler.active_cell_iterators()) {

    // assamble the rhs of the interface condition
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {
      unsigned int face_id = cell->face(face)->boundary_id();
      if (cell->face(face)->at_boundary() == false || face_id < 2)
        continue;

      fe_face_values.reinit(cell, face);
      cell_rhs = 0;

      // compute the normal
      std::vector<Tensor<1, dim>> normal( n_face_q_points );
      for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
        normal[q_point] = fe_face_values.normal_vector(q_point);
      }

      for ( const unsigned int i : fe_face_values.dof_indices() ) {
      
        if ( fe.has_support_on_face(i, face) == false )
          continue;

        const unsigned int block_index_i = fe.system_to_block_index(i).first;

        for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
          // TODO: an dieser Stelle g_in.value erstezen durch das Ergebniss vom Neuronalen Netz.
          update_value[q_point] = g_in.value(face_id - 2, domain_id)[face_n[face_id - 2]][q_point];

          Tensor<1, dim> phi_tangential_i = 
            CrossProduct::trace_tangential(fe_face_values[vec[block_index_i]].value(i, q_point), normal[q_point]);

          std::complex<double> s_tmp = 
            (update_value[q_point] * phi_tangential_i) * fe_face_values.JxW(q_point);

          if ( cell->is_locally_owned() )  {
            if ( block_index_i == 0) { 
              cell_rhs[i] += s_tmp.real();
            } else {
              cell_rhs[i] += s_tmp.imag();
            }
          } // fi: locally owned

        } // rof: q_point

        cell->get_dof_indices(local_dof_indices);
        system_rhs(local_dof_indices[i]) -= cell_rhs[i];

      } // rof: dof_i

      face_n[face_id - 2]++;
    } // rof: face

  } // rof: cell

  system_rhs += rhs_backup;
  system_rhs.compress(VectorOperation::add);

  pcout << "done!" << std::endl;
  timer.leave_subsection();
}

template <int dim>
void MaxwellProblem<dim>::solve() {
  timer.enter_subsection ("Solve");
  pcout << "solve the linear system with MUMPS...";

  SolverControl solver_control( // Alternative:  ReductionControl 
    500,
    1e-6 * system_rhs.l2_norm()
  );

  TrilinosWrappers::SolverDirect::AdditionalData additional_data(
    false,         /* output_solver_details */
    // "Amesos_Mumps" /* solver type */
    "Amesos_Klu" /* solver type <-- if MUMPS is not available */
  );

  TrilinosWrappers::SolverDirect direct_solver(
    solver_control, 
    additional_data
  );

  TrilinosWrappers::MPI::Vector distributed_solution(
    locally_owned_dofs,
    mpi_communicator
  );

  constraints.set_zero(distributed_solution);

  direct_solver.solve(system_matrix, distributed_solution, system_rhs);
  
  constraints.distribute(distributed_solution);
  
  solution = distributed_solution;

  pcout << " done!" << std::endl;
  timer.leave_subsection();
}


template <int dim>
void MaxwellProblem<dim>::output_results() const {
  pcout << "write the results...";

  //Define objects of our ComputeIntensity class
  ComputeIntensity<dim> intensities;
  DataOut<dim>          data_out;
  
  // and a DataOut object:
  data_out.attach_dof_handler(dof_handler);

  const std::string filename = prm.get_string("Output parameters", "Output file");
  const std::string format   = ".vtu";
  const std::string outfile  =  filename + "-" + std::to_string(domain_id) + format;

  std::vector<std::string> solution_names;
  solution_names.emplace_back("Re_E1");
  solution_names.emplace_back("Re_E2");
  if( dim == 3 ) {
    solution_names.emplace_back("Re_E3");
  }
  solution_names.emplace_back("Im_E1");
  solution_names.emplace_back("Im_E2");
  if( dim == 3 ) {
    solution_names.emplace_back("Im_E3");
  }

  // mark Re_E1, Re_E2 (, Re_E3) as one vector
  std::vector<DataComponentInterpretation::DataComponentInterpretation> 
    data_component_interpretation (
      dim, DataComponentInterpretation::component_is_part_of_vector
    );

  // mark Im_E1, Im_E2 (, Im_E3) as one vector
  for ( unsigned int i = 0; i < dim; i++ ) 
    data_component_interpretation.push_back (
      DataComponentInterpretation::component_is_part_of_vector
    );

  data_out.add_data_vector(
    solution, 
    solution_names,
    DataOut<dim>::type_dof_data,
    data_component_interpretation
  );
  data_out.add_data_vector(solution, intensities);

  data_out.build_patches();
  data_out.write_vtu_in_parallel(outfile.c_str(), mpi_communicator);

  pcout << " done!" << std::endl;
}


// === public functions === 
template <int dim>
void MaxwellProblem<dim>::solution_to_file(std::string name) {

    std::ofstream file;
    file.open(name);

    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    FEFaceValues<dim> fe_face_values(
      fe, face_quadrature_formula,
      update_values
      | update_quadrature_points
      | update_normal_vectors
    );

    const FEValuesExtractors::Vector E_re(0);
    const FEValuesExtractors::Vector E_im(dim);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators()) {

        cell->get_dof_indices(local_dof_indices);

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; face++) {

          unsigned int face_id = cell->face(face)->boundary_id();
          if (cell->face(face)->at_boundary() && face_id >= 2) {

              fe_face_values.reinit(cell, face);

              std::vector<Tensor<1, dim>> normal( n_face_q_points );
              for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {
                normal[q_point] = fe_face_values.normal_vector(q_point);
              }

              std::vector<Tensor<1,dim>>  E_value_re(n_face_q_points);
              std::vector<Tensor<1,dim>>  E_value_im(n_face_q_points);
              fe_face_values[E_re].get_function_values(solution, E_value_re);
              fe_face_values[E_im].get_function_values(solution, E_value_im);

              for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++) {

                  file << ( CrossProduct::trace_tangential(E_value_re[q_point], normal[q_point])
                            + (imag_i * CrossProduct::trace_tangential(E_value_im[q_point], normal[q_point])) );

                  if (q_point != n_face_q_points - 1) {
                      file << "; ";
                  } else {
                      file << "\n";
                  }
              }
          }
        }
    }
    file.close();
}


template <int dim>
void MaxwellProblem<dim>::initialize() {
  setup_system();
  assemble_system();
  //assemble_rhs();
}

template <int dim>
void MaxwellProblem<dim>::print_results() const {
  output_results();
  timer.print_summary();
}

// return functions:
template <int dim>
parallel::shared::Triangulation<dim>& MaxwellProblem<dim>::return_triangulation() {
  return triangulation;
}

template <int dim>
SurfaceCommunicator<dim> MaxwellProblem<dim>::return_g_in() {
  return g_in;
}

// update
template<int dim>
void MaxwellProblem<dim>::update_g_in( SurfaceCommunicator<dim> g) {
  g_in = g;
}

// compile the tamplate with certain parameters
template class MaxwellProblem<2>;
//template class MaxwellProblem<3>;

} // namespace: KirasFM

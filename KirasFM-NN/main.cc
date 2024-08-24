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

// === C++ Includes ===
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <algorithm>
#include <iostream>

// === my includes ===
#include <grid_generator.h>
#include <maxwell_solver.h>
#include <refinement_communicator.h>
#include <surface_communicator.h>

// Grid
#include <deal.II/grid/grid_out.h>

namespace KirasFM
{
  using namespace dealii;

  bool
  doesFileExist(const std::string &filename)
  {
    std::ifstream file(filename.c_str());
    return file.good();
  }

  template <int dim>
  void
  make_grid(Triangulation<dim> &triangulation,
            unsigned            domain,
            unsigned            n_domains,
            unsigned int        refinements)
  {
    double start;
    double end;

    unsigned int              actual_refinements;
    std::vector<unsigned int> repetitions;
    if (n_domains == 2)
      {
        start = 1.0 * ((1.0 * domain));
        end   = 1.0 * ((1.0 * domain) + 1.0);

        repetitions        = (dim == 2) ? std::vector<unsigned int>{1, 1} :
                                          std::vector<unsigned int>{1, 1, 1};
        actual_refinements = refinements;
      }
    else if (n_domains == 4)
      {
        start = 0.5 * ((1.0 * domain));
        end   = 0.5 * ((1.0 * domain) + 1.0);

        repetitions        = (dim == 2) ? std::vector<unsigned int>{2, 1} :
                                          std::vector<unsigned int>{2, 1, 2};
        actual_refinements = refinements - 1;
      }
    else if (n_domains == 8)
      {
        start = 0.25 * ((1.0 * domain));
        end   = 0.25 * ((1.0 * domain) + 1.0);

        repetitions        = (dim == 2) ? std::vector<unsigned int>{4, 1} :
                                          std::vector<unsigned int>{4, 1, 4};
        actual_refinements = refinements - 2;
      }
    else if (n_domains == 16)
      {
        start = 0.125 * ((1.0 * domain));
        end   = 0.125 * ((1.0 * domain) + 1.0);

        repetitions        = (dim == 2) ? std::vector<unsigned int>{4, 1} :
                                          std::vector<unsigned int>{4, 1, 4};
        actual_refinements = refinements - 3;
      }
    else
      {
        Assert(false, ExcInternalError());
      }


    const Point<dim> left_edge =
      (dim == 2) ? Point<dim>(0.0, start) : Point<dim>(0.0, start, 0.0);

    const Point<dim> right_edge =
      (dim == 2) ? Point<dim>(1.0, end) : Point<dim>(1.0, end, 1.0);

    std::cout << "On domain " << domain << ": " << left_edge << ", "
              << right_edge << std::endl;

    // create the rectangle
    GridGenerator::subdivided_hyper_rectangle(
      triangulation, repetitions, left_edge, right_edge, true);

    // refine the grid
    triangulation.refine_global(actual_refinements);


    for (auto &cell : triangulation.active_cell_iterators())
      {
        cell->set_material_id(0);

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(0);
      }

    for (auto &cell : triangulation.active_cell_iterators())
      {
        Point<3>                  center(0.5, 0.5, 0.5);
        double                    distance_from_center = 0.0;
        std::vector<unsigned int> axis                 = {0, 2};
        for (unsigned int i = 0; i < dim - 1; ++i)
          distance_from_center +=
            std::pow(cell->center()[axis[i]] - center[axis[i]], 2.0);
        distance_from_center = std::sqrt(distance_from_center);

        double radius = 0.4;
        if (distance_from_center < radius)
          cell->set_material_id(1);

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            if (!cell->face(face)->at_boundary())
              continue;

            if (domain == 0)
              {
                if (cell->face(face)->center()[1] < start + 1e-8)
                  cell->face(face)->set_boundary_id(1);
                else if (cell->face(face)->center()[1] > end - 1e-8)
                  cell->face(face)->set_boundary_id(domain + 2 + 1);
                else
                  cell->face(face)->set_boundary_id(0);
              }
            else if (domain == (n_domains - 1))
              {
                if (cell->face(face)->center()[1] < start + 1e-8)
                  cell->face(face)->set_boundary_id(domain + 2 - 1);
                else if (cell->face(face)->center()[1] > end - 1e-8)
                  cell->face(face)->set_boundary_id(0);
                else
                  cell->face(face)->set_boundary_id(0);
              }
            else
              {
                if (cell->face(face)->center()[1] < start + 1e-8)
                  cell->face(face)->set_boundary_id(domain + 2 - 1);
                else if (cell->face(face)->center()[1] > end - 1e-8)
                  cell->face(face)->set_boundary_id(domain + 2 + 1);
                else
                  cell->face(face)->set_boundary_id(0);
              }
          }
      }

    // std::string name = "Grid-" + std::to_string(domain) + ".vtk";
    // std::ofstream output_file(name);
    // GridOut().write_vtk(triangulation, output_file);
  }


  template <int dim>
  class DDM
  {
  public:
    // constructor:
    DDM(boost::mpi::communicator world,
        ParameterReader         &prm,
        ConditionalOStream       pcout,
        TimerOutput              timer,
        const unsigned int       cpus_per_domain,
        const unsigned int       slizes);

    void
    initialize(std::vector<std::vector<unsigned int>> connectivity);
    void
    step(std::vector<std::vector<unsigned int>> connectivity,
         unsigned int                           step);

    // refinement methods:
    void
    prepare_refine(std::vector<std::vector<unsigned int>> connectivity);
    void
    refine();

    void
    print_result() const;

  private:
    // MPI
    boost::mpi::communicator world;
    const unsigned int       world_size;
    const unsigned int       world_rank;

    // Parameters
    ParameterReader prm;

    std::vector<unsigned int> proc_id_list;
    const unsigned int        cpus_per_domain;

    // information about the system
    const unsigned int slizes;
    const unsigned int size;

    // MPI communicatior
    MPI_Comm                               mpi_local_comm;
    std::vector<unsigned int>              owned_problems;
    std::vector<std::vector<unsigned int>> domain_map;

    // the system:
    std::vector<MaxwellProblem<dim>> thm;
    SurfaceCommunicator<dim>         g_out, g_in;

    RefinementCommunicator<dim> r_in;

    // help functions
    bool
    proc_first() const;
    unsigned int
    proc_previous() const;
    unsigned int
    proc_size() const;
  };



  template <int dim>
  DDM<dim>::DDM(boost::mpi::communicator world,
                ParameterReader         &prm,
                ConditionalOStream       pcout,
                TimerOutput              timer,
                const unsigned int       cpus_per_domain,
                const unsigned int       slizes)
    :

    world(world)
    , world_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , world_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , prm(prm)
    , cpus_per_domain(cpus_per_domain)
    , slizes(slizes)
    , size(slizes)
    , domain_map(std::vector<std::vector<unsigned int>>(size))
    , g_out(SurfaceCommunicator<dim>(size))
    , g_in(SurfaceCommunicator<dim>(size))
    , r_in(RefinementCommunicator<dim>(size))

  {
    if (world_size < cpus_per_domain)
      {
        pcout << "WARNING: You tried to assaign to each domain "
              << cpus_per_domain << " threads, but there are only "
              << world_size << "threads available!" << std::endl;
      }

    if (world_size > size * cpus_per_domain)
      {
        pcout << "WARNING: There are " << size
              << " tasks to solve, but there are " << world_size
              << " threads available, therefore " << world_size - size
              << " remain unused!" << std::endl;
      }

    // distribute the indicies fot the CPUs
    for (unsigned int proc_id = 0; proc_id * cpus_per_domain < world_size;
         proc_id++)
      {
        for (unsigned int i = 0; i < cpus_per_domain; i++)
          {
            proc_id_list.push_back(proc_id);
          }
      }

    // assign each CPU to a certain domain
    MPI_Comm_split(MPI_COMM_WORLD,
                   proc_id_list[world_rank],
                   world_rank,
                   &mpi_local_comm);

    // number of domains that can be evaluated in parallel:
    unsigned int n_chunks = world_size / cpus_per_domain;
    if (world_size % cpus_per_domain != 0)
      n_chunks++;

    // number of domains each chunk has to solve
    for (unsigned int i = 0; i < size; i++)
      {
        if (i % n_chunks == proc_id_list[world_rank])
          {
            unsigned int pos =
              (i / n_chunks) * n_chunks; // since i and n_chunks are integers,
                                         // (i / n_chunks) is adjusted downards.
            owned_problems.push_back(pos + proc_id_list[world_rank]);
          }
      }

    // create the domain map (a vector that contains on position domain_id the
    // process_id that owns this domain)
    std::vector<std::vector<unsigned int>> gathered_problems =
      Utilities::MPI::all_gather(MPI_COMM_WORLD, owned_problems);
    for (unsigned int i = 0; i < world_size; i++)
      {
        for (unsigned int j : gathered_problems[i])
          domain_map[j].push_back(i);
      }

    // create list of Maxwell Problems:
    for (unsigned int i : owned_problems)
      thm.push_back(
        MaxwellProblem<dim>(prm, pcout, timer, i, size, mpi_local_comm));
  }

  template <int dim>
  void
  DDM<dim>::initialize(std::vector<std::vector<unsigned int>> connectivity)
  {
    // set up grid:
    const unsigned int refinements =
      prm.get_integer("Mesh & geometry parameters", "Number of refinements");

    for (unsigned int i = 0; i < owned_problems.size(); i++)
      {
        DDM_Grid_Generator::DDMGridGenerator<dim> ddm_gg(owned_problems[i], size, refinements);

        // Simple cube:
        ddm_gg.make_simpleblock(thm[i].return_triangulation());

	      //// Simple waveguide:
        //ddm_gg.make_waveguide(thm[i].return_triangulation(),
        //                      owned_problems[i],
        //                      size,
        //                      refinements);
      }

    // initalize the maxwell problems:
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].initialize();

    // we already solve it the first time here:
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].solve();
  }

  template <int dim>
  void
  DDM<dim>::step(std::vector<std::vector<unsigned int>> connectivity,
                 unsigned int                           step)
  {
    // update the interfaces (compute g_in)
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].update_interface_rhs();

    // Gather g_in:
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      g_in.update(thm[i].return_g_in(), owned_problems[i]);

    // distribute g_in to all neighbours:
    if (proc_first())
      {
        for (unsigned int problem : owned_problems)
          {
            for (unsigned int dest : connectivity[problem])
              {
                const unsigned int       i = domain_map[dest][0];
                SurfaceCommunicator<dim> g_tmp(size);

                // send & recieve
                boost::mpi::request reqs[2];
                reqs[0] = world.isend(i, dest, g_in);
                reqs[1] = world.irecv(i, problem, g_tmp);
                boost::mpi::wait_all(reqs, reqs + 2);
                g_in.update(g_tmp, dest);
              }
          }

        for (unsigned int i = 1; i < proc_size(); i++)
          world.send(world_rank + i, world_rank + i, g_in);
      }
    else
      {
        SurfaceCommunicator<dim> g_tmp(size);
        world.recv(proc_previous(), world_rank, g_tmp);
        g_in = g_tmp;
      }


    bool nn_files_exists = false;
    {
      std::string file_name_1 =
        "interface_" + std::to_string(1) + "_" + std::to_string(0) + "_nn.csv";
      std::string file_name_2 =
        "interface_" + std::to_string(0) + "_" + std::to_string(1) + "_nn.csv";
      if (doesFileExist(file_name_1) && doesFileExist(file_name_2))
        nn_files_exists = true;
    }

    if (nn_files_exists)
      {
        std::cout << "Use Interface computed by NN:" << std::endl;

        for (unsigned int i = 0; i < owned_problems.size(); i++)
          {
            std::string file_name = "interface_";

            // --- Vorbereitung ---
            // Zum wieder einlesen der Datei benoetigst du:
            //  - Wie viele Zahlen stehen in einer Reihe, das ist besimmt durch:
            //  n_face_q_points
            //  - Wie viele Zeilen hat die Matrix: ich nenne die Variable hier:
            //  n_rows (beim schreiben der Datei, nennst du die Variable
            //  subfaces)
            SurfaceCommunicator<dim> g_in_nn(size,
                                             file_name,
                                             8 /*n_face_q_points*/,
                                             16 /*n_rows*/);
            thm[i].update_g_in(g_in_nn);
          }

        std::cout << "done!" << std::endl;

        // Debugging:
        // for( unsigned int i = 0; i < size; i++ )
        //  for( unsigned int j = 0; j < size; j++ ) {
        //    std::string file_name2 = "interface_" + std::to_string(i) + "_" +
        //    std::to_string(j) + "_" + "control" + ".csv";
        //    g_in_nn.to_file(file_name2, i, j);
        //  }
      }
    else
      {
        // Print to file
        if (proc_first())
          {
            for (unsigned int i = 0; i < size; i++)
              for (unsigned int j = 0; j < size; j++)
                {
                  std::string file_name = "interface_" + std::to_string(i) +
                                          "_" + std::to_string(j) + "_" +
                                          std::to_string(step) + ".csv";
                  g_in.to_file(file_name, i, j);
                }
          }

        // Distribute g_in:
        for (unsigned int i = 0; i < owned_problems.size(); i++)
          thm[i].update_g_in(g_in);
      }


    // compute with the new g_in the interface right hand side (and compute
    // g_out)
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].assemble_interface_rhs();

    // Solve the system
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].solve();

    if (size == 2)
      {
        for (unsigned int i = 0; i < owned_problems.size(); i++)
          {
            std::string file_name = "solution_" +
                                    std::to_string(owned_problems[i]) + "_" +
                                    std::to_string(step) + ".csv";
            thm[i].solution_to_file(file_name);
          }
      }
    else
      {
        std::string file_name0 =
          "solution_" + std::to_string(0) + "_" + std::to_string(step) + ".csv";
        thm[0].solution_to_file(file_name0);
        std::string file_name1 = "solution_" + std::to_string(size - 1) + "_" +
                                 std::to_string(step) + ".csv";
        thm[size - 1].solution_to_file(file_name1);
        for (unsigned int i = 1; i < size - 1; i++)
          {
            std::string file_name_above = "solution_" + std::to_string(i) +
                                          "_" + std::to_string(step) +
                                          "_above.csv";
            thm[i].solution_to_file(file_name_above);
            std::string file_name_below = "solution_" + std::to_string(i) +
                                          "_" + std::to_string(step) +
                                          "_below.csv";
            thm[i].solution_to_file(file_name_below);
          }
      }
  }

  template <int dim>
  void
  DDM<dim>::prepare_refine(std::vector<std::vector<unsigned int>> connectivity)
  {
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].prepare_mark_interface_for_refinement();

    // Gather refinement:
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      r_in.update(thm[i].return_refine(), owned_problems[i]);

    // distribute r_in to all neighbours:
    if (proc_first())
      {
        for (unsigned int problem : owned_problems)
          {
            for (unsigned int dest : connectivity[problem])
              {
                const unsigned int          i = domain_map[dest][0];
                RefinementCommunicator<dim> r_tmp(size);

                // send & recieve
                boost::mpi::request reqs[2];
                reqs[0] = world.isend(i, dest, r_in);
                reqs[1] = world.irecv(i, problem, r_tmp);
                boost::mpi::wait_all(reqs, reqs + 2);
                r_in.update(r_tmp, dest);
              }
          }

        for (unsigned int i = 1; i < proc_size(); i++)
          world.send(world_rank + i, world_rank + i, r_in);
      }
    else
      {
        RefinementCommunicator<dim> r_tmp(size);
        world.recv(proc_previous(), world_rank, r_tmp);
        r_in = r_tmp;
      }

    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].update_refine(r_in);

    // mark the cells at the interface for refinment:
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].apply_mark_interface_for_refinement();
  }

  template <int dim>
  void
  DDM<dim>::refine()
  {
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].refine();
  }

  template <int dim>
  void
  DDM<dim>::print_result() const
  {
    for (unsigned int i = 0; i < owned_problems.size(); i++)
      thm[i].print_results();
  }

  template <int dim>
  bool
  DDM<dim>::proc_first() const
  {
    if (world_rank == 0)
      return true;

    if (proc_id_list[world_rank] == proc_id_list[world_rank - 1])
      return false;

    return true;
  }

  template <int dim>
  unsigned int
  DDM<dim>::proc_previous() const
  {
    unsigned int out = 0;
    for (unsigned int i : proc_id_list)
      {
        if (i == proc_id_list[world_rank])
          return out;
        out++;
      }

    return out;
  }

  template <int dim>
  unsigned int
  DDM<dim>::proc_size() const
  {
    unsigned int out = 0;
    for (unsigned int i : proc_id_list)
      if (i == proc_id_list[world_rank])
        out++;

    return out;
  }

} // namespace KirasFM



int
main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace KirasFM;
      using namespace DDM_Grid_Generator;

      // deal.II MPI interface
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      // boost MPI interface
      boost::mpi::environment  env(argc, argv);
      boost::mpi::communicator world;

      ParameterHandler param;
      ParameterReader  prm(param);
      prm.read_parameters("options.prm");

      const unsigned int dim =
        prm.get_integer("Mesh & geometry parameters", "Dimension");

      const unsigned int cpus_per_domain =
        prm.get_integer("MPI parameters", "CPUs per domain");

      const unsigned int slizes = prm.get_integer("MPI parameters", "slizes");

      // Initialize an output, that will be only printed once over all threads,
      // we are initializing this already here, for the late use with domain
      // decomposition
      ConditionalOStream pcout(
        std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

      TimerOutput timer(pcout, TimerOutput::never, TimerOutput::wall_times);

      // create the connectivity map
      const unsigned int                     size = slizes;
      std::vector<std::vector<unsigned int>> connectivity(size);

      for (unsigned int i = 0; i < size; i++)
        {
          if (i != 0)
            connectivity[i].push_back(i - 1);

          if (i != size - 1)
            connectivity[i].push_back(i + 1);
        }

      //    // --- For debugging ---
      //    pcout << "Connectivity map:" << std::endl;
      //    for(std::vector<unsigned int> bla : connectivity) {
      //      for(unsigned int i : bla)
      //        pcout << i << " ";
      //      pcout << std::endl;
      //    }

      //    // --- For debugging ---
      //    // print the grid:
      //    Triangulation<3> tria;

      //    KirasFM_Grid_Generator::KirasFMGridGenerator<3> ddm_gg(0, 3 * 8, 2);
      //    std::vector<double> layer_thickness = {2.00, 1.4, 0.6};
      //    ddm_gg.create_parted_nano_particle(
      //      tria,
      //      1.0, /*  radius of the silver ball */
      //      layer_thickness
      //    );

      //    std::ofstream output_file1("Grid1.vtk");
      //    GridOut().write_vtk(tria, output_file1);


      pcout
        << " __    __  __                               ________  __       __ \n"
        << "/  |  /  |/  |                             /        |/  \\     /  |\n"
        << "$$ | /$$/ $$/   ______   ______    _______ $$$$$$$$/ $$  \\   /$$ |\n"
        << "$$ |/$$/  /  | /      \\ /      \\  /       |$$ |__    $$$  \\ /$$$ |\n"
        << "$$  $$<   $$ |/$$$$$$  |$$$$$$  |/$$$$$$$/ $$    |   $$$$  /$$$$ |\n"
        << "$$$$$  \\  $$ |$$ |  $$/ /    $$ |$$      \\ $$$$$/    $$ $$ $$/$$ |\n"
        << "$$ |$$  \\ $$ |$$ |     /$$$$$$$ | $$$$$$  |$$ |      $$ |$$$/ $$ |\n"
        << "$$ | $$  |$$ |$$ |     $$    $$ |/     $$/ $$ |      $$ | $/  $$ |\n"
        << "$$/   $$/ $$/ $$/       $$$$$$$/ $$$$$$$/  $$/       $$/      $$/ \n"
        << std::endl;

      switch (dim)
        {
          case 2:
            {
              DDM<2> problem(world, prm, pcout, timer, cpus_per_domain, slizes);
              pcout
                << "=================================================================="
                << std::endl;
              pcout << "INITIALIZE:" << std::endl;
              problem.initialize(connectivity);
              problem.print_result();

              for (unsigned int i = 0;
                   i < prm.get_integer("Mesh & geometry parameters",
                                       "Number of global iterations");
                   i++)
                {
                  pcout
                    << "=================================================================="
                    << std::endl;
                  pcout << "STEP " << i + 1 << ":" << std::endl;
                  problem.step(connectivity, i);
                  problem.print_result();
                }
              pcout
                << "=================================================================="
                << std::endl;
            }
            break;

          default:
            Assert(false, ExcNotImplemented());
            break;
        }
    }

  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}

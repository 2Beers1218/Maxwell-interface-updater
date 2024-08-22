  // === C++ Includes ===
#include <iostream>

  // === my includes ===
#include <surface_communicator.h>
#include <grid_generator.h>
#include <maxwell_solver.h>

namespace KirasFM {
  using namespace dealii;

  
  template<int dim>
  class DDM {
    public:
      //constructor:
      DDM(
        ParameterReader &prm,
        MPI_Comm mpi_local,
        ConditionalOStream pcout,
        TimerOutput timer,
        const unsigned int slizes
      );

      void initialize();
      void step(int steps);
      void print_result() const;

    private:

      // Parameters
      ParameterReader prm;

      // MPI communicatior
      MPI_Comm mpi_local_comm;

      // information about the system
      const unsigned int size;

      // the system:
      std::vector<MaxwellProblem<dim>> thm;
      SurfaceCommunicator<dim> g_in;
  };

  template<int dim>
  DDM<dim>::DDM (
    ParameterReader &prm,
    MPI_Comm mpi_local,
    ConditionalOStream pcout,
    TimerOutput timer,
    const unsigned int slizes
  ) :

    prm(prm), 

    mpi_local_comm(mpi_local),

    size(slizes),

    g_in(SurfaceCommunicator<dim>(size))

  {
    for( unsigned int i = 0; i < size; i++)  
      thm.push_back( MaxwellProblem<dim>(prm, pcout, timer, i, size , mpi_local_comm) );
    
  }

  template<int dim>
  void DDM<dim>::initialize() {
    // set up grid:
    const unsigned int refinements = prm.get_integer("Mesh & geometry parameters", "Number of refinements");

    for( unsigned int i = 0; i < size; i++ ) {
      // Simple Block Benchmark (2D & 3D)
        DDM_Grid_Generator::DDMGridGenerator<dim> ddm_gg(i, size, refinements);
        ddm_gg.make_simpleblock( thm[i].return_triangulation() );
    }

    // initalize the maxwell problems:
    for( unsigned int i = 0; i < size; i++ ) {
      thm[i].initialize();
    }

    // we already solve it the first time here:
    for( unsigned int i = 0; i < size; i++ ) 
      thm[i].solve();

  }

  template<int dim>
  void DDM<dim>::step(int steps) {
    
    // update the interfaces (compute g_in)
    for( unsigned int i = 0; i < size; i++ ) 
      thm[i].update_interface_rhs();

    // Gather g_in:
    for( unsigned int i = 0; i < size; i++ ) 
      g_in.update(thm[i].return_g_in(), i);

    // Print to file
    for( unsigned int i = 0; i < size; i++ ) 
      for( unsigned int j = 0; j < size; j++ ) {
        std::string file_name = "interface_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(steps) + ".csv";
        g_in.to_file(file_name, i, j);
      }

    // -->At this point you should also read in the output from a neural network

    // Distribute g_in:
    for( unsigned int i = 0; i < size; i++ ) 
      thm[i].update_g_in(g_in);

    // compute with the new g_in the interface right hand side (and compute g_out)
    for( unsigned int i = 0; i < size; i++ ) 
      thm[i].assemble_interface_rhs();

    // Solve the system
    for( unsigned int i = 0; i < size; i++ ) 
      thm[i].solve();

    if ( steps % 4 == 0 ) {
      unsigned int tmp_step = step / 4;
      for( unsigned int i = 0; i < size; i++ ) {
          std::string file_name = "solution_" + std::to_string(i) + "_" + std::to_string(steps) + ".csv";
          thm[i].solution_to_file(file_name);
      }

  }

  template<int dim>
  void DDM<dim>::print_result() const {
    for( unsigned int i = 0; i < size; i++ ) 
      thm[i].print_results();

    // for debugging: print the grid:
    //for( unsigned int i = 0; i < owned_problems.size(); i++ ) {
    //  std::string out_name = "Grid-";
    //  out_name += std::to_string(owned_problems[i]);
    //  out_name += ".vtk";
    //  std::ofstream output_file1(out_name);
    //  GridOut().write_vtk(thm[i].return_triangulation(), output_file1);
    //}

  }

} // namespace: KirasFM

int main(int argc, char *argv[]) {

  try {
    using namespace dealii;
    using namespace KirasFM;
    using namespace DDM_Grid_Generator;

    // deal.II MPI interface
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    ParameterHandler param;
    ParameterReader  prm(param);
    prm.read_parameters("options.prm");

    const unsigned int dim = prm.get_integer (
        "Mesh & geometry parameters", 
        "Dimension" );

    const unsigned int slizes = prm.get_integer (
        "MPI parameters",
        "slizes" );


    // Initialize an output, that will be only printed once over all threads,
    // we are initializing this already here, for the late use with domain decomposition
    ConditionalOStream pcout(
      std::cout, 
      (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
    );

    TimerOutput timer (
      pcout, 
      TimerOutput::never,
      TimerOutput::wall_times
    );

    pcout << " __    __  __                               ________  __       __ \n" 
          << "/  |  /  |/  |                             /        |/  \\     /  |\n" 
          << "$$ | /$$/ $$/   ______   ______    _______ $$$$$$$$/ $$  \\   /$$ |\n" 
          << "$$ |/$$/  /  | /      \\ /      \\  /       |$$ |__    $$$  \\ /$$$ |\n" 
          << "$$  $$/   $$ |/$$$$$$  |$$$$$$  |/$$$$$$$/ $$    |   $$$$  /$$$$ |\n"
          << "$$$$$  \\  $$ |$$ |  $$/ /    $$ |$$      \\ $$$$$/    $$ $$ $$/$$ |\n"  
          << "$$ |$$  \\ $$ |$$ |     /$$$$$$$ | $$$$$$  |$$ |      $$ |$$$/ $$ |\n" 
          << "$$ | $$  |$$ |$$ |     $$    $$ |/     $$/ $$ |      $$ | $/  $$ |\n" 
          << "$$/   $$/ $$/ $$/       $$$$$$$/ $$$$$$$/  $$/       $$/      $$/ \n" 
          << std::endl;

    switch ( dim ) {
      case 2: {
          DDM<2> problem(prm, MPI_COMM_WORLD, pcout, timer, slizes);
          pcout << "==================================================================" << std::endl;
          pcout << "INITIALIZE:" << std::endl;
          problem.initialize();
          for( unsigned int i = 0; i < prm.get_integer("Mesh & geometry parameters", "Number of global iterations"); i++ ) {
            pcout << "==================================================================" << std::endl;
            pcout << "STEP " << i + 1 << ":" << std::endl;
            problem.step(i);
          }
          pcout << "==================================================================" << std::endl;
          problem.print_result();
        } break;

      default:
        Assert(false, ExcNotImplemented());
        break;
    }

  }

  catch (std::exception &exc) {
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
  catch (...) {
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

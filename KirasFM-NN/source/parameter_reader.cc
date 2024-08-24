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
#include <parameter_reader.h>

namespace KirasFM
{

  // === The Parameter Reader class (derived from ParameterHandler) ===
  ParameterReader::ParameterReader(ParameterHandler &param)
    : prm(param)
  {}

  // === declare the parameters ===
  void
  ParameterReader::declare_parameters()
  {
    // === mesh & geometry ===
    prm.enter_subsection("Mesh & geometry parameters");
    {
      prm.declare_entry("Dimension",
                        "2",
                        Patterns::Integer(2, 3),
                        "The dimension"
                        "of the used Nedelec elements");

      prm.declare_entry("Number of refinements",
                        "6",
                        Patterns::Integer(0),
                        "Number of globar mesh refinement steps"
                        "applied to inital coarse grid");

      prm.declare_entry("Polynomial degree",
                        "1",
                        Patterns::Integer(0, 10),
                        "Polynomial degree"
                        "of the used Nedelec elements");

      prm.declare_entry("Number of global iterations",
                        "1",
                        Patterns::Integer(0),
                        "Number of outer jacobi fixed point iterations"
                        "of the used Nedelec elements");

      prm.declare_entry("Size of grid",
                        "1.0",
                        Patterns::Double(0),
                        "Scaling factor of the grid");
    }
    prm.leave_subsection();


    prm.enter_subsection("Incident field");
    {
      prm.declare_entry(
        "beam width real",
        "100.0",
        Patterns::Double(0),
        "Width of the incident Gauss beam (used in the real part)");

      prm.declare_entry(
        "beam width imag",
        "100.0",
        Patterns::Double(0),
        "Width of the incident Gauss beam (used in the imaginary part)");


      prm.declare_entry("beam center real",
                        "0.5",
                        Patterns::Double(0),
                        "Center of the Gauss beam (used in the real part)");

      prm.declare_entry(
        "beam center imag",
        "0.5",
        Patterns::Double(0),
        "Center of the Gauss beam (used in the imaginary part)");

      prm.declare_entry("beam intensity real",
                        "1.0",
                        Patterns::Double(0),
                        "Intensity (of the real part)");

      prm.declare_entry("beam intensity imag",
                        "0.0",
                        Patterns::Double(0),
                        "Intensity (of the imaginary part)");
    }
    prm.leave_subsection();


    // === MPI parameter ===
    prm.enter_subsection("MPI parameters");
    {
      prm.declare_entry("CPUs per domain",
                        "1",
                        Patterns::Integer(0),
                        "Number of threads assigned to each domain");

      prm.declare_entry("slizes",
                        "1",
                        Patterns::Integer(0),
                        "Number of slizes (length) that will be created");
    }
    prm.leave_subsection();


    // === material constants ===
    prm.enter_subsection("Material parameters");
    {
      prm.declare_entry("List of refrective indicies",
                        "1.0",
                        Patterns::Anything(),
                        "A list of values seperated by whitespaces ");
    }
    prm.leave_subsection();

    // === physical constants ===
    prm.enter_subsection("Physical constants");
    {
      prm.declare_entry("lambda",
                        "1.0",
                        Patterns::Double(0),
                        "Frequency we want to observe");
    }
    prm.leave_subsection();

    // === filename & format ===
    prm.enter_subsection("Output parameters");
    {
      prm.declare_entry("Output file",
                        "solution",
                        Patterns::Anything(),
                        "Name of the output file (without extension)");

      // declare parameters handels all the parameters needed for a certain
      // return type
      DataOutInterface<1>::declare_parameters(prm);
    }
    prm.leave_subsection();
  }



  // ==== read the refrective index list ===
  void
  ParameterReader::interpret_refrectiv_index()
  {
    std::vector<std::string> sub_sec;
    sub_sec.push_back("Material parameters");
    std::string refrectiv_index =
      prm.get(sub_sec, "List of refrective indicies");

    std::istringstream string_stream(refrectiv_index);
    double             tmp_refrectiv_index;

    std::vector<double> tmp_refrectiv_index_list;
    while (string_stream >> tmp_refrectiv_index)
      tmp_refrectiv_index_list.push_back(tmp_refrectiv_index);

    for (unsigned int i = 0; i < tmp_refrectiv_index_list.size() / 2; i++)
      refrectiv_index_list.push_back(
        std::complex<double>(tmp_refrectiv_index_list[2 * i + 0],
                             tmp_refrectiv_index_list[2 * i + 1]));

    Assert(refrectiv_index_list.size() != 0,
           ExcMessage(
             "The list of refrectives indicies has the wrong format, the"
             "list must be a list of float numbers seperated by white spaces, "
             "for example: "
             "1.0 1.51 1.41"));
  }



  // === read parameters ===
  void
  ParameterReader::read_parameters(const std::string &parameter_file)
  {
    // read the parameters
    declare_parameters();
    prm.parse_input(parameter_file);

    // read the refrectiv index list
    interpret_refrectiv_index();
  }



  // === print parameters to console ===
  void
  ParameterReader::print_parameters()
  {
    // TODO: need some rework, not now

    prm.enter_subsection("Mesh & geometry parameters");
    const unsigned int poly_degree   = prm.get_integer("Polynomial degree");
    const unsigned int n_refinements = prm.get_integer("Number of refinements");
    prm.leave_subsection();

    prm.enter_subsection("Physical constants");
    const double omega = prm.get_double("omega");
    const double c     = prm.get_double("c");
    prm.leave_subsection();

    deallog << "Using the following parameters: " << std::endl;
    deallog << "\tMesh & geometry parameters: " << std::endl;
    deallog << "\t\tPolynomial degree:\t" << poly_degree << std::endl;
    deallog << "\t\tNumber of refinements:\t" << n_refinements << std::endl;
    deallog << "\tPhysical constants: " << std::endl;
    deallog << "\t\tomega:\t" << omega << std::endl;
    deallog << "\t\tc:\t" << c << std::endl;
  }



  unsigned int
  ParameterReader::get_integer(const std::string &entry_subsection_path,
                               const std::string &entry_string) const
  {
    std::vector<std::string> sub_sec;
    sub_sec.push_back(entry_subsection_path);

    return prm.get_integer(sub_sec, entry_string);
  }



  double
  ParameterReader::get_double(const std::string &entry_subsection_path,
                              const std::string &entry_string) const
  {
    std::vector<std::string> sub_sec;
    sub_sec.push_back(entry_subsection_path);

    return prm.get_double(sub_sec, entry_string);
  }



  std::string
  ParameterReader::get_string(const std::string &entry_subsection_path,
                              const std::string &entry_string) const
  {
    std::vector<std::string> sub_sec;
    sub_sec.push_back(entry_subsection_path);

    return prm.get(sub_sec, entry_string);
  }



  std::complex<double>
  ParameterReader::get_refrective_index(const unsigned int i) const
  {
    Assert(i < refrectiv_index_list.size(),
           ExcMessage(
             "You requested a refrective index, for a material id that"
             "does not exists. This can have several resons:\n"
             "\t1. The input string in the option file for the refrectiv"
             "index has the wrong format.\n"
             "\t2. You requesting a material ID, that does not exists, "
             "which is a sign, that your input geometry and input for the "
             "refrective index do not correspond."));

    return refrectiv_index_list[i];
  }



  double
  ParameterReader::get_wavenumber() const
  {
    std::string sub_sec      = "Physical constants";
    std::string entry_string = "lambda";

    const double omega = get_double(sub_sec, entry_string);
    return (2 * numbers::PI) / omega;
  }



  std::complex<double>
  ParameterReader::get_wavenumber(const unsigned int i) const
  {
    Assert(i < refrectiv_index_list.size(),
           ExcMessage(
             "You requested a refrective index, for a material id that"
             "does not exists. This can have several resons:\n"
             "\t1. The input string in the option file for the refrectiv"
             "index has the wrong format.\n"
             "\t2. You requesting a material ID, that does not exists, "
             "which is a sign, that your input geometry and input for the "
             "refrective index do not correspond."));

    std::string sub_sec      = "Physical constants";
    std::string entry_string = "lambda";

    const double omega = get_double(sub_sec, entry_string);
    return refrectiv_index_list[i] * (2 * numbers::PI) / omega;
  }

} // namespace KirasFM

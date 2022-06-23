#ifndef WRITE_DATA_TABLE_H
#define WRITE_DATA_TABLE_H

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>

namespace KirasFM {
using namespace dealii;

class WriteDataTable {
  public:
    // constructor
    WriteDataTable();

    void add_entery(const std::string name, const double in);
    void add_entery(const std::string name, const int in);
    void add_entery(const std::string name, const std::string &in);

    int write(const std::string &file_name, const std::string &format = "csv");

  private:
    MPI_Comm                    mpi_communicator;

    // internal function
    void print_csv(const std::string &file_name, bool headline = true);
    void print_hr(const std::string &file_name, bool headline = true);
    void print_cml();

    std::vector< std::string >  intern_string;

    std::vector< std::string >  string_names;

};

} // namespace: KirasFM
#endif

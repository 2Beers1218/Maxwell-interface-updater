#ifndef PARAMETER_HANDLER_H 
#define PARAMETER_HANDLER_H 

// === deal.II includes ===
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/data_out.h>
// === C++ includes ===
#include <iostream>

namespace KirasFM {
  using namespace dealii;

  class ParameterReader : public Subscriptor {
    public: 
      ParameterReader(ParameterHandler &);

      void read_parameters(const std::string &);

      void print_parameters();

      // return functions
      unsigned int get_integer( const std::string &entry_subsection_path, 
                       const std::string &entry_string ) const;

      double get_double( const std::string &entry_subsection_path, 
                         const std::string &entry_string ) const;

      std::string get_string( const std::string &entry_subsection_path, 
                              const std::string &entry_string ) const;

      double get_refrective_index(const unsigned int i) const;

      double get_wavenumber() const;


    private:
      // internal functions
      void declare_parameters();
      void interpret_refrectiv_index( );

      // the PRM itself
      ParameterHandler    &prm;

      // the list of refrectiv indecies
      std::vector< double > refrectiv_index_list;
      
  };


} // namespace KirasFM

#endif

#include <write_data_table.h>

namespace KirasFM {
  using namespace dealii;

  // help functions:
  std::string fill_with_witespaces(std::string in, unsigned int length) {
    std::string out = in;
    for( unsigned int i = out.length(); i <= length; i++ ) {
      out += " ";
    }

    return out;
  }

  // constructor
  WriteDataTable::WriteDataTable():
    mpi_communicator(MPI_COMM_WORLD)
  {}

  void WriteDataTable::add_entery(const std::string name, const double in) {
    if(Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {
      intern_string.push_back(std::to_string(in));
      string_names.push_back(name);
    }
  }

  void WriteDataTable::add_entery(const std::string name, const int in) {
    if(Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {
      intern_string.push_back(std::to_string(in));
      string_names.push_back(name);
    }
  }

  void WriteDataTable::add_entery(const std::string name, const std::string &in) {
    if(Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {
      intern_string.push_back(in);
      string_names.push_back(name);
    }
  }

  int WriteDataTable::write(
    const std::string &file_name, 
    const std::string &format
  ) {

    if(Utilities::MPI::this_mpi_process(mpi_communicator) != 0)
      return 0;

    Assert( 
      intern_string.size() == 0,
      ExcMessage("You tried to print an empty list, try to add some conente first")
    );
   
    if ( format.compare("csv") == 0) {
      print_csv(file_name);
    } else if ( format.compare("cml") == 0 ) {
      print_cml();
    } else if ( format.compare("hr") == 0 ) {
      print_hr(file_name);
    }

    else {
      Assert(
        false,
        ExcMessage("Unkown file format:" + format + "Valid formats are: csv")
      );
    }

    return 1;
  }

  void WriteDataTable::print_csv(const std::string &file_name, bool headline) {
    std::string content_string;
    std::string headline_string;

    Assert(
      content_string.size() == 0,
      ExcMessage("You tried to print an empty list")
    )

    for ( std::string s : intern_string ) {
      content_string += s;
      content_string += ",";
    }
    // remove the "," at the end of the string:
    content_string.pop_back();

    content_string += "\n";

    for ( std::string s : string_names ) {
      headline_string += s;
      headline_string += ",";
    }
    // remove the "," at the end of the string:
    headline_string.pop_back();

    headline_string += "\n";

    // check if the file already exists
    std::ifstream file;
    file.open(file_name);
    bool file_empty = ( file.peek() == std::ifstream::traits_type::eof() );
    file.close();

    // write the output into the file
    std::ofstream file_write;
    file_write.open(file_name, std::ios::app);
    if( file_empty && headline ) {
      file_write << headline_string;
    }
    file_write << content_string;
    file_write.close();
  }

  void WriteDataTable::print_hr(const std::string &file_name, bool headline) {

    std::vector<unsigned int> length_list;

    // get the longest arguments
    for ( unsigned int i = 0; i < intern_string.size(); i++ ) {
      if ( intern_string[i].length() < string_names[i].length() ) {
        length_list.push_back( string_names[i].length() );
      } else {
        length_list.push_back( intern_string[i].length() );
      }
    }

    std::string seperator = "+";
    std::string headline_seperator = "+";
    for( unsigned int i : length_list ) {
      for( unsigned int j = 0; j <= i; j++ ) {
        seperator += "-";
        headline_seperator += "=";
      }
      seperator += "+";
      headline_seperator += "+";
    }

    std::string out;

    // check if the file already exists
    std::ifstream file;
    file.open(file_name);
    bool file_empty = ( file.peek() == std::ifstream::traits_type::eof() );
    file.close();

    // write the output into the file
    std::ofstream file_write;
    file_write.open(file_name, std::ios::app);
    if( file_empty && headline ) {
      out += seperator;
      out += "\n";

      // print the headline
      out += "|";
      for( unsigned int i = 0; i < string_names.size(); i++) {
        out += fill_with_witespaces(string_names[i], length_list[i]);
        out += "|";
      }
      out += "\n";
      out += headline_seperator;
      out += "\n";
    }
    // print the content
    out += "|";
    for( unsigned int i = 0; i < string_names.size(); i++) {
      out += fill_with_witespaces(intern_string[i], length_list[i]); 
      out += "|";
    }
    out += "\n";
    out += seperator;
    out += "\n";

    file_write << out;
    file_write.close();
  }

  void WriteDataTable::print_cml() {

    unsigned int max_length_string = 0;
    for( std::string s : intern_string ) {
      if ( s.length() > max_length_string )
        max_length_string = s.size();
    }

    // get the longest arguments
    unsigned int max_length_name = 0;
    for( std::string s : string_names ) {
      if ( s.length() > max_length_name )
        max_length_name = s.size();
    }

    std::string seperator = "+";
    for( unsigned int i = 0; i <= max_length_name; i++)
      seperator += "-";
    seperator += "+";

    for( unsigned int i = 0; i < max_length_name; i++)
      seperator += "-";
    seperator += "+";

    std::cout << seperator << std::endl;
    for ( unsigned int i = 0; i < intern_string.size(); i++ ) {
      std::cout << "|" 
                << fill_with_witespaces(string_names[i], max_length_name) 
                << "|"
                << fill_with_witespaces(intern_string[i], max_length_string) 
                << "|" << std::endl;
      std::cout << seperator << std::endl;
    }


  }

}//  namespace: KirasFM

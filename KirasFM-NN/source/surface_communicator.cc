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
#include <surface_communicator.h>

namespace KirasFM
{
  using namespace dealii;

  template <int dim>
  SurfaceCommunicator<dim>::SurfaceCommunicator()
  {}

  template <int dim>
  SurfaceCommunicator<dim>::SurfaceCommunicator(unsigned int n_domains_)
  {
    n_domains = n_domains_;

    // initialize value
    std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> tmp_value;

    std::vector<std::vector<std::complex<double>>> tmp_curl;

    for (unsigned int i = 0; i < n_domains * n_domains; i++)
      {
        value_data.push_back(tmp_value);
        curl_data.push_back(tmp_curl);
      }
  }

  // Constructor: Read from file
  template <int dim>
  SurfaceCommunicator<dim>::SurfaceCommunicator(unsigned int n_domains_,
                                                std::string  filename,
                                                unsigned int n_face_q_points,
                                                unsigned int n_rows)
  {
    n_domains = n_domains_;

    // initialize value
    std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> tmp_value;

    for (unsigned int i = 0; i < n_domains * n_domains; i++)
      {
        value_data.push_back(tmp_value);
      }

    // --- Vorbereitung ---
    // Zum wieder einlesen der Datei benoetigst du:
    //  - Wie viele Zahlen stehen in einer Reihe, das ist besimmt durch:
    //  n_face_q_points
    //  - Wie viele Zeilen hat die Matrix: ich nenne die Variable hier: n_rows
    //  (beim schreiben der Datei, nennst du die Variable subfaces

    // Erstelle einen Vektor von Tensor<1, dim, std::complex>,
    // Remark: A rank-1 tensor is a vector with dim components
    // const unsigned int n_face_q_points = 8; // TODO
    std::vector<Tensor<1, dim, std::complex<double>>> update_value(
      n_face_q_points);

    // length n_rows
    // const unsigned int n_rows = 16; // TODO
    std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> g_value(
      n_rows);

    // --- Einlesen ---
    // lese die csv als matrix ein, z.B.: als
    // std::vector<std::vector<std::complex<double>>> csv_input;
    for (unsigned int i = 0; i < n_domains; i++)
      {
        for (unsigned int j = 0; j < n_domains; j++)
          {
            int diff_i_j = i - j;
            if (std::abs(diff_i_j) == 1)
              {
                std::string filename_i_j = filename + std::to_string(i) + "_" +
                                           std::to_string(j) + "_nn" + ".csv";
                // std::cout << filename_i_j << std::endl;
                std::ifstream file(filename_i_j);
                std::vector<std::vector<std::complex<double>>>
                  csv_input; // TODO

                std::string line;


                while (getline(file, line))
                  {
                    std::vector<std::complex<double>> vec;
                    int                               pos;

                    while ((pos = line.find(";")) >= 0)
                      {
                        std::string field = line.substr(0, pos);
                        line              = line.substr(pos + 1);
                        // std::cout << "field= " << field << std::endl;
                        std::istringstream   is(field);
                        std::complex<double> field_num;
                        is >> field_num;
                        // std::cout << "field_num= " << field_num << std::endl;
                        vec.push_back(field_num);
                      }

                    // for(unsigned int index = 0; index < vec.size(); index++)
                    //   std::cout << vec.at(index) << ' ';
                    csv_input.push_back(vec);
                  }

                file.close();

                // Schleife ueber die Zeilen der Matrix
                for (unsigned int row = 0; row < n_rows; row++)
                  {
                    // Schleife uber die Spalten der Matrix
                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                         q_point++)
                      {
                        // Aufbau von jedem Elemente: (Re, Im) (Re, Im) (fuer
                        // dim = 2) hier: (a, b) (c, d)
                        // TODO: Ich bin mir nicht sicher ob die Syntax so
                        // funktioniert
                        std::complex<double> ab = csv_input[row][2 * q_point];
                        std::complex<double> cd =
                          csv_input[row][2 * q_point + 1];
                        update_value[q_point] =
                          Tensor<1, dim, std::complex<double>>({ab, cd});

                      } // rof: q_point

                    g_value[row] = update_value;
                  } // rof: row

                // --- Schreibe die eingelesenen Werte in value_data ---
                value_data[(i * n_domains) + j] = g_value;
              } // fi

          } // rof: j
      } // rof: i
  }

  // Copy Constructor (Rule of 3: 2/3)
  template <int dim>
  SurfaceCommunicator<dim>::SurfaceCommunicator(
    const SurfaceCommunicator<dim> &copy)
  {
    n_domains = copy.n_domains;
    n_faces   = copy.n_faces;

    value_data = copy.value_data;
    curl_data  = copy.curl_data;
  }

  // return functions
  template <int dim>
  std::vector<std::vector<Tensor<1, dim, std::complex<double>>>>
  SurfaceCommunicator<dim>::value(unsigned int i, unsigned int j)
  {
    AssertIndexRange(i, n_domains);
    AssertIndexRange(j, n_domains);
    return value_data[(i * n_domains) + j];
  }

  template <int dim>
  std::vector<std::vector<std::complex<double>>>
  SurfaceCommunicator<dim>::curl(unsigned int i, unsigned int j)
  {
    AssertIndexRange(i, n_domains);
    AssertIndexRange(j, n_domains);
    return curl_data[(i * n_domains) + j];
  }

  template <int dim>
  void
  SurfaceCommunicator<dim>::value(
    std::vector<std::vector<Tensor<1, dim, std::complex<double>>>> in,
    unsigned int                                                   i,
    unsigned int                                                   j)
  {
    AssertIndexRange(i, n_domains);
    AssertIndexRange(j, n_domains);
    value_data[(i * n_domains) + j] = in;
  }

  template <int dim>
  void
  SurfaceCommunicator<dim>::curl(
    std::vector<std::vector<std::complex<double>>> in,
    unsigned int                                   i,
    unsigned int                                   j)
  {
    AssertIndexRange(i, n_domains);
    AssertIndexRange(j, n_domains);
    curl_data[(i * n_domains) + j] = in;
  }


  template <int dim>
  void
  SurfaceCommunicator<dim>::update(SurfaceCommunicator<dim> sc, unsigned int i)
  {
    AssertIndexRange(i, n_domains);
    for (unsigned int j = 0; j < n_domains; j++)
      {
        value_data[(i * n_domains) + j] = sc.value(i, j);
        curl_data[(i * n_domains) + j]  = sc.curl(i, j);
      }
  }


  // Assaignment operator (Rule of 3: 3/3)
  template <int dim>
  SurfaceCommunicator<dim> &
  SurfaceCommunicator<dim>::operator=(const SurfaceCommunicator<dim> &copy)
  {
    if (this == &copy)
      return *this;

    n_domains = copy.n_domains;
    n_faces   = copy.n_faces;

    value_data = copy.value_data;
    curl_data  = copy.curl_data;

    return *this;
  }

  template <int dim>
  void
  SurfaceCommunicator<dim>::to_file(std::string  name,
                                    unsigned int i,
                                    unsigned int j)
  {
    AssertIndexRange(i, n_domains);
    AssertIndexRange(j, n_domains);

    unsigned int subfaces = value_data[(i * n_domains) + j].size();
    if (subfaces != 0)
      {
        unsigned int n_face_q_points =
          value_data[(i * n_domains) + j][0].size();
        std::ofstream file;
        file.open(name);
        // loop over all subfaces
        for (unsigned int k = 0; k < subfaces; k++)
          {
            // loop over all quadrature points
            for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++)
              {
                file << value_data[(i * n_domains) + j][k][q_point];
                if (q_point != n_face_q_points - 1)
                  file << "; ";
              }
            if (k != n_faces - 1)
              file << "\n";
          }
        file.close();

      } // fi
  }


  template class SurfaceCommunicator<2>;
} // namespace KirasFM

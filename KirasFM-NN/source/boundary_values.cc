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

#include <boundary_values.h>

namespace KirasFM
{
  template <>
  void
  DirichletBoundaryValues<2>::vector_value(const Point<2> &p,
                                           Vector<double> &values) const
  {
    //values(0) = intensity_real *
    //            std::exp(-(std::pow(p(0) - center_real, 2) * width_real));
    //values(1) = 0.0;
    //values(2) = intensity_imag *
    //            std::exp(-(std::pow(p(0) - center_imag, 2) * width_imag));
    //values(3) = 0.0;

    // Boundary conditions for the control set:
    values(0) = intensity_real * ( sin(numbers::PI * p(0) * width_real) );
    values(1) = 0.0;
    values(2) = intensity_imag * ( sin(numbers::PI * p(0) * width_real) );
    values(3) = 0.0;
  }

  template <>
  void
  DirichletBoundaryValues<3>::vector_value(const Point<3> & /*p*/,
                                           Vector<double> & /*values*/
  ) const
  {
    Assert(false, ExcNotImplemented());
  }
} // namespace KirasFM

#ifndef BOUNDARY_VALUES_H 
#define BOUNDARY_VALUES_H 

// === deal.II includes ===
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

// === C++ includes ===
#include <iostream>
#include <cmath>
#include <complex>

namespace KirasFM {
  using namespace dealii;

// === DirichletBoundaryValues ===
/*
 * DirichletBoundaryValues defines the function E_{inc} in the
 * equation:
 *      trace( E ) = trace( E_{inc} ) on \Gamma_{inc}
 */
  template <int dim>
  class DirichletBoundaryValues : public Function<dim> {
    public:
      DirichletBoundaryValues():
        Function<dim>(2 * dim)
      {}

      virtual void vector_value(
        const Point<dim> & p, 
        Vector<double> &values
      ) const ;

      virtual void vector_value_list(
        const std::vector<Point<dim>> &points, 
        std::vector<Vector<double>> &value_list
      ) const override {
        Assert(
          value_list.size() == points.size(), 
          ExcDimensionMismatch(value_list.size(), 
          points.size())
        );

        for (unsigned int p = 0; p < points.size(); p++) {
          DirichletBoundaryValues<dim>::vector_value(
            points[p], 
            value_list[p]
          );
        }
      }
  };



// === CurlRHS ===
/* 
 * this defines the function f(x) in:
 *      curl ( curl ( E ) ) - \omega^2 E = f( x ) on \Omega
 */
  template <int dim>
  class CurlRHS : public Function<dim> {
    public:
      CurlRHS():
        Function<dim>(2 * dim)
      {}

      virtual void vector_value(
        const Point<dim> & p, 
        Vector<double> &values
      ) const ;

      virtual void vector_value_list(
        const std::vector<Point<dim>> &points, 
        std::vector<Vector<double>> &value_list
      ) const override {
        Assert(
          value_list.size() == points.size(), 
          ExcDimensionMismatch(value_list.size(), 
          points.size())
        );

        for (unsigned int p = 0; p < points.size(); p++) {
          CurlRHS<dim>::vector_value(
            points[p], 
            value_list[p]
          );
        }
      }
  };


} // namespace KirasFM 
#endif

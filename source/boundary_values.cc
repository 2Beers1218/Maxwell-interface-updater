#include <boundary_values.h>

namespace KirasFM {
  const double PI  = 3.141592653589793;
  const double PI2 = 9.869604401089359;

  // here we are using quite artificial boundary values,
  // the benefit is, that we know the exact solution for the
  // electric field E, i.e. this is also the exact solution 
  // to the partial differential equation we aim to solve
  template<>
  void DirichletBoundaryValues<2>::vector_value (
    const Point<2> &p,
    Vector<double> &values
  ) const {
    // standart RHS (analythical solution is not known)
    //values(0) = exp(- std::pow( p(0) - 0.5, 2 ) / 0.003 );
    //values(1) = 0.0;
    //values(2) = 0.0;
    //values(3) = 0.0;

    //values(0) = exp(- std::pow( p(0) - 0.7, 2 ) / 0.008 );
    //values(1) = 0.0;
    //values(2) = 0.0;
    //values(3) = 0.0;

    //values(0) = exp(- std::pow( p(0) - 0.2, 2 ) / 0.002 );
    //values(1) = 0.0;
    //values(2) = 1.0;
    //values(3) = 0.0;

    //values(0) = exp(- std::pow( p(0) - 0.7, 2 ) / 0.003 );
    //values(1) = 0.0;
    //values(2) = 1.0;
    //values(3) = 0.0;

    //values(0) = exp(- std::pow( p(0) - 0.8, 2 ) / 0.003 );
    //values(1) = 0.0;
    //values(2) = sin( PI2 * p[0] );
    //values(3) = 0.0;

    //values(0) = exp(- std::pow( p(0) - 0.5, 2 ) / 0.003 );
    //values(1) = 0.0;
    //values(2) = cos( PI2 * p[0] );
    //values(3) = 0.0;

    // analyticall test case where the RHS is zero
    //const double PI2 = 9.869604401089359;
    values(0) = cos( PI2 * p[1] );
    values(1) = sin( PI2 * p[0] );
    values(2) = sin( PI2 * p[1] );
    values(3) = 0.5 * cos( PI2 * p[0] );

    //values(0) = sin( PI2 * p[0] );
    //values(1) = sin( PI2 * p[0] );
    //values(2) = sin( PI2 * p[1] );
    //values(3) = 0.5 * cos( PI2 * p[0] );

    //values(0) = sin( PI2 * p[0] );
    //values(1) = sin( PI2 * p[0] );
    //values(2) = sin( PI2 * p[0] );
    //values(3) = 0.5 * cos( PI2 * p[0] );

    //values(0) = sin( PI2 * p[1] );
    //values(1) = sin( PI2 * p[0] );
    //values(2) = sin( PI2 * p[0] );
    //values(3) = 0.5 * cos( PI2 * p[0] );

    //values(0) = cos( PI2 * p[1] );
    //values(1) = sin( PI2 * p[0] );
    //values(2) = cos( PI2 * p[0] );
    //values(3) = 0.5 * cos( PI2 * p[0] );

    //values(0) = cos( PI2 * p[0] );
    //values(1) = sin( PI2 * p[0] );
    //values(2) = cos( PI2 * p[0] );
    //values(3) = 0.5 * cos( PI2 * p[0] );

    //values(0) = cos( PI2 * p[0] );
    //values(1) = sin( PI2 * p[0] );
    //values(2) = cos( PI2 * p[1] );
    //values(3) = 0.5 * cos( PI2 * p[0] );

    //values(0) = cos( PI2 * p[1] );
    //values(1) = sin( PI2 * p[0] );
    //values(2) = cos( PI2 * p[1] );
    //values(3) = 0.5 * cos( PI2 * p[0] );
    
    // analyticall test case where the RHS is not zero
    //values(0) = cos( PI * p(0) ) * sin( PI * p(1) );
    //values(1) = 0.0;
    //values(2) = 0.0;
    //values(3) = 0.0;
  }

  template<>
  void DirichletBoundaryValues<3>::vector_value (
    const Point<3> &p,
    Vector<double> &values
  ) const {
    double center = 8;
    values(0) = exp( - (std::pow(p(0) - center, 2) + std::pow(p(1) - center, 2)) / 7.0 );
    values(1) = 0.0; //exp( - (std::pow(p(0) - center, 2) + std::pow(p(1) - center, 2)) / 5.0 );
    values(2) = exp( - (std::pow(p(0) - center, 2) + std::pow(p(1) - center, 2)) / 7.0 );
    values(3) = exp( - (std::pow(p(0) - center, 2) + std::pow(p(1) - center, 2)) / 7.0 );
    values(4) = 0.0; //exp( - (std::pow(p(0) - center, 2) + std::pow(p(1) - center, 2)) / 5.0 );
    values(5) = exp( - (std::pow(p(0) - center, 2) + std::pow(p(1) - center, 2)) / 7.0 );


    // 3D anallytical test case
    //std::complex<double> imag_i(0.0, 1.0);
    //std::complex<double> tmp1 = ( imag_i * exp(imag_i * p(1) * PI2) ) + ( exp(imag_i * p(2) * PI2) );
    //std::complex<double> tmp2 = ( imag_i * exp(imag_i * p(0) * PI2) ) + ( exp(imag_i * p(2) * PI2) );
    //std::complex<double> tmp3 = ( imag_i * exp(imag_i * p(0) * PI2) ) + ( exp(imag_i * p(1) * PI2) );
    //values(0) = tmp1.real();
    //values(1) = tmp2.real();
    //values(2) = tmp3.real();
    //values(3) = tmp1.imag(); 
    //values(4) = tmp2.imag(); 
    //values(5) = tmp3.imag(); 
  }

  template<>
  void CurlRHS<2>::vector_value (
    [[maybe_unused]] const Point<2> &p,
    Vector<double> &values
  ) const {
    values(0) = 0.0;
    //values(1) = - PI2 * sin( PI * p(0) ) * cos( PI * p(1) ); // needed for the analyticall testcase where the RHS is not zero
    values(1) = 0.0;
    values(2) = 0.0;
    values(3) = 0.0;
  }

  template<>
  void CurlRHS<3>::vector_value (
    [[maybe_unused]] const Point<3> &p,
    Vector<double> &values
  ) const {
    values(0) = 0.0;
    values(1) = 0.0; //- PI2 * sin( PI * p(0) ) * cos( PI * p(1) );
    values(2) = 0.0;
    values(3) = 0.0;
    values(4) = 0.0;
    values(5) = 0.0;
  }




} // namespace KirasFM

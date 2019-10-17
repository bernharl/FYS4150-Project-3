# include "catch.hpp"
# include "integration.cpp"



TEST_CASE("Comparing Gauss-Laguerre quadrature to known result")
{ 
  double N = 20;
  double alpha = 2;
  double tolerance = 1e-7;
  // Correction due to change of variables to spherical coordinates in integral function
  double correction_factor = (32 * pow(alpha, 5)); 
  double analytical = 64 * PI * PI;
  double numerical;
  numerical = gauss_quad_improved( N, alpha, &debug_integrand) * correction_factor;

  REQUIRE(std::fabs(numerical - analytical) == Approx(0).epsilon(tolerance));
}

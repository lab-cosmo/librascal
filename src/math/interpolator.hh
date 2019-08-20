#ifndef SRC_MATH_INTERPOLATOR_HH_
#define SRC_MATH_INTERPOLATOR_HH_

#include <functional>
#include <forward_list>
#include <iostream>
#include <limits>
#include <cassert>
#include "math_utils.hh"

namespace rascal {
  namespace math {
    using Vector_Ref = typename Eigen::Ref<const Vector_t>;
    using Matrix_Ref = typename Eigen::Ref<const Matrix_t>;

    enum class GridType_t {Uniform};
    enum class RefinementMethod_t {Exponential, Linear, Adaptive};

    template <GridType_t Type, RefinementMethod_t Method>
    struct GridRational {
    };

    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive> {
      GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>() : grid_meshes{},grid_meshes_error{}, grid_size{0}, max_it{} {}

      constexpr static GridType_t GridType { GridType_t::Uniform };
      constexpr static RefinementMethod_t RefinementMethod { RefinementMethod_t::Adaptive };
      // These are possible improvements for this grid rational:
      // - for this grid rational it makes sense to save the test grid,
      // because it is one step ahead and we can add more points 
      // - make hyperparameter k highest error
      // - make grid rational init compute, increase_finness() and then remove x1,x2,fineness
      Vector_t compute_grid(double x1, double x2, int fineness) {
        if (fineness == 0) {
          this->grid_meshes = {x1, x2};
          this->grid_size = 2;
          return this->grid_from_meshes();
        }
        this->refine();
        return this->grid_from_meshes();
      }

      void refine() {
        auto next_it{this->max_it};
        next_it++;
        double cur{*this->max_it};
        double mid_point{cur+(*next_it-cur)/2};
        this->grid_meshes.emplace_after(this->max_it, mid_point);
        this->grid_size++;
      }

      Vector_t grid_from_meshes() {
        Vector_t grid = Vector_t::Zero(this->grid_size);
        int i{0};
        for (auto it=this->grid_meshes.begin(); it!=this->grid_meshes.end(); ++it) {
          grid(i) = (*it);
          i++;
        }
        return grid;
      }
      
      Vector_t compute_test_grid(double, double, int) {
        Vector_t test_grid = Vector_t::Zero(this->grid_size-1);
        int i{0};
        auto next_it{this->grid_meshes.begin()};
        next_it++;
        // if auto does not work use: std::forward_list<double>::iterator  
        for (auto it=this->grid_meshes.begin(); next_it!=this->grid_meshes.end(); ++it) {
          double cur{*it};
          double mid_point{cur+(*next_it-cur)/2};
          test_grid(i) = mid_point;

          next_it++;
          i++;
        }
        return test_grid;
      }

      void update_errors(Vector_Ref error_grid) {
        if (error_grid.size() != this->grid_size -1) {
          std::runtime_error("Gridsize does not match with error grid");
        }
        double max{0};
        auto it{this->grid_meshes.begin()};
        for (int i{0}; i < error_grid.size(); i++) {
          if (std::abs(error_grid(i)) > max) {
            max = std::abs(error_grid(i));          
            this->max_it = it;
          }
          it++;
        }
      }

      // We store the grid sorted in two ways:
      // - for the refinement key which is the error
      // - for constructing of the grid x1 of each mesh
      std::forward_list<double> grid_meshes;
      // One could make a tree with meshes as nodes (x1, error), each mesh
      // only owns the x1 and error. 
      std::forward_list<double> grid_meshes_error;
      int grid_size;
      std::forward_list<double>::iterator max_it;
    };

    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::Linear> {
      constexpr static GridType_t GridType {  GridType_t::Uniform };
      constexpr static RefinementMethod_t RefinementMethod { RefinementMethod_t::Linear };

      Vector_t compute_grid(double x1, double x2, int fineness) {
        double nb_grid_points = fineness+2;
        return Vector_t::LinSpaced(nb_grid_points, x1, x2);
      }
      Vector_t compute_test_grid(double x1, double x2, int fineness) {
        double nb_grid_points = fineness+2;
        // half of a step size in the grid to get the inner grid points
        // step size = (x2-x1)/(nb_grid_points-1)
        double offset{(x2-x1)/(2*(nb_grid_points-1))};
        return Vector_t::LinSpaced(nb_grid_points-1, x1+offset, x2-offset);
      }

      void update_errors(Vector_Ref) {} 
      int grid_size{0};
      // One could extend this with a configurable linear slope factor
    };

    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential> {
      constexpr static GridType_t GridType {  GridType_t::Uniform };
      constexpr static RefinementMethod_t RefinementMethod { RefinementMethod_t::Exponential };

      Vector_t compute_grid(double x1, double x2, int fineness) {
        double nb_grid_points = 2 << fineness;
        return Vector_t::LinSpaced(nb_grid_points, x1, x2);
      }
      Vector_t compute_test_grid(double x1, double x2, int fineness) {
        double nb_grid_points = 2 << fineness;
        // half of a step size in the grid to get the inner grid points
        // step size = (x2-x1)/(nb_grid_points-1)
        double offset{(x2-x1)/(2*(nb_grid_points-1))};
        return Vector_t::LinSpaced(nb_grid_points-1, x1+offset, x2-offset);
      }

      void update_errors(Vector_Ref) {} 
      int grid_size{0};
      // One could extend this with a configurable exponential basis factor
      // default would be 2
    };

    enum class InterpolationMethod_t {CubicSpline, CubicSplineVectorized};

    template <InterpolationMethod_t Type>
    struct InterpolationMethod{};

    // An implementation of the cubic spline method, the implementation is adapted from "Numerical Recipes" [1] and optimized for uniform grids.
    //
    // [1] Press, William H., et al. Numerical recipes 3rd edition: The art of scientific computing. Cambridge university press, 2007.
    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSpline> {
     public:
      constexpr static InterpolationMethod_t Method { InterpolationMethod_t::CubicSpline };

      InterpolationMethod<InterpolationMethod_t::CubicSpline>() {}

      void initialize(const Vector_Ref & grid,
          const Vector_Ref & evaluated_grid){
        this->h = grid(1) - grid(0);
        this->h_sq_6 = this->h*this->h/6.0;
        this->compute_second_derivative(evaluated_grid);
      }

      void initialize(const Vector_Ref & grid,
          const Vector_Ref & evaluated_grid, double yd1, double ydn){
        this->h = grid(1) - grid(0);
        this->h_sq_6 = this->h*this->h/6.0;
        this->compute_second_derivative(evaluated_grid, yd1, ydn);
      }

      inline double interpolate(const Vector_Ref & grid,
          const Vector_Ref & evaluated_grid,
          double x, int nearest_grid_index_to_x) {
        return this->rawinterp(grid, evaluated_grid,
            nearest_grid_index_to_x, x);
      }
      inline double interpolate_derivative(const Vector_Ref & grid,
          const Vector_Ref & evaluated_grid,
          double x, int nearest_grid_index_to_x) {
        return this->rawinterp_derivative(grid, evaluated_grid,
            nearest_grid_index_to_x, x);
      }

     private:
      // The two functions below implement the tridiagonal matrix algorithm
      // to solve the linear system in
      // http://mathworld.wolfram.com/CubicSpline.html
      // under different boundary conditions. We describe with s the spline function and with y the targeted function for interpolation.
      // 
      // Reference:
      // https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
      
      // clamped boundary conditions when first derivative for boundaries is known
      // s'(x_0) = y'(x_0), s'(x_n) = y'(x_n)
      inline void compute_second_derivative(const Vector_Ref & yv, double yd1, double ydn) {
        // Bad xa input to routine splint
        assert (this->h != 0.0);
        int n{static_cast<int>(yv.size())};
        Vector_t y2 = Vector_t::Zero(n);
        Vector_t u = Vector_t::Zero(n);
        double p;

        // forward sweeping 
        y2(0) = -0.5;
        u(0)=(3.0/this->h)*((yv(1)-yv(0))/this->h -yd1);
        for (int i{1}; i<n-1; i++) {
          p=0.5*y2(i-1)+2.0;
          y2(i)=-0.5/p;
          u(i)=(yv(i+1)-2*yv(i)+yv(i-1))/this->h;
          u(i)=(3.0*u(i)/this->h -0.5*u(i-1))/p;
        }

        // back substitution
        p=0.5; // qn=0.5 but we just reuse p because it is not needed anyway
        u(n-1) = (3.0/this->h)*(ydn-(yv(n-1)-yv(n-2))/this->h);
        y2(n-1)=(u(n-1)-p*u(n-2))/(p*y2(n-2)+1.0);
        for (int k{n-2};k>=0;k--) {
          y2(k)=y2(k)*y2(k+1)+u(k);
        }
        this->second_derivative_h_sq_6 = y2*this->h_sq_6 ;
      }

      // natural/simple boundary conditions s''(x_0) = s''(x_n) = 0
      inline void compute_second_derivative(const Vector_Ref & yv) {
        // Bad xa input to routine splint
        assert (this->h != 0.0);
        int n{static_cast<int>(yv.size())};
        Vector_t y2 = Vector_t::Zero(n);
        Vector_t u = Vector_t::Zero(n);
        double p;

        // forward sweeping 
        y2(0) = 0.0;
        u(0) = 0.0;
        for (int i{1}; i<n-1; i++) {
          p=0.5*y2(i-1)+2.0;
          y2(i)=-0.5/p;
          u(i)=(yv(i+1)-2*yv(i)+yv(i-1))/this->h;
          u(i)=(3.0*u(i)/this->h -0.5*u(i-1))/p;
        }

        // back substitution
        y2(n-1) = 0.0; // (u(n-1)-p*u(n-2))/(p*y2(n-2)+1.0);
        for (int k{n-2};k>0;k--) {
          y2(k)=y2(k)*y2(k+1)+u(k);
        }
        y2(0)= 0.0; // y2(0)*y2(1)+u(0);
        this->second_derivative_h_sq_6 = y2*this->h_sq_6 ;
      }

      inline double rawinterp(const Vector_Ref & xx, const Vector_Ref & yy,
          const int & j1, const double & x) {
        assert (this->h != 0.0);
        const int klo{j1}, khi{j1+1};
        // a+b=1
        const double a{(xx(khi)-x)/this->h};
        double b{1-a};
        return a*( yy(klo)+(a*a-1)*this->second_derivative_h_sq_6(klo) ) +
            b*( yy(khi)+(b*b-1)*this->second_derivative_h_sq_6(khi) );
      }

      inline double rawinterp_derivative(const Vector_Ref & xx, const Vector_Ref & yy,
          const int & j1, const double & x) {
        assert (this->h != 0.0);
        const int klo{j1}, khi{j1+1};
        // It is a+b=1
        const double a{(xx(khi)-x)/this->h};
        const double b{1-a};
        return ( yy(khi)-yy(klo) - (3*a*a-1) *this->second_derivative_h_sq_6(klo) + (3*b*b-1)*this->second_derivative_h_sq_6(khi) ) / this->h;
      }

      // grid step size
      double h{0};
      // h*h/6.0
      double h_sq_6{0};
      // This term is the solution vector of the linear system of the tridiagonal toeplitz matrix derived by the conditions of cubic spline. The linear system can be found in  http://mathworld.wolfram.com/CubicSpline.html
      // The term is not the second derivative, but called in numerical recipes like this. Because of consistency we keep the naming. In addition we already multiply them with coefficient
      // y2 *h*h/6
      Vector_t second_derivative_h_sq_6{};
    };


    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized> {
     public:
      constexpr static InterpolationMethod_t Method { InterpolationMethod_t::CubicSplineVectorized };

      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>() {}

      void initialize(const Vector_Ref & grid,
          const Matrix_Ref & evaluated_grid){
        this->h = grid(1) - grid(0);
        this->h_sq_6 = this->h*this->h/6.0;
        this->compute_second_derivative(evaluated_grid);
      }

      inline Vector_t interpolate(const Vector_Ref & grid,
          const Matrix_Ref & evaluated_grid,
          double x, int nearest_grid_index_to_x) {
        return this->rawinterp(grid, evaluated_grid,
            nearest_grid_index_to_x, x);
      }
      inline Vector_t interpolate_derivative(const Vector_Ref & grid,
          const Matrix_Ref & evaluated_grid,
          double x, int nearest_grid_index_to_x) {
        return this->rawinterp_derivative(grid, evaluated_grid,
            nearest_grid_index_to_x, x);
      }
     private:
      inline void compute_second_derivative(const Matrix_Ref & yv) {
        int n{static_cast<int>(yv.rows())};
        Matrix_t y2 = Matrix_t::Zero(n,yv.cols());
        Matrix_t u = Matrix_t::Zero(n,yv.cols());
        const double sig=0; //wtf this is more accurate than the correct value 0.5 TODO(alex)
        Vector_t p = Vector_t::Zero(n);
        y2.row(0) = Vector_t::Zero(yv.cols());
        u.row(0) = Vector_t::Zero(yv.cols());
        for (int i{1}; i<n-1; i++) {
          p=sig*y2.row(i-1).array()+2.0;
          y2.row(i) = (sig-1.0)/p.array();
          // noalias keework prevents storage of temporaries
          u.row(i).noalias() = ( 3.*( yv.row(i+1).array()-2*yv.row(i).array()+yv.row(i-1).array() )/(this->h*this->h) ).matrix();
          u.row(i).array() -= sig*u.row(i-1).array();
          u.row(i).array() /= p.array();
        }
        y2.row(n-1) = Vector_t::Zero(y2.cols());
        for (int k{n-2};k>0;k--) {
          y2.row(k) = y2.row(k).array()*y2.row(k+1).array()+u.row(k).array();
        }
        y2.row(0) = Vector_t::Zero(y2.cols());
        //this->second_derivative = y2;
        this->second_derivative_h_sq_6 = y2*this->h_sq_6 ;
      }

      inline Vector_t rawinterp(const Vector_Ref & xx, const Matrix_Ref & yy,
          const int & j1, const double & x) {
        int klo{j1}, khi{j1+1};
        // Bad xa input to routine splint
        assert (h != 0.0);
        // a+b=1
        double a{(xx(khi)-x)/this->h};
        double b{1-a};
        return ( a*( yy.row(klo).array()+(a*a-1)*this->second_derivative_h_sq_6.row(klo).array() ) +
            b*( yy.row(khi).array()+(b*b-1)*this->second_derivative_h_sq_6.row(khi).array() ) ).matrix();
      }

      inline Vector_t rawinterp_derivative(const Vector_Ref & xx, const Matrix_Ref & yy,
          const int & j1, const double & x) {
        int klo{j1}, khi{j1+1};
        // Bad xa input to routine splint
        assert (h != 0.0); 
        // a+b=1
        double a{(xx(khi)-x)/this->h};
        double b{1-a};
        return ( ( yy.row(khi).array()-yy.row(klo).array() - (3*a*a-1) *this->second_derivative_h_sq_6.row(klo).array() + (3*b*b-1)*this->second_derivative_h_sq_6.row(khi).array() ) / this->h ).matrix();
      }

      // The gap between two grid points
      double h{0};
      // h*h/6.0
      double h_sq_6{0};
      // second_derivative*h*h/6
      Matrix_t second_derivative_h_sq_6{};
    };

    enum class ErrorMetric_t {Absolute, Relative};
    enum class ErrorNorm_t {Mean, Max};

    template <ErrorMetric_t Metric>
    struct ErrorMetric{};

    template <ErrorMetric_t Metric, ErrorNorm_t Norm>
    struct ErrorMethod{};

    // Works for functions with large value outside of the range [0,1] well
    template <>      
    struct ErrorMetric<ErrorMetric_t::Relative> {
      static Vector_t compute_entrywise_error(const Vector_Ref & values, const Vector_Ref & references) {
        return 2*((values - references).array().abs()/
          (std::numeric_limits< double >::min()+values.array().abs() + references.array().abs()));
      }

      static Matrix_t compute_entrywise_error(const Matrix_Ref & values, const Matrix_Ref & references) {
        return 2*((values - references).array().abs()/
          (std::numeric_limits< double >::min()+values.array().abs() + references.array().abs()));
      }
    };

    // Works for functions within the range [0,1] well
    template <>
    struct ErrorMetric<ErrorMetric_t::Absolute> {
      static Vector_t compute_entrywise_error(
          const Vector_Ref & values, const Vector_Ref & references) {
        return (values - references).array().abs();
      }

      static Matrix_t compute_entrywise_error(
          const Matrix_Ref & values, const Matrix_Ref & references) {
        return (values - references).array().abs();
      }
    };


    template <ErrorMetric_t Metric>
    struct ErrorMethod<Metric, ErrorNorm_t::Mean> : public ErrorMetric<Metric> {
      using Parent = ErrorMetric<Metric>;

      static double compute_global_error(const Vector_Ref & values, const Vector_Ref & references) {
        return Parent::compute_entrywise_error(values, references).mean();
      }

      // for the case of a (grid_size, result_size) the maxmimal mean error is used
      static double compute_global_error(const Matrix_Ref & values, const Matrix_Ref & references) {
        return Parent::compute_entrywise_error(values, references).colwise().mean().maxCoeff();
      }

      // takes a (grid_size) Vector
      static double compute_global_error(
          const Vector_Ref & errors) {
        return errors.mean();
      }

      // takes a (grid_size, result_size) matrix
      static double compute_global_error(
          const Matrix_Ref & errors) {
        return errors.colwise().mean().maxCoeff();
      }
    };

    template <ErrorMetric_t Metric>
    struct ErrorMethod<Metric, ErrorNorm_t::Max> : public ErrorMetric<Metric> {
      using Parent = ErrorMetric<Metric>;
      static double compute_global_error(const Vector_Ref & values, const Vector_Ref & references) {
        return Parent::compute_entrywise_error(values, references).maxCoeff();
      }

      // for the case of a (grid_size, result_size) the maxmimal mean error is used
      static double compute_global_error(const Matrix_Ref & values, const Matrix_Ref & references) {
        return Parent::compute_entrywise_error(values, references).maxCoeff();
      }

      // takes a (grid_size) Vector
      static double compute_global_error(
          const Vector_Ref & errors) {
        return errors.maxCoeff();
      }

      // takes a (grid_size, result_size) matrix
      static double compute_global_error(
          const Matrix_Ref & errors) {
        return errors.maxCoeff();
      }
    };


    enum class SearchMethod_t {Hunt, Locate, Uniform};

    template <SearchMethod_t Method>
    struct SearchMethod{};

    template <>
    struct SearchMethod<SearchMethod_t::Uniform> {

      constexpr static SearchMethod_t Method { SearchMethod_t::Uniform };

      void initialize(const Vector_Ref & grid){
        // nb_grid_points/unit
        this->nb_grid_points_per_unit = grid.size()/(grid(grid.size()-1)-grid(0));
        this->x1 = grid(0);
        this->grid_size = grid.size();
        this->search_size = this->grid_size - this->nb_support_points;
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref &) {
        return std::max(0,std::min(this->search_size, static_cast<int>((x-this->x1) * this->nb_grid_points_per_unit)-1));

      }
      // the number of support methods the interpolation method uses
      double nb_grid_points_per_unit{0};
      double x1{0};      
      // asserts that we only use it with CubicSpline
      size_t nb_support_points{2};
      size_t grid_size{0};
      int search_size{0};
    };


    // Improvements:
    // - make nb_support_points const
    template <>
    struct SearchMethod<SearchMethod_t::Locate> {
      constexpr static SearchMethod_t Method { SearchMethod_t::Locate };

      SearchMethod<SearchMethod_t::Locate>() : 
          nb_support_points{2} {} 

      void initialize(const Vector_Ref & ){
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref & grid) {
        return this->locate(x, grid);
      }

      int locate(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};

        int ju,jm,jl;
        // size error
        assert( not(n < 2 || mm < 2 || mm > n) );
        bool ascnd=(xx[n-1] >= xx[0]);
        jl=0;
        ju=n-1;
        while (ju-jl > 1) {
          jm = (ju+jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl=jm;
          else
            ju=jm;
        }
        return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
      }

      // the number of support methods the interpolation method uses
      size_t nb_support_points;
    };



    // Improvements/Problems of hunt:
    // - Hunt gives indices outside of the range of the grid when compiled in Debug or RelWithDebugInfo
    // - make nb_support_points const
    template <>
    struct SearchMethod<SearchMethod_t::Hunt> {

      constexpr static SearchMethod_t Method { SearchMethod_t::Hunt };
      SearchMethod<SearchMethod_t::Hunt>() : correlated{false},
          nb_support_points{2}, last_accessed_index{0}, dj{0} {} 

      void initialize(const Vector_Ref & grid){
        this->dj = std::min(1, 
            static_cast<int>(std::round(std::sqrt(std::sqrt(grid.size())))));
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref & grid) {
        return this->correlated ? this->hunt(x, grid) : this->locate(x, grid);
      }

      int locate(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int jsav{static_cast<int>(this->last_accessed_index)};

        int ju,jm,jl;
        assert( not(n < 2 || mm < 2 || mm > n) );
        bool ascnd=(xx[n-1] >= xx[0]);
        jl=0;
        ju=n-1;
        while (ju-jl > 1) {
          jm = (ju+jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl=jm;
          else
            ju=jm;
        }
        this->correlated = abs(jl-jsav) > this->dj ? 0 : 1;
        jsav = jl;

        this->last_accessed_index = jsav;
        return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
      }

      int hunt(double x, const Vector_Ref & xx){
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int jsav{static_cast<int>(this->last_accessed_index)};

        int jl=jsav, jm, ju, inc=1;
        // hunt size error
        assert( not(n < 2 || mm < 2 || mm > n) );

        bool ascnd=(xx[n-1] >= xx[0]);
        if (jl < 0 || jl > n-1) {
          jl=0;
          ju=n-1;
        } else {
          if ((x >= xx[jl]) == ascnd) {
            for (;;) {
              ju = jl + inc;
              if (ju >= n-1) { ju = n-1; break;}
              else if ((x < xx[ju]) == ascnd) break;
              else {
                jl = ju;
                inc += inc;
              }
            }
          } else {
            ju = jl;
            for (;;) {
              jl = jl - inc;
              if (jl <= 0) { jl = 0; break;}
              else if ((x >= xx[jl]) == ascnd) break;
              else {
                ju = jl;
                inc += inc;
              }
            }
          }
        }
        while (ju-jl > 1) {
          jm = (ju+jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl=jm;
          else
            ju=jm;
        }
        this->correlated = abs(jl-jsav) > this->dj ? 0 : 1;
        jsav = jl;
        this->last_accessed_index = jsav;
        return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
      }

      bool correlated;
      // the number of support methods the interpolation method uses
      size_t nb_support_points;
      size_t last_accessed_index;
      // parameter used to determine if search requests are correlated
      int dj;
    };

    // Untemplated Interpolator class used to upcast the interpolator and store it in memory structures
    class InterpolatorBase {
     public:
      InterpolatorBase() {}; 

      void initialize(double x1, double x2, double error_bound, int max_grid_points = 10000000, int fineness=5) {
        if (x2<x1) { throw std::runtime_error("x2 must be greater x1"); }
        this->x1 = x1;
        this->x2 = x2;
        this->error_bound = error_bound;
        this->max_grid_points = max_grid_points;
        this->fineness = fineness;
      }
     
      // the function which is interpolated
      std::function<double(double)> function{};
      // The boundary points of the range of interpolation
      double x1{0};
      double x2{1};
      double error_bound{1e-5};
      double grid_error{0};

      // If the first derivative of the function at the boundary points is known, one can initialize the interpolator with these values to obtain a higher accuracy with the interpolation method.
      double clamped_boundary_conditions{false}; 
      // y'(x1)
      double yd1{0.};
      // y'(x2)
      double ydn{0.};
      // Describes the finess of the grid, higher value means finer grids. It is up to the GridRational how to use the fineness.
      int fineness{0};
      int max_grid_points{10000000}; //1e7
    };

    // Templated interpolator class used to
    template<class InterpolationMethod, class GridRational, class SearchMethod, class ErrorMethod>
    class InterpolatorTyped : public InterpolatorBase {
     public:
      InterpolatorTyped() : InterpolatorBase(), intp_method{InterpolationMethod()}, grid_rational{GridRational()}, search_method{SearchMethod()}, error_method{ErrorMethod()} {};

      InterpolationMethod intp_method;
      GridRational grid_rational;
      SearchMethod search_method;
      ErrorMethod error_method;
    };

    template<class InterpolationMethod, class GridRational, class SearchMethod, class ErrorMethod=ErrorMethod<ErrorMetric_t::Absolute, ErrorNorm_t::Mean>>
    class Interpolator : public InterpolatorTyped<InterpolationMethod, GridRational, SearchMethod, ErrorMethod> {
      // We specialized the spline methods for uniform grid. The spline methods can be also adapted for general grids, but for uniform grids the computation can much more optimized.
      static_assert (not(InterpolationMethod::Method == InterpolationMethod_t::CubicSpline) || GridRational::GridType == GridType_t::Uniform);
      static_assert (not(SearchMethod::Method == SearchMethod_t::Uniform) || GridRational::GridType == GridType_t::Uniform);
     public: 
      using Parent = InterpolatorTyped<InterpolationMethod, GridRational, SearchMethod, ErrorMethod>;
      Interpolator() : Parent() {}

      void initialize(std::function<double(double)> function, double x1, double x2, double error_bound, int max_grid_points = 10000000, int fineness=5, bool is_benchmarked=false) {
        Parent::initialize(x1, x2, error_bound, max_grid_points, fineness);
        this->function = function;
        this->clamped_boundary_conditions = false;
        this->is_benchmarked = is_benchmarked=false;

        this->initialize();
      }

      void initialize(std::function<double(double)> function, double x1, double x2, double error_bound, double yd1, double ydn, int max_grid_points = 10000000, int fineness=5, bool is_benchmarked=false) {
        Parent::initialize(x1, x2, error_bound, max_grid_points, fineness);
        this->function = function;
        this->clamped_boundary_conditions = true;
        this->yd1 = yd1;
        this->ydn = ydn;
        this->is_benchmarked = is_benchmarked=false;

        this->initialize();
      }

      // Initialization function with an already precomputed grid. For optimization purposes important
      void initialize(std::function<double(double)> function, Vector_t grid) {
        this->function = function;
        this->x1 = grid(0);
        this->x2 = grid(grid.size()-1);
        this->evaluated_grid = this->eval(this->grid);

        this->intp_method.initialize(this->grid, this->evaluated_grid);
        this->search_method.initialize(this->grid);
      }

      double eval(double x) {return this->function(x);}

      // We use evaluate when the function is used and interpolate when the
      // interpolation method is used
      Vector_t eval(const Vector_Ref & grid) {
        Vector_t evaluated_grid = Vector_t::Zero(grid.size());
        for (int i{0}; i<evaluated_grid.size(); i++) {
          evaluated_grid(i) = this->eval(grid(i));
        }
        return evaluated_grid;
      }

      double interpolate(double x) {
        assert (x>=this->x1 && x<=this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate(
            this->grid, this->evaluated_grid,
            x, nearest_grid_index_to_x);
      }

      Vector_t interpolate(const Vector_Ref & points) {
        Vector_t interpolated_points = Vector_t::Zero(points.size());
        for (int i{0}; i<points.size(); i++) {
          interpolated_points(i) = this->interpolate(points(i));
        }
        return interpolated_points;
      }

      double interpolate_derivative(double x) {
        assert (x>=this->x1 && x<=this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate_derivative(
            this->grid, this->evaluated_grid,
            x, nearest_grid_index_to_x);
      }
     
      Vector_t grid{};
      // evaluated by the function to interpolate
      Vector_t evaluated_grid{}; 
      double mean_grid_error{0.};
      double max_grid_error{0.};
      bool is_benchmarked{false};

     private:
      // Should be called after the parameters have been correctly initialized.
      void initialize() {
        this->compute_grid_error();
        while (this->grid_error > this->error_bound && this->grid.size() < this->max_grid_points) {
          this->fineness++;
          this->compute_grid_error();
        }
        if(this->is_benchmarked){
          this->compute_paratemeters_for_evaluation();
        }
      }

      void compute_grid_error() {
        this->grid = 
            this->grid_rational.compute_grid(this->x1,this->x2, this->fineness);
        this->evaluated_grid = this->eval(this->grid);

        if (this->clamped_boundary_conditions) {
          this->intp_method.initialize(this->grid, this->evaluated_grid, this->yd1, this->ydn);
        } else {
          this->intp_method.initialize(this->grid, this->evaluated_grid);
        }
        this->search_method.initialize(this->grid);

        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        Vector_t test_grid_interpolated{this->interpolate(test_grid)};
        Vector_t test_grid_evaluated{this->eval(test_grid)}; 
        Vector_t error_grid{this->error_method.compute_entrywise_error(Vector_Ref(test_grid_interpolated), Vector_Ref(test_grid_evaluated))};
        this->grid_rational.update_errors(Vector_Ref(error_grid));
        this->grid_error = this->error_method.compute_global_error(Vector_Ref(error_grid));
        //if (this->grid.size() % 1==0) {
        //  std::cout << "grid_size=" << this->grid.size() << std::endl;
        //  std::cout << "mean_grid_error=" << this->mean_grid_error << std::endl;
        //  std::cout << "max error=" << this->max_grid_error << std::endl;
        //  std::cout << "error_grid=" << error_grid.head(3) << std::endl;
        //  std::cout << "test_grid_eval=" << test_grid_evaluated.head(3) << std::endl;
        //  std::cout << "test_grid_intp=" << test_grid_interpolated.head(3) << std::endl;
        //}
        //if (grid_rational.grid_size % 50 == 0) {
        //  std::cout << "fineness=" << this->fineness << std::endl;
        //  std::cout << "mean error=" << this->mean_grid_error << std::endl;
        //  std::cout << "max error=" << this->max_grid_error << std::endl;
        //}
      }

      void compute_paratemeters_for_evaluation() {
        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        Vector_t test_grid_interpolated{this->interpolate(test_grid)};
        Vector_t test_grid_evaluated{this->eval(test_grid)}; 
        Vector_t error_grid{this->error_method.compute_entrywise_error(Vector_Ref(test_grid_interpolated), Vector_Ref(test_grid_evaluated))};
        this->max_grid_error = error_grid.maxCoeff();
        this->mean_grid_error = error_grid.mean();
      }

    };

    template<class InterpolationMethod, class GridRational, class SearchMethod, class ErrorMethod=ErrorMethod<ErrorMetric_t::Absolute, ErrorNorm_t::Mean>>
    class InterpolatorVectorized : public InterpolatorTyped<InterpolationMethod, GridRational, SearchMethod, ErrorMethod> {
      static_assert (not(InterpolationMethod::Method == InterpolationMethod_t::CubicSplineVectorized) || GridRational::GridType == GridType_t::Uniform);
      static_assert (not(SearchMethod::Method == SearchMethod_t::Uniform) || GridRational::GridType == GridType_t::Uniform);
     public: 
      using Parent = InterpolatorTyped<InterpolationMethod, GridRational, SearchMethod, ErrorMethod>;
      InterpolatorVectorized() : Parent() {}

      void initialize(std::function<Matrix_t(double)> function, double x1, double x2, double error_bound, int max_grid_points = 10000000, int fineness=5, bool is_benchmarked=false) {
        Parent::initialize(x1, x2, error_bound, max_grid_points, fineness);
        this->function = function;
        Matrix_t result = function(0);
        this->cols = result.cols();
        this->rows = result.rows();
        this->matrix_size = this->cols*this->rows;
        this->is_benchmarked = is_benchmarked;

        this->initialize();
      }

      // Allows the initialization with an already defined grid 
      void initialize(std::function<Matrix_t(double)> function, Vector_t grid) {
        this->function = function;
        Matrix_t result = function(0);
        this->cols = result.cols();
        this->rows = result.rows();
        this->matrix_size = this->cols*this->rows;
        this->x1 = grid(0);
        this->x2 = grid(grid.size()-1);
        this->grid = grid;
        this->evaluated_grid = this->eval(this->grid);

        this->intp_method.initialize(this->grid, this->evaluated_grid);
        this->search_method.initialize(this->grid);
      }


      Vector_t eval(double x) {return Eigen::Map<Vector_t>(this->function(x).data(),this->matrix_size);}

      Matrix_t eval(const Vector_Ref & grid) {
        Matrix_t evaluated_grid = Matrix_t::Zero(grid.size(), this->matrix_size);
        for (int i{0}; i<grid.size(); i++) {
          evaluated_grid.row(i) = this->eval(grid(i));
        }
        return evaluated_grid;
      }


      // Possible source of improvement, if we store the result in the
      // interpolator and return a Eigen::Map, some checks do not have to be made
      Matrix_t interpolate(double x) {
        return Eigen::Map<Matrix_t>(this->interpolate_raw(x).data(), this->rows, this->cols);
      }

      // The vectorized interpolator stores the output of the function as vector,
      //it interpolation results is by default a vector.
      Vector_t interpolate_raw(double x) {
        // x is outside of range
        assert (x>=this->x1 && x<=this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate(
            this->grid, this->evaluated_grid,
            x, nearest_grid_index_to_x);
      }

      Matrix_t interpolate_raw(const Vector_Ref & points) {
        Matrix_t interpolated_points = Matrix_t::Zero(points.size(), this->matrix_size);
        for (int i{0}; i<points.size(); i++) {
          interpolated_points.row(i) = this->interpolate_raw(points(i));
        }
        return interpolated_points;      
      }

      Matrix_t interpolate_derivative(double x) {
        return Eigen::Map<Matrix_t>(this->interpolate_derivative_raw(x).data(), this->rows, this->cols);
      }

      Vector_t interpolate_derivative_raw(double x) {
        // x is outside of range
        assert (x>=this->x1 && x<=this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate_derivative(
            this->grid, this->evaluated_grid,
            x, nearest_grid_index_to_x);
      } 

      // This function exists for benchmark purposes, to compute the memory
      // and function call overhead.
      Matrix_t interpolate_optimal(const double x) const {
        return x * Matrix_t::Ones(this->rows, this->cols);
      }

      
      int rows{0};
      int cols{0};
      int matrix_size{0};
      std::function<Matrix_t(double)> function{};
      Vector_t grid{};
      Matrix_t evaluated_grid{};

      // Sets if additional parameter are calculated for benchmark purposes
      bool is_benchmarked{false};
      // the max error of all entries
      double max_grid_error{0.};
      // the mean error of all entries
      double mean_grid_error{0.};

     private:
      void initialize() {
        this->compute_grid_error();
        while (this->grid_error > this->error_bound && this->grid.size() < this->max_grid_points) {
          this->fineness++;
          this->compute_grid_error();
        }
        if(this->is_benchmarked){
          this->compute_paratemeters_for_evaluation();
        }
      }

      void compute_grid_error() {
        this->grid = 
            this->grid_rational.compute_grid(this->x1,this->x2, this->fineness);
        this->evaluated_grid = this->eval(this->grid);

        this->intp_method.initialize(this->grid, this->evaluated_grid);
        this->search_method.initialize(this->grid);

        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        // (grid_size, row*col)
        Matrix_t test_grid_interpolated{this->interpolate_raw(test_grid)};
        Matrix_t test_grid_evaluated{this->eval(test_grid)};
        this->grid_error = this->error_method.compute_global_error(Matrix_Ref(test_grid_interpolated), Matrix_Ref(test_grid_evaluated));
        // TODO(alex) remove when benchmark finished
        //if (this->grid.size() % 1==0) {
        //  std::cout << "grid_size=" << this->grid.size() << std::endl;
        //  std::cout << "mean_grid_error=" << this->mean_grid_error << std::endl;
        //  std::cout << "max error=" << this->max_grid_error << std::endl;
        //  //std::cout << "max error at=" << maxRow << maxCol << std::endl;
        //  std::cout << "error_mat=" << error_mat.col(0).head(3).transpose() << std::endl;
        //  std::cout << "test_grid_eval=" << test_grid_evaluated.col(0).head(3).transpose() << std::endl;
        //  std::cout << "test_grid_intp=" << test_grid_interpolated.col(0).head(3).transpose() << std::endl;
        //  std::cout << std::endl;
        //}
        //if (grid_rational.grid_size % 50 == 0) {
        //  std::cout << "fineness=" << this->fineness << std::endl;
        //  std::cout << "mean error=" << this->mean_grid_error << std::endl;
        //  std::cout << "max error=" << this->max_grid_error << std::endl;
        //}
      }

      // This function exists for benchmark purposes
      void compute_paratemeters_for_evaluation() {
        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        // (grid_size, row*col)
        Matrix_t test_grid_interpolated{this->interpolate_raw(test_grid)};
        Matrix_t test_grid_evaluated{this->eval(test_grid)};
        Matrix_t error_mat = this->error_method.compute_entrywise_error(Matrix_Ref(test_grid_interpolated), Matrix_Ref(test_grid_evaluated));
        this->max_grid_error = error_mat.maxCoeff();
        this->mean_grid_error = error_mat.mean();
      }

    };

    using InterpolatorVectorized_t = InterpolatorVectorized<
          InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
          GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
          SearchMethod<SearchMethod_t::Uniform>
        >;


  }  // namespace math
}  // namespace rascal
#endif  // SRC_MATH_INTERPOLATOR_HH_

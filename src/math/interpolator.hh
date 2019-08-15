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

    // TODO(alex) make plots of hyp1f1 normalized
    // TODO(alex) look at the graphs again, and make a grid type which is similar
    // to the shape
    // is similar to the graph

    // TODO(alex) currently the grid rational could be static
    template <GridType_t Type, RefinementMethod_t Method>
    struct GridRational {
      constexpr static GridType_t GridType { Type };
      constexpr static RefinementMethod_t RefinementMethod { Method };
    };

    //TODO(alex) adaptive: tree structure blocks{start,end,error,blocks}
    // refine where k largest error exist
    // get grid from (tree -> grid)
    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive> {
      GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>() : grid_meshes{},grid_meshes_error{}, grid_size{0}, max_it{} {}

      //}
      // TODO for this grid rational it makes sense to save the test grid,
      // because it is one step ahead and we can add more points 
      // TODO make hyperparameter k highest error
      // TODO make grid rational init compute, increase_finness() and then remove x1,x2,fineness
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
        //if (this->grid_size % 1000 == 0) {
        //  std::cout << this->grid_size << std::endl;
        //}

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

      // two sortings
      // for refinement key = error
      // for constructing grid x1
      std::forward_list<double> grid_meshes;
      // TODO priority queue does not work well because splitting meshes does.
      // One could make Mesh only owning the x1 and error but the. For the
      // current test grid we need some clever binary tree structure where
      // leafs have errors and their parents add up, Node{x1,error}
      //std::priority_queue<Mesh> grid_meshes_error;
      std::forward_list<double> grid_meshes_error;
      int grid_size;
      std::forward_list<double>::iterator max_it;
    };

    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::Linear> {
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
      // TODO(alex) use grid_size
      int grid_size{0};
    };

    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential> {
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
    };


    // Search, grid can be kept. CubicSpline needs to be vectorized, intp error has to be a bit adapted


    enum class InterpolationMethod_t {CubicSpline, CubicSplineVectorized};

    template <InterpolationMethod_t Type>
    struct InterpolationMethod{};

    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSpline> {
     public:
      InterpolationMethod<InterpolationMethod_t::CubicSpline>() {}

      void initialize(const Vector_Ref & grid,
          const Vector_Ref & evaluated_grid){
        this->compute_second_derivatives_on_grid(grid, evaluated_grid);
        this->h = grid(1) - grid(0); // TODO(alex) asserts uniform grid
        assert (this->h != 0.0); // Bad xa input to routine splint
        this->h_sq_6 = this->h*this->h/6.0;
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
      // TODO(felix) the numerical recipes gives the option to set the first
      // derivative's starting and end point,
      // for now I did not include this option, I do not see now where
      // we would use it
      // TODO(alex) reference numerical recipes
      void compute_second_derivatives_on_grid(
          const Vector_Ref & grid, const Vector_Ref & evaluated_grid) {
        this->second_derivatives = std::move(this->sety2(grid, evaluated_grid));
        //std::cout << "this->second_derivatives" << this->second_derivatives.head(3) << std::endl; 
      }

      // This is done to be close to the numerical recipes implementation in
      // naming while making it more readable.
      // TODO(alex) reference numerical recipes
      inline Vector_t sety2(const Vector_Ref & xv, const Vector_Ref & yv) {
        int n{static_cast<int>(xv.size())};
        Vector_t y2 = Vector_t::Zero(n);
        Vector_t u = Vector_t::Zero(n);
        size_t sig;
        double p;
        y2(0) = 0.0;
        u(0) = 0.0;
        for (int i{1}; i<n-1; i++) {
          sig=(xv(i)-xv(i-1))/(xv(i+1)-xv(i-1));
          p=sig*y2(i-1)+2.0;
          y2(i)=(sig-1.0)/p;
          u(i)=(yv(i+1)-yv(i))/(xv(i+1)-xv(i)) - (yv(i)-yv(i-1))/(xv(i)-xv(i-1));
          u(i)=(6.0*u(i)/(xv(i+1)-xv(i-1))-sig*u(i-1))/p;
        }
        //u(n-1) = 0.0;
        //p=0.0;
        y2(n-1) = 0.0; // (u(n-1)-p*u(n-2))/(p*y2(n-2)+1.0);
        for (int k{n-2};k>0;k--) {
          y2(k)=y2(k)*y2(k+1)+u(k);
        }
        y2(0)= 0.0; // y2(0)*y2(1)+u(0);
        return y2;
      }

      //inline double rawinterp(const Vector_Ref & xx, const Vector_Ref & yy,
      //    size_t j1, double x) {
      //  size_t klo{j1}, khi{j1+1};
      //  const Vector_Ref && y2 = std::move(Vector_Ref(this->second_derivatives));
      //  double h{xx(khi)-xx(klo)};
      //  DEBUG_IF (h == 0.0) { throw std::runtime_error ("Bad xa input to routine splint");}
      //  double a{(xx(khi)-x)/h};
      //  double b{(x-xx(klo))/h};
      //  return a*yy(klo)+b*yy(khi)+((a*a*a-a)*y2(klo) +(b*b*b-b)*y2(khi))*(h*h)/6.0;
      //}

      inline double rawinterp(const Vector_Ref & xx, const Vector_Ref & yy,
          const int & j1, const double & x) {
        int klo{j1}, khi{j1+1};
        // a+b=1
        double a{(xx(khi)-x)/this->h};
        //double b{1-a};
        double b{(x-xx(klo))/this->h};
        return a*yy(klo)+b*yy(khi)+((a*a*a-a)*this->second_derivatives(klo) +(b*b*b-b)*this->second_derivatives(khi))*h_sq_6;
      }

      inline double rawinterp_derivative(const Vector_Ref & xx, const Vector_Ref & yy,
          const int & j1, const double & x) {
        int klo{j1}, khi{j1+1};
        assert (h != 0.0); // Bad xa input to routine splint
        // a+b=1
        double a{(xx(khi)-x)/this->h};
        double b{1-a};
        //double b{(x-xx(klo))/this->h};
        // (yy(khi)-yy(klo))/this->h can be precomputed
        return (yy(khi)-yy(klo))/this->h - ( (3*a*a-1) *this->second_derivatives(klo) + (3*b*b-1)*this->second_derivatives(khi) ) * this->h/6.;
      }

      double h{0};
      double h_sq_6{0}; // h*h/6.0
      Vector_t second_derivatives{};
    };


    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized> {
     public:
      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>() {}

      void initialize(const Vector_Ref & grid,
          const Matrix_Ref & evaluated_grid){
        this->compute_second_derivatives_on_grid(grid, evaluated_grid);
        this->h = grid(1) - grid(0);
        this->h_sq_6 = this->h*this->h/6.0;
      }

      inline Vector_t interpolate(const Vector_Ref & grid,
          const Matrix_Ref & evaluated_grid,
          double x, int nearest_grid_index_to_x) {
        return this->rawinterp(grid, evaluated_grid,
            nearest_grid_index_to_x, x);
      }
     private:
      // TODO(felix) the numerical recipes gives the option to set the first
      // derivative's starting and end point,
      // for now I did not include this option, I do not see now where
      // we would use it
      // TODO(alex) reference numerical recipes
      void compute_second_derivatives_on_grid(
          const Vector_Ref & grid, const Matrix_Ref & evaluated_grid) {
        // TODO(alex) optimize
        this->second_derivatives = std::move(this->sety2(grid, evaluated_grid));
        //std::cout << "this->second_derivatives" << this->second_derivatives.col(0).head(3).transpose() << std::endl; 
      }

      // This is done to be close to the numerical recipes implementation in
      // naming while making it more readable.
      // TODO(alex) reference numerical recipes
      inline Matrix_t sety2(const Vector_Ref & xv, const Matrix_Ref & yv) {
        int n{static_cast<int>(xv.size())};
        Matrix_t y2 = Matrix_t::Zero(n,yv.cols());
        Matrix_t u = Matrix_t::Zero(n,yv.cols());
        size_t sig;
        Vector_t p = Vector_t::Zero(n);
        y2.row(0) = Vector_t::Zero(yv.cols());
        u.row(0) = Vector_t::Zero(yv.cols());
        for (int i{1}; i<n-1; i++) {
          sig=(xv(i)-xv(i-1))/(xv(i+1)-xv(i-1));
          p=sig*y2.row(i-1).array()+2.0;
          y2.row(i)=(sig-1.0)/p.array();
          u.row(i).array() = (yv.row(i+1).array()-yv.row(i).array())/(xv(i+1)-xv(i)) - (yv.row(i).array()-yv.row(i-1).array())/(xv(i)-xv(i-1));
          u.row(i).array() =(6.0*u.row(i).array()/(xv(i+1)-xv(i-1))-sig*u.row(i-1).array())/p.array();
        }
        //u(n-1) = 0.0;
        //p=0.0;
        y2.row(n-1) = Vector_t::Zero(y2.cols()); //(u(n-1)-p*u(n-2))/(p*y2.row(n-2).array()+1.0);
        for (int k{n-2};k>0;k--) {
          y2.row(k).array()=y2.row(k).array()*y2.row(k+1).array()+u.row(k).array();
        }
        y2.row(0) = Vector_t::Zero(y2.cols()); // y2.row(0)*y2.row(1)+u(0);
        return y2;
      }

      inline Vector_t rawinterp(const Vector_Ref & xx, const Matrix_Ref & yy,
          const int & j1, const double & x) {
        int klo{j1}, khi{j1+1};
        assert (h != 0.0); // Bad xa input to routine splint
        // a+b=1
        double a{(xx(khi)-x)/this->h};
        double b{1-a};
        //double b{(x-xx(klo))/this->h};
        // TODO(alex)
        // h_sq_6 * sec_der can be stored to save one multiplication
        return a*yy.row(klo).array()+b*yy.row(khi).array()+((a*a*a-a)*this->second_derivatives.row(klo).array() +(b*b*b-b)*this->second_derivatives.row(khi).array())*h_sq_6;
      }

      inline Vector_t rawinterp_derivative(const Vector_Ref & xx, const Matrix_Ref & yy,
          const int & j1, const double & x) {
        int klo{j1}, khi{j1+1};
        assert (h != 0.0); // Bad xa input to routine splint
        // a+b=1
        double a{(xx(khi)-x)/this->h};
        double b{1-a};
        //double b{(x-xx(klo))/this->h};
        // yy(khi)-yy(klo)/this->h can be precomputed
        // TODO(alex)
        return (yy.row(khi).array()-yy.row(klo).array())/this->h - ( (3*a*a-1) *this->second_derivatives.row(klo).array() + (3*b*b-1)*this->second_derivatives.row(khi).array() ) * this->h/6;
      }

      double h{0};
      double h_sq_6{0}; // h*h/6.0
      Matrix_t second_derivatives{};
    };



    enum class SearchMethod_t {Hunt, Locate, Uniform};

    template <SearchMethod_t Type>
    struct SearchMethod{};

    // TODO(alex) assert that it only can be used with Uniform grids
    template <>
    struct SearchMethod<SearchMethod_t::Uniform> {

      SearchMethod<SearchMethod_t::Uniform>() {} 

      void initialize(const Vector_Ref & grid){
        this->nb_grid_points_per_unit = grid.size()/(grid(grid.size()-1)-grid(0));
        this->x1 = grid(0);
        this->grid_size = grid.size();
        this->search_size = this->grid_size - this->nb_support_points;
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref &) {
        // TODO(alex) make this work for general grids
        // nb_grid_points/unit
        //TODO(alex) save this
        // for heap_based this is less costly
        // (x-grid(0)) * nb_grid_points_per_unit >> 1
        //int raw_index = static_cast<int>(std::floor((x-grid(0)) * nb_grid_points_per_unit)-1);
        // if x1 is zero and x2-x1 = 2**n then we could save some computation time
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


    template <>
    struct SearchMethod<SearchMethod_t::Locate> {

      // TODO(alex) initilize nb_support_points with parameters from the interpolation mtehod
      SearchMethod<SearchMethod_t::Locate>() : 
          nb_support_points{2} {} 

      void initialize(const Vector_Ref & ){
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref & grid) {
        return this->locate(x, grid);
      }

      // TODO(alex) move this to a Base class if we want to implement more
      // search methods
      // TODO(alex) ref numerical recipes
      int locate(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};

        int ju,jm,jl;
        //TODO(alex) activate in debug mode
        //if (n < 2 || mm < 2 || mm > n) throw("locate size error");
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



    template <>
    struct SearchMethod<SearchMethod_t::Hunt> {

      // TODO(alex) initilize nb_support_points with parameters from the interpolation mtehod
      SearchMethod<SearchMethod_t::Hunt>() : correlated{false},
          nb_support_points{2}, last_accessed_index{0}, dj{0} {} 

      void initialize(const Vector_Ref & grid){
        this->dj = std::min(1, 
            static_cast<int>(std::round(std::sqrt(std::sqrt(grid.size())))));
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref & grid) {
        //return this->locate(x, grid);
        return this->correlated ? this->hunt(x, grid) : this->locate(x, grid);
      }

      // TODO(alex) move this to a Base class if we want to implement more
      // search methods
      // TODO(alex) ref numerical recipes
      int locate(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int jsav{static_cast<int>(this->last_accessed_index)};

        // TODO(alex) is this faster than pow(n, 0.25) ?
        int ju,jm,jl;
        //TODO(alex) activate in debug mode
        //if (n < 2 || mm < 2 || mm > n) throw("locate size error");
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
        assert ( (n < 2 || mm < 2 || mm > n) ); // hunt size error

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

    template<class InterpolationMethod, class GridRational, class SearchMethod>
    class Interpolator {
     public: 
      Interpolator() : intp_method{InterpolationMethod()}, grid_rational{GridRational()}, search_method{SearchMethod()} {}

      void initialize(std::function<double(double)> function, double x1, double x2, double precision, int max_grid_points = 10000000) {
        if (x2<x1) {
          throw std::runtime_error("x2 must be greater x1");
        }
        this->function = function;
        this->x1 = x1;
        this->x2 = x2;
        this->precision = precision;
        this->max_grid_points = max_grid_points;

        this->initialize_interpolator();
        this->search_method.initialize(this->grid);
      }

      // TODO(alex) initialize_parameter? precompute  

      // Initialization function given an alread precomputed grid. For optimization purposes
      void initialize(std::function<double(double)> function, double x1, double x2, Vector_t grid) {
        this->function = function;
        this->x1 = x1;
        this->x2 = x2;
        this->grid = grid;
        this->evaluated_grid = this->eval(this->grid);
        this->intp_method.initialize(this->grid, this->evaluated_grid);
        this->search_method.initialize(this->grid);
      }


      void initialize_interpolator() {
        // Fineness starts with zero and is incremently increased
        // this definition is arbitrary but make computation more readable
        this->fineness = 5;
        double error{this->compute_grid_error()};
        // TODO(alex) add some procedure to not get locked if precision is too
        // high
        while (error > this->precision && this->grid.size() < this->max_grid_points) {
          this->fineness++;
          error = this->compute_grid_error();
        }
      }

      // TODO(alex) if I use temporary variables instead of this, does the
      // compiler optimize this?
      double compute_grid_error() {
        this->grid = 
            this->grid_rational.compute_grid(this->x1,this->x2, this->fineness);
        this->evaluated_grid = this->eval(this->grid);

        this->intp_method.initialize(this->grid, this->evaluated_grid);
        this->search_method.initialize(this->grid);

        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        Vector_t test_grid_interpolated{this->interpolate(test_grid)};
        Vector_t test_grid_evaluated{this->eval(test_grid)}; 
        // relative error
        //Vector_t error_grid{2*((test_grid_interpolated - test_grid_evaluated).array()/
        //  (std::numeric_limits< double >::min()+test_grid_interpolated.array().abs() + test_grid_evaluated.array().abs())).abs()};
        // absolute error
        Vector_t error_grid{(test_grid_interpolated - test_grid_evaluated).array().abs()};
        grid_rational.update_errors(Vector_Ref(error_grid));        
        this->max_error = error_grid.maxCoeff();
        this->mean_error = error_grid.mean();
        //if (this->grid.size() % 1==0) {
        //  std::cout << "grid_size=" << this->grid.size() << std::endl;
        //  std::cout << "mean_error=" << this->mean_error << std::endl;
        //  std::cout << "max error=" << this->max_error << std::endl;
        //  std::cout << "error_grid=" << error_grid.head(3) << std::endl;
        //  std::cout << "test_grid_eval=" << test_grid_evaluated.head(3) << std::endl;
        //  std::cout << "test_grid_intp=" << test_grid_interpolated.head(3) << std::endl;
        //}
        //if (grid_rational.grid_size % 50 == 0) {
        //  std::cout << "fineness=" << this->fineness << std::endl;
        //  std::cout << "mean error=" << this->mean_error << std::endl;
        //  std::cout << "max error=" << this->max_error << std::endl;
        //}
        return this->mean_error;
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

      // TODO(alex) test this assumption:
      // this should save one copy operation instead of the above used for 
      // this->evaluated_grid 
      //void eval() {
      //  for (size_t i{0}; i<this->evaluated_grid.size(); i++) {
      //    this->evaluated_grid(i) = this->function(grid(i));
      //  }
      //}

      double interpolate(double x) {
        // TODO(alex) throw runtime error, what is diff?
        assert (x>=this->x1 && x<=this->x2); // x is outside of range, below x1
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate(
            this->grid, this->evaluated_grid,
            x, nearest_grid_index_to_x);
      }

      double interpolate_derivative(double x) {
        // TODO(alex) throw runtime error, what is diff?
        assert (x>=this->x1 && x<=this->x2); // x is outside of range, below x1
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate_derivative(
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
     
      std::function<double(double)> function{};
      double x1{0};
      double x2{1};
      double precision{1e-5};
      double mean_error{0};
      double max_error{0};
      int fineness{0};
      int max_grid_points{10000000}; //1e7
      Vector_t grid{};
      Vector_t evaluated_grid{}; // map matrix function to Vector_t to decrease ressources

      InterpolationMethod intp_method;
      GridRational grid_rational;
      SearchMethod search_method;
    };

    template<class InterpolationMethod, class GridRational, class SearchMethod>
    class InterpolatorVectorized : public Interpolator<InterpolationMethod, GridRational, SearchMethod> {
     public: 
      using Parent = Interpolator<InterpolationMethod, GridRational, SearchMethod>;
      InterpolatorVectorized() : Parent() {}

      void initialize(std::function<Matrix_t(double)> function, double x1, double x2, double precision, int max_grid_points = 10000000) {
        if (x2<x1) {
          throw std::runtime_error("x2 must be greater x1");
        }
        this->function = function;
        Matrix_t result = function(0);
        this->cols = result.cols();
        this->rows = result.rows();
        this->matrix_size = this->cols*this->rows;
        //std::cout << "Matrix shape ("<< this->rows << "," << this->cols << ")";
        //std::cout << std::endl;
        this->x1 = x1;
        this->x2 = x2;
        this->precision = precision;
        this->max_grid_points = max_grid_points;

        this->initialize_interpolator();
        this->search_method.initialize(this->grid);
      }

      void initialize(std::function<Matrix_t(double)> function, double x1, double x2, Vector_t grid) {
        this->function = function;
        Matrix_t result = function(0);
        this->cols = result.cols();
        this->rows = result.rows();
        this->matrix_size = this->cols*this->rows;
        this->x1 = x1;
        this->x2 = x2;
        this->grid = grid;
        this->evaluated_grid = this->eval(this->grid);
        this->intp_method.initialize(this->grid, this->evaluated_grid);
        this->search_method.initialize(this->grid);
      }

      void initialize_interpolator() {
        // Fineness starts with zero and is incremently increased
        // this definition is arbitrary but make computation more readable
        this->fineness = 5;
        double error{this->compute_grid_error()};
        // TODO(alex) add some procedure to not get locked if precision is too
        // high
        while (error > this->precision && this->grid.size() < this->max_grid_points) {
          this->fineness++;
          error = this->compute_grid_error();
        }
      }

      // TODO(alex) if I use temporary variables instead of this, does the
      // compiler optimize this?
      double compute_grid_error() {
        this->grid = 
            this->grid_rational.compute_grid(this->x1,this->x2, this->fineness);
        this->evaluated_grid = this->eval(this->grid);

        this->intp_method.initialize(this->grid, this->evaluated_grid);
        this->search_method.initialize(this->grid);

        // calcu
        this->max_error = 0.;
        this->mean_error = 0.;
        // (grid_size, row*col)
        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        Matrix_t test_grid_interpolated{this->interpolate_raw(test_grid)};
        Matrix_t test_grid_evaluated{this->eval(test_grid)};
        // computes the relative error
        //Matrix_t error_mat = 2*(test_grid_interpolated-test_grid_evaluated).array().abs()/
        //  (std::numeric_limits< double >::min()+test_grid_interpolated.array().abs()+test_grid_evaluated.array().abs());
        
        Matrix_t error_mat = (test_grid_interpolated-test_grid_evaluated).array().abs();
        Matrix_t::Index maxRow, maxCol;
        this->max_error = error_mat.maxCoeff(&maxRow, &maxCol);
        // maximal mean error of each col
        this->mean_error = error_mat.colwise().mean().maxCoeff();
        // print mean error for interest
        //if (this->grid.size() % 1==0) {
        //  std::cout << "grid_size=" << this->grid.size() << std::endl;
        //  std::cout << "mean_error=" << this->mean_error << std::endl;
        //  std::cout << "max error=" << this->max_error << std::endl;
        //  //std::cout << "max error at=" << maxRow << maxCol << std::endl;
        //  std::cout << "error_mat=" << error_mat.col(0).head(3).transpose() << std::endl;
        //  std::cout << "test_grid_eval=" << test_grid_evaluated.col(0).head(3).transpose() << std::endl;
        //  std::cout << "test_grid_intp=" << test_grid_interpolated.col(0).head(3).transpose() << std::endl;
        //  std::cout << std::endl;
        //}
        //if (grid_rational.grid_size % 50 == 0) {
        //  std::cout << "fineness=" << this->fineness << std::endl;
        //  std::cout << "mean error=" << this->mean_error << std::endl;
        //  std::cout << "max error=" << this->max_error << std::endl;
        //}
        return this->mean_error;
      }

      Vector_t eval(double x) {return Eigen::Map<Vector_t>(this->function(x).data(),this->matrix_size);}

      Matrix_t eval(const Vector_Ref & grid) {
        Matrix_t evaluated_grid = Matrix_t::Zero(grid.size(), this->matrix_size);
        for (int i{0}; i<grid.size(); i++) {
          evaluated_grid.row(i) = this->eval(grid(i));
        }
        return evaluated_grid;
      }

      // TODO(alex) 
      // make a function availabe to return as Matrix directly without a map
      // The problem here is that returning a map does not store the underlying
      // matrix, by returning only the map, it is like returning a pointer without
      // storing the referenced structure resulting in segfault.
      // Solution: We need to return a Matrix_t, we cannot
      //
      // TODO(alex) returning the type Matrix_t without makinngEigen::Map<Matrix_t> executes faster, but does result in errors in the tests
      Matrix_t interpolate(double x) {
        return Eigen::Map<Matrix_t>(this->interpolate_raw(x).data(), this->rows, this->cols);
      }

      // This function is only for optimization purposes used.
      // To extract the overhead created by the class
      Matrix_t interpolate_optimal(const double & x) const {
        return x * Matrix_t::Ones(this->rows, this->cols);
      }

      // should be private
      Vector_t interpolate_raw(double x) {
        // TODO(alex) throw runtime error, what is diff?
        assert (x>=this->x1 && x<=this->x2); // x is outside of range
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
     
      int rows{0};
      int cols{0};
      int matrix_size{0};
      std::function<Matrix_t(double)> function{};
      Vector_t grid{};
      Matrix_t evaluated_grid{};
    };

    using InterpolatorVectorized_t = InterpolatorVectorized<
          InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
          GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
          SearchMethod<SearchMethod_t::Uniform>
        >;

  // TODO(alex) make a CRTP calculator and check if this can be merged with the GradientCalutaro stuff

  }  // namespace math
}  // namespace rascal
#endif  // SRC_MATH_INTERPOLATOR_HH_

#ifndef SRC_MATH_INTERPOLATOR_HH_
#define SRC_MATH_INTERPOLATOR_HH_

#include <functional>
#include "math_utils.hh"

namespace rascal {
  namespace math {
    using Vector_Ref = typename Eigen::Ref<const Vector_t>;

    enum class GridType_t {Uniform};
    enum class RefinementMethod_t {HeapBased, Uniform};

    // TODO(alex) currently the grid rational could be static
    template <GridType_t Type, RefinementMethod_t Method>
    struct GridRational {
      constexpr static GridType_t GridType { Type };
      constexpr static RefinementMethod_t RefinementMethod { Method };
    };

    // I think they all can be static
    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::HeapBased> {
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
    };

    enum class InterpolationMethod_t {CubicSpline};

    template <InterpolationMethod_t Type>
    struct InterpolationMethod{};

    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSpline> {
     public:
      InterpolationMethod<InterpolationMethod_t::CubicSpline>(){}

      void initialize(Vector_Ref grid, Vector_Ref evaluated_grid){
        this->compute_second_derivatives_on_grid(grid, evaluated_grid);
      }

      // TODO(felix) the numerical recipes gives the option to set the first
      // derivative's starting and end point,
      // for now I did not include this option, I do not see now where
      // we would use it
      // TODO(alex) reference numerical recipes
      void compute_second_derivatives_on_grid(
          Vector_Ref grid, Vector_Ref evaluated_grid) {
        this->second_derivatives = this->sety2(grid, evaluated_grid);
      }

      double interpolate(Vector_Ref grid, Vector_Ref evaluated_grid,
          double x, size_t nearest_grid_index_to_x) {
        return this->rawinterp(grid, evaluated_grid, nearest_grid_index_to_x, x);
      }

     private:
      // This is done to be close to the numerical recipes implementation in
      // naming while making it more readable.
      // TODO(alex) reference numerical recipes
      Vector_t sety2(Vector_Ref xv, Vector_Ref yv) {
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
        u(n-1) = 0.0;
        p=0.0;
        y2(n-1)=(u(n-1)-p*u(n-2))/(p*y2(n-2)+1.0);
        // beautiful bug when k=0 is included in for loop, 
        // k=-1 will result in a positive number, thus the for loop continues
        // and a segmentation fault in Eigen happens
        // TODO now it is long so we can merge this again
        for (int k{n-2};k>0;k--) {
          y2(k)=y2(k)*y2(k+1)+u(k);
        }
        y2(0)=y2(0)*y2(1)+u(0);
        return y2;
      }

      double rawinterp(Vector_Ref xx, Vector_Ref yy, size_t j1, double x){
        size_t klo{j1}, khi{j1+1};
        Vector_Ref y2 = Vector_Ref(this->second_derivatives);
        double h{xx(khi)-xx(klo)};
        if (h == 0.0) { throw ("Bad xa input to routine splint");}
        double a{(xx(khi)-x)/h};
        double b{(x-xx(klo))/h};
        return a*yy(klo)+b*yy(khi)+((a*a*a-a)*y2(klo)
          +(b*b*b-b)*y2(khi))*(h*h)/6.0;
      }

      Vector_t second_derivatives{};
    };

    enum class SearchMethod_t {Hunt};

    template <SearchMethod_t Type>
    struct SearchMethod{};

    template <>
    struct SearchMethod<SearchMethod_t::Hunt> {

      SearchMethod<SearchMethod_t::Hunt>() : correlated{false},
          nb_support_points{2}, last_accessed_index{0} {} 

      // If the requests to locate seem correlated, then the heuristic is used
      size_t search(double x, Vector_Ref grid) {
        return this->correlated ? this->hunt(x, grid) : this->locate(x, grid);
      } 

      // TODO(alex) move this to a Base class if we want to implement more
      // search methods
      // TODO(alex) ref numerical recipes
      size_t locate(double x, Vector_Ref xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int jsav{static_cast<int>(this->last_accessed_index)};

        // TODO(alex) is this faster than pow(n, 0.25) ?
        int dj = std::min(1, 
            static_cast<int>(std::round(std::sqrt(std::sqrt(n)))));
        int ju,jm,jl;
        if (n < 2 || mm < 2 || mm > n) throw("locate size error");
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
        this->correlated = abs(jl-jsav) > dj ? 0 : 1;
        jsav = jl;

        this->last_accessed_index = jsav;
        return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
      }

      size_t hunt(double x, Vector_Ref xx){
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int dj = std::min(1, 
            static_cast<int>(std::round(std::sqrt(std::sqrt(n)))));
        int jsav{static_cast<int>(this->last_accessed_index)};

        int jl=jsav, jm, ju, inc=1;
        if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
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
        this->correlated = abs(jl-jsav) > dj ? 0 : 1;
        jsav = jl;
        return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
      }

      bool correlated;
      size_t nb_support_points;
      size_t last_accessed_index;
    };

    //template<InterpolationMethod_t InterpolationMethod, class GridRational, SearchMethod_t SearchMethod>
    template<class InterpolationMethod, class GridRational, class SearchMethod>
    class Interpolator {
     public: 
      Interpolator() : intp_method{InterpolationMethod()}, grid_rational{GridRational()}, search_method{SearchMethod()} {}

      void initalize(std::function<double(double)> function, double x1, double x2, double precision) {
        if (x2<x1) {
          throw std::runtime_error("x2 must be greater x1");
        }
        this->function = function;
        this->x1 = x1;
        this->x2 = x2;
        this->precision = precision;

        this->initialize_interpolator();
      }

      void initialize_interpolator() {
        // Fineness starts with zero and is incremently increased
        // this definition is arbitrary but make computation more readable
        this->fineness = 0;
        double error{this->compute_grid_error()};
        // TODO(alex) add some procedure to not get locked if precision is too
        // high
        while (error > this->precision) {
          this->fineness++;
          error = this->compute_grid_error();
        }
      }

      // TODO(alex) if I use temporary variables instead of this, does the
      // compiler optimize this?
      double compute_grid_error() {
        this->grid = 
            this->grid_rational.compute_grid(this->x1,this->x2, this->fineness);
        this->evaluated_grid = this->eval(Vector_Ref(this->grid));

        this->intp_method.initialize(this->grid, this->evaluated_grid);

        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        Vector_t test_grid_interpolated{this->interpolate(Vector_Ref(test_grid))};
        Vector_t test_grid_evaluated{this->eval(Vector_Ref(test_grid))}; 
        return (test_grid_interpolated - test_grid_evaluated).norm();
      }

      double eval(double x) {return this->function(x);}

      // We use evaluate when the function is used and interpolate when the
      // interpolation method is used
      Vector_t eval(Vector_Ref grid) {
        Vector_t evaluated_grid = Vector_t::Zero(grid.size());
        for (int i{0}; i<evaluated_grid.size(); i++) {
          evaluated_grid(i) = this->function(grid(i));
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
        if (x<this->x1) { throw ("x is outside of range, below x1"); }
        if (x>this->x2) { throw ("x is outside of range, above x2"); }
        size_t nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return intp_method.interpolate(
            Vector_Ref(this->grid), Vector_Ref(this->evaluated_grid),
            x, nearest_grid_index_to_x);
      }

      Vector_t interpolate(Vector_Ref points) {
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
      int fineness{0};
      Vector_t grid{};
      Vector_t evaluated_grid{};

      InterpolationMethod intp_method;
      GridRational grid_rational;
      SearchMethod search_method;
    };

  }  // namespace math
}  // namespace rascal
#endif  // SRC_MATH_INTERPOLATOR_HH_

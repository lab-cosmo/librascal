#include "utils/sparsify_utilities.hh"

#define DO_TIMING 1
#ifdef DO_TIMING
  #include <iostream>
  #include <chrono>
  using hrclock = std::chrono::high_resolution_clock;
#endif

namespace rascal {
namespace utils {
  /**
  * Farthest Point Sampling selection of points given the feature matrix
  * 
  * @param feature_matrix is a NxD matrix containing N inputs with D 
  *        features each. Defaults to numpy row-major type
  * @param n_sparse specifies how many points should be selected. The default 
  *        value of zero sorts the whole set of inputs
  * @param i_first_point indicates the index of the first FPS selection. 
  *        Defaults to zero
  * @return a tuple containing the list of n_sparse indices of the selected 
  *        points, and a list of the maximum minimum distance obtained at each
  *        stage
  */
  std::tuple<Eigen::ArrayXi, Eigen::ArrayXd> select_fps(const 
    Eigen::Ref<const RowMatrixXd>& feature_matrix, int n_sparse, 
    int i_first_point) {
    
    int n_inputs = feature_matrix.rows(); // number of inputs 
    
    if (n_sparse == 0) n_sparse = n_inputs;  // defaults to full sorting of the inputs
    if (n_sparse > n_inputs) 
      throw std::runtime_error("Cannot FPS more inputs than those provided");  //TODO <- use the exception mechanism for librascal whatever it is


    /* return arrays */
    auto sparse_indices = Eigen::ArrayXi(n_sparse);     // FPS indices
    auto sparse_minmax_d2 = Eigen::ArrayXd(n_sparse);   // minmax distances (squared)
    
    auto feature_x2 = Eigen::ArrayXd(n_inputs); // square moduli of inputs
    auto list_new_d2 = Eigen::ArrayXd(n_inputs); // list of square distances to latest FPS point
    auto list_min_d2 = Eigen::ArrayXd(n_inputs); // list of minimum square distances to each input
    int i_new{}; 
    double d2max_new{};
    
    for (ssize_t i=0; i<n_inputs; ++i) {  // computes the squared modulus of input points
      feature_x2(i) = feature_matrix.row(i).squaredNorm(); 
    }
    
    // initializes arrays taking the first point provided in input
    sparse_indices(0) = i_first_point;  
    //  distance square to the selected point
    list_new_d2 = feature_x2 + feature_x2(i_first_point) - 
                  2*(feature_matrix*feature_matrix.row(i_first_point).transpose()).array();  
    list_min_d2 = list_new_d2;  // we only have this point....
      
    for (ssize_t i=1 ; i<n_sparse; ++i) {
      d2max_new = list_min_d2.maxCoeff(&i_new);  // picks max dist and its index
      sparse_indices(i) = i_new;
      sparse_minmax_d2(i-1) = d2max_new;    
      // compute distances^2 to the new point
      list_new_d2 = feature_x2 + feature_x2(i_new) - 
          2*(feature_matrix*feature_matrix.row(i_new).transpose()).array();
      
      // this actually returns a list with the element-wise minimum between 
      // list_min_d2(i) and list_new_d2(i)
      list_min_d2 = list_min_d2.min(list_new_d2);
    }
    sparse_minmax_d2(n_sparse-1) = 0;
    
    return std::make_tuple(sparse_indices, sparse_minmax_d2);
  }

  /**
  * Farthest Point Sampling selection of points given the feature matrix. Uses
  * a Voronoi cell algorithm that can be faster when selecting many points, or the
  * dimensionality is relatively low.
  * 
  * @param feature_matrix is a NxD matrix containing N inputs with D 
  *        features each. Defaults to numpy row-major type
  * @param n_sparse specifies how many points should be selected. The default 
  *        value of zero sorts the whole set of inputs
  * @param i_first_point indicates the index of the first FPS selection. 
  *        Defaults to zero
  * @return a tuple containing the list of n_sparse indices of the selected 
  *        points, a list of the maximum minimum distance obtained at each
  *        stage, a list of the assignment of each of the input points to the
  *        FPS selected inputs, and the radius of each Voronoi cell
  */
  std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXi, Eigen::ArrayXd> 
       select_fps_voronoi(const Eigen::Ref<const RowMatrixXd>& feature_matrix, 
                   int n_sparse, int i_first_point) {
    
    int n_inputs = feature_matrix.rows(); // number of inputs 
    int n_features = feature_matrix.cols(); // number of features
  
    // defaults to full sorting of the inputs
    if (n_sparse == 0) n_sparse = n_inputs;  
    if (n_sparse > n_inputs) 
      throw std::runtime_error("Cannot FPS more inputs than those provided");  //TODO <- use the exception mechanism for librascal whatever it is    
    
    // return arrays
    auto sparse_indices = Eigen::ArrayXi(n_sparse);   // FPS indices
    auto sparse_minmax_d2 = Eigen::ArrayXd(n_sparse); // minmax distances^2
    auto voronoi_r2 = Eigen::ArrayXd(n_sparse);       // size^2 of Voronoi cells
    auto voronoi_indices = Eigen::ArrayXi(n_inputs);  // assignment of points to Voronoi cells

    // work arrays
    auto voronoi_i_far = Eigen::ArrayXd(n_sparse);  // index of maximum-d point in each cell
    auto feature_x2 = Eigen::ArrayXd(n_inputs); // square moduli of inputs
    auto list_new_d2 = Eigen::ArrayXd(n_inputs); // list of distances^2 to latest FPS point
    auto list_min_d2 = Eigen::ArrayXd(n_inputs); // list of minimum distances^2 to each input
    auto f_active = Eigen::ArrayXi(n_sparse);  // flags for "active" cells
    auto list_sel_d2q = Eigen::ArrayXd(n_sparse); // list of dist^2/4 to previously selected points
    auto feature_new = Eigen::VectorXd(n_features);  // feaures of the latest FPS point
    auto feature_sel = RowMatrixXd(n_sparse, n_features); // matrix of the features for the active point selection
    
    int i_new{}; 
    double d2max_new{};
    
    for (ssize_t i=0; i<n_inputs; ++i) {  // computes the squared modulus of input points
      feature_x2(i) = feature_matrix.row(i).squaredNorm(); 
    }
    
    // initializes arrays taking the first point provided in input
    sparse_indices(0) = i_first_point;  
    //  distance square to the selected point
    list_new_d2 = feature_x2 + feature_x2(i_first_point) - 
                  2*(feature_matrix*feature_matrix.row(i_first_point).transpose()).array();  
    list_min_d2 = list_new_d2;  // we only have this point....
    
    voronoi_r2 = 0.0;
    voronoi_indices = 0; 
    // picks the initial Voronoi radius and the farthest point index
    voronoi_r2(0) = list_min_d2.maxCoeff( &voronoi_i_far(0) );
    
    feature_sel.row(0) = feature_matrix.row(i_first_point);

#ifdef DO_TIMING    
    // timing code    
    double tmax {0}, tactive{0}, tloop{0};
    long ndist_eval {0}, npoint_skip{0}, ndist_active{0};
    auto gtstart = hrclock::now();
#endif
    for (int i=1 ; i<n_sparse; ++i) {
#ifdef DO_TIMING         
      auto tstart = hrclock::now();
#endif    
      // find the maximum minimum distance and the corresponding point. this is our next FPS
      // the maxmin point must be one of the voronoi radii. So we pick it from this
      // smaller array. Note we only act on the first i items as the array is 
      // filled incrementally
      d2max_new = voronoi_r2.head(i).maxCoeff(&i_new);   // picks max dist and index of the cell
      i_new = voronoi_i_far(i_new);   // the actual index of the fartest point
#ifdef DO_TIMING    
      auto tend = hrclock::now();
      tmax += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();
#endif    
      // store properties of the new FP selection
      sparse_indices(i) = i_new;
      sparse_minmax_d2(i-1) = d2max_new;
      feature_new = feature_matrix.row(i_new);
      // we store indices of the selected features because we can then compute
      // some of the distances with contiguous array operations
      feature_sel.row(i) = feature_new;   
          
      // now we find the "active" Voronoi cells, i.e. those 
      // that might change due to the new selection. 
      f_active = 0;

#ifdef DO_TIMING      
      tstart = hrclock::now();
      ndist_active += i;      
#endif      
      // must compute distance of the new point to all the previous FPS. 
      // some of these might have been computed already, but bookkeeping 
      // could be worse that recomputing (TODO: verify!)
      list_sel_d2q.head(i) = feature_x2(i_new) - 
            2*(feature_sel.topRows(i)*feature_new).array();    
      for (ssize_t j=0; j<i; ++j) list_sel_d2q(j) += feature_x2(sparse_indices(j));
      list_sel_d2q.head(i) *= 0.25;  // triangle inequality reads voronoi_r < d/2 
    
      for (ssize_t j=0; j<i; ++j) {
        // computes distances to previously selected points and uses 
        // triangle inequality to find which voronoi sets might be affected 
        // by the newly selected point divide by four so we don't have to do 
        // that later to speed up later on the bound on distance to the new point
        if (list_sel_d2q(j) < voronoi_r2(j)) {
          f_active(j) = 1;
          voronoi_r2(j) = 0;  // size of active cells will have to be recomputed          
        }
#ifdef DO_TIMING
        else {
          ++npoint_skip;
        }
#endif
    }

#ifdef DO_TIMING    
    tend = hrclock::now();
    tactive += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();    
    
    tstart = hrclock::now();
#endif    
    
    for (ssize_t j=0; j<n_inputs; ++j) {// only considers "active" points      
      int voronoi_idx_j = voronoi_indices(j); 
      if (f_active(voronoi_idx_j) > 0) {  
        // check if we can skip this check for point j. this is a 
        // tighter bound on the distance, since |x_j-x_sel|<rvoronoi_sel
        if (list_sel_d2q(voronoi_idx_j) < list_min_d2(j)) {  
          double d2_j = feature_x2(i_new) + feature_x2(j) 
              - 2*feature_new.dot(feature_matrix.row(j));
#ifdef DO_TIMING
          ndist_eval ++;
#endif
          // we have to reassign point j to the new selection. also,
          // the voronoi center is actually that of the new selection
          if ( d2_j < list_min_d2(j) ) {
            list_min_d2(j) = d2_j;   
            voronoi_indices(j) = voronoi_idx_j = i;            
          }  
        }
        // also must update the voronoi radius
        if ( list_min_d2(j) > voronoi_r2(voronoi_idx_j) ) {
          voronoi_r2(voronoi_idx_j) = list_min_d2(j);
          voronoi_i_far(voronoi_idx_j) = j;   // stores the index of the FP of the cell
        }
      }    
    } 
    
#ifdef DO_TIMING
    tend = hrclock::now();
    tloop += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();    
#endif
       
  }
  sparse_minmax_d2(n_sparse-1) = 0;

#ifdef DO_TIMING  
  auto gtend = hrclock::now();
  
  std::cout<<"Skipped "<<npoint_skip<<" FPS centers of "<< n_sparse*(n_sparse-1)/2<< " - "
    <<npoint_skip*100./(n_sparse*(n_sparse-1)/2)<<"%\n";  
  std::cout<<"Computed "<<ndist_eval<<" distances rather than "<< n_inputs*n_sparse << " - "
    <<ndist_eval*100./(n_inputs*n_sparse)<<" %\n";
    
  std::cout<<"Time total "<<std::chrono::duration_cast<std::chrono::nanoseconds>(gtend-gtstart).count()*1e-9<<"\n";
  std::cout<<"Time looking for max "<<tmax*1e-9<<"\n";
  std::cout<<"Time looking for active "<<tactive*1e-9<<" with "<< ndist_active << " distances\n";
  std::cout<<"Time general loop "<<tloop*1e-9<<"\n";
#endif
  
  return std::make_tuple(sparse_indices, sparse_minmax_d2, 
                         voronoi_indices, voronoi_r2);
}

} // namespace utils
} // namespace rascal

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <iostream>
#include <vector>
#include <chrono>
using hrclock = std::chrono::high_resolution_clock;
//using hrnano = std::chrono::duration_cast<std::chrono::nanoseconds>; 

namespace py = pybind11;
using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

namespace rascal {
/**
* FPS selection of inputs given the feature matrix  
* @param X is a NxD matrix containing N inputs with D features each. Defaults to numpy row-major type
* @param nfps specifies how many points should be selected
* or ...
*/
std::tuple<Eigen::ArrayXi, Eigen::ArrayXd> fps(const Eigen::Ref<const RowMatrixXd>& X, int nfps=0, int ifirst=0) {
    
  auto n = X.rows(); // number of inputs 
  //auto d = X.cols(); // number of features
  
  if (nfps == 0) nfps = n;  // defaults to full sorting of the inputs
  if (nfps > n) 
    throw std::runtime_error("Cannot FPS more inputs than those provided");


  /* return arrays */
  auto ifps = Eigen::ArrayXi(nfps);  // FPS indices
  auto dfps = Eigen::ArrayXd(nfps);  // minmax distances
  
  auto lx2 = Eigen::ArrayXd(n); // square moduli of inputs
  auto ld2 = Eigen::ArrayXd(n); // list of square distances to latest FPS point
  auto ldmin2 = Eigen::ArrayXd(n); // list of minimum square distances to each input
  int imax; double dmax;
  
  for (ssize_t i=0; i<n; ++i) {
    lx2(i) = X.row(i).squaredNorm(); 
  }
  
  // initializes taking first point provided in input
  ifps(0) = ifirst;  
  ld2 = lx2 + lx2(ifirst) - 2*(X*X.row(ifirst).transpose()).array();  //  distance square to the point (minus the norm of the point, which we actually don't need!)
  ldmin2 = ld2;
    
  for (ssize_t i=1 ; i<nfps; ++i) {
    dmax = ldmin2.maxCoeff(&imax);   // picks max dist and its index
    ifps(i) = imax;
    dfps(i-1) = dmax;    
    ld2 = lx2 + lx2(imax) - 2*(X*X.row(imax).transpose()).array(); // compute distances to the new point
    
    ldmin2 = ldmin2.min(ld2);
  }
  dfps(nfps-1) = 0;
  
  return std::make_tuple(ifps,dfps);
}

/**
* FPS selection of inputs given the feature matrix. 
* Computes voronoi assignment of all points to the selected ones, and uses that to accelerate the processing
* @param X is a NxD matrix containing N inputs with D features each
* @param nfps specifies how many points should be selected
* or ...
*/
std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXi, Eigen::ArrayXd> 
  fps_voronoi(const Eigen::Ref<const RowMatrixXd>& X, int nfps=0, int ifirst=0) {
    
  auto n = X.rows(); // number of inputs 
  auto d = X.cols(); // number of features
  
  if (nfps == 0) nfps = n;  // defaults to full sorting of the inputs
  if (nfps > n) 
    throw std::runtime_error("Cannot FPS more inputs than those provided");


  /* return arrays */
  auto ifps = Eigen::ArrayXi(nfps);  // FPS indices
  auto dfps = Eigen::ArrayXd(nfps);  // minmax distances (squared)
  auto rvor = Eigen::ArrayXd(nfps);  // voronoi radii
  auto irvor = Eigen::ArrayXi(nfps); // index associated with the voronoi radius
  auto ivor = Eigen::ArrayXi(n);     // indices of the voronoi assignments
  
  
  auto lx2 = Eigen::ArrayXd(n); // square moduli of inputs
  auto ld2 = Eigen::ArrayXd(n); // list of square distances to latest FPS point
  auto ldmin2 = Eigen::ArrayXd(n); // list of minimum square distances to each input
  auto sdy = Eigen::ArrayXd(nfps); // list of distances to previously selected points
  auto iactive = Eigen::ArrayXi(nfps); // list of "active" FPS points
  auto xmax = Eigen::VectorXd(d);  // temporary matrix with all FPS features to compute fast distances
  
  int imax {0}; double dmax;
  auto Xy = RowMatrixXd(nfps, d); // matrix of the features for the active point selection
  
  for (ssize_t i=0; i<n; ++i) {
    lx2(i) = X.row(i).squaredNorm(); 
  }
  
  // initializes taking first point provided in input
  ifps(0) = ifirst;
  dfps(0) = 0;
  ld2 = lx2 + lx2(ifirst) - 2*(X*X.row(ifirst).transpose()).array();  //  distance square to the point
  ldmin2 = ld2; 
  rvor = 0; ivor = 0;
  rvor(0) = ldmin2.maxCoeff(&irvor(0));
  Xy.row(0) = X.row(ifirst);
  
  double tmax {0}, tactive{0}, tloop{0};
  
  long ndist_eval {0}, npoint_skip{0}, ndist_active{0};
  auto gtstart = hrclock::now();
  for (ssize_t i=1 ; i<nfps; ++i) {
    auto tstart = hrclock::now();
    /*
    dmax = 0;     
    for (ssize_t j=0; j<i; ++j) { // the maxmin point is one of the voronoi radii 
      if (rvor(j) > dmax) {
        dmax = rvor(j); imax = j;
      }
    }*/
    
    // find the maximum minimum distance and the corresponding point. this is our next FPS
    // the maxmin point must be one of the voronoi radii. So we pick it
    dmax = rvor.head(i).maxCoeff(&imax);   // picks max dist and its index
    imax = irvor(imax);                    // we stored the indices corresponding to the farthest point within each cell
    
    auto tend = hrclock::now();
    tmax += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();
    //dmax = ldmin2.maxCoeff(&imax);   // picks max dist and its index
    
    // store properties of the new FP selection
    ifps(i) = imax;
    dfps(i-1) = dmax;
    xmax = X.row(imax);
    Xy.row(i) = xmax;
    
    
    // now we find the "active" Voronoi cells, i.e. those that might change due to
    // the new selection. 
    int nactive = 0;
    iactive = 0;
    tstart = hrclock::now();
    
    // must compute distance of the new point to all the previous FPS. some of these
    // might have been computed already, but bookkeeping would be worse that recomputing (TODO: verify!)
    sdy.head(i) = lx2(imax) - 2*(Xy.topRows(i)*xmax).array();
    for (ssize_t j=0; j<i; ++j) sdy(j) += lx2(ifps(j));
    sdy.head(i) *= 0.25;  // rvor < d/2 is the test we have to do to find the active points    
    ndist_active += i;
    
    for (ssize_t j=0; j<i; ++j) {
        /*sdy(j) = lx2(imax) + lx2(ifps(j)) - 2*xmax.dot(X.row(ifps(j)));
        sdy(j) = dmax + sdy(j) - sqrt(four_dmax*sdy(j));
        if ( rvor(j) > sdy(j) ) {
        */
        //sdy(j) = (lx2(imax) + lx2(ifps(j)) - 2*xmax.dot(X.row(ifps(j))))
        //sdcombo = dmax + sdy(j) - rvor(j);
        //if (four_dmax*sdy(j) > sdcombo*sdcombo) {
        
        // computes distances to previously selected points and uses triangle inequality to find which 
        // voronoi sets might be affected by the newly selected point
        // divide by four so we don't have to do that later to speed up later on the bound on distance to the new point
        if (sdy(j) < rvor(j)) {
          iactive(j) = 1;
          rvor(j) = 0;
          nactive += 1;
        }        
    }
    
    tend = hrclock::now();
    tactive += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();    
    npoint_skip += i-nactive;
    
    tstart = hrclock::now();
    
    /*
      int np = 0;
      for (ssize_t j=0; j<n; ++j) {// only considers "active" points      
      int ivorj = ivor(j); // does nothing 
      if ( iactive(ivorj) == 1) np++;
    }*/
    
    for (ssize_t j=0; j<n; ++j) {// only considers "active" points      
      int ivorj = ivor(j);
      if (iactive(ivorj) == 1) {       
        if (sdy(ivorj) < ldmin2(j)) {  // tighter bound on the distance, since ldmin2(j)<rvor(ivorj)
          double dj2 = lx2(imax) + lx2(j) - 2*xmax.dot(X.row(j));
          ndist_eval ++;
          if ( dj2 < ldmin2(j) ) {
            ivor(j) = ivorj  = i;
            ldmin2(j) = dj2;   
          }  
        }
        // also must update the voronoi radius
        if ( ldmin2(j) > rvor(ivorj) ) {
          rvor(ivorj) = ldmin2(j);
          irvor(ivorj) = j;           
        }
      }    
    } 
    
    tend = hrclock::now();
    tloop += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();    
       
  }
  dfps(nfps-1) = 0;
  auto gtend = hrclock::now();
  
  std::cout<<"Skipped "<<npoint_skip<<" FPS centers of "<< nfps*(nfps-1)/2<< " - "
    <<npoint_skip*100./(nfps*(nfps-1)/2)<<"%\n";  
  std::cout<<"Computed "<<ndist_eval<<" distances rather than "<< n*nfps << " - "
    <<ndist_eval*100./(n*nfps)<<" %\n";
    
  std::cout<<"Time total "<<std::chrono::duration_cast<std::chrono::nanoseconds>(gtend-gtstart).count()*1e-9<<"\n";
  std::cout<<"Time looking for max "<<tmax*1e-9<<"\n";
  std::cout<<"Time looking for active "<<tactive*1e-9<<" with "<< ndist_active << " distances\n";
  std::cout<<"Time general loop "<<tloop*1e-9<<"\n";
  
  return std::make_tuple(ifps,dfps,ivor,rvor);
}

/*
std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXi, Eigen::ArrayXd> 
  fps_list(const Eigen::Ref<const RowMatrixXd>& X, int nfps=0, int ifirst=0) {
    
  auto n = X.rows(); // number of inputs 
  auto d = X.cols(); // number of features
  
  if (nfps == 0) nfps = n;  // defaults to full sorting of the inputs
  if (nfps > n) 
    throw std::runtime_error("Cannot FPS more inputs than those provided");


  
  auto ifps = Eigen::ArrayXi(nfps);  // FPS indices
  auto dfps = Eigen::ArrayXd(nfps);  // minmax distances (squared)
  auto rvor = Eigen::ArrayXd(nfps);  // voronoi radii
  auto irvor = Eigen::ArrayXi(nfps); // index associated with the voronoi radius
  auto ivor = Eigen::ArrayXi(n);     // indices of the voronoi assignments
  auto lvor = std::vector<std::vector<int> >(nfps);
  
  auto lx2 = Eigen::ArrayXd(n); // square moduli of inputs
  auto ld2 = Eigen::ArrayXd(n); // list of square distances to latest FPS point
  auto ldmin2 = Eigen::ArrayXd(n); // list of minimum square distances to each input
  auto sdy = Eigen::ArrayXd(nfps); // list of distances to previously selected points
  auto iactive = Eigen::ArrayXi(nfps); // list of "active" FPS points
  auto xmax = Eigen::VectorXd(d);
  
  int imax {0}; double dmax;
  
  for (ssize_t i=0; i<n; ++i) {
    lx2(i) = X.row(i).squaredNorm(); 
  }
  
  // initializes taking first point provided in input
  ifps(0) = ifirst;
  dfps(0) = 0;
  ld2 = lx2 + lx2(ifirst) - 2*(X*X.row(ifirst).transpose()).array();  //  distance square to the point (minus the norm of the point, which we actually don't need!)
  ldmin2 = ld2; 
  rvor = 0; ivor = 0;
  rvor(0) = ldmin2.maxCoeff(&irvor(0));
  lvor[0].resize(n);
  for (ssize_t i=0; i<n; ++i) lvor[0][i] = i;
  
  double tmax {0}, tactive{0}, tloop{0};
  
  int ndist_eval {0}, npoint_skip{0};
  auto gtstart = hrclock::now();
  for (ssize_t i=1 ; i<nfps; ++i) {
    dmax = 0; 
    auto tstart = hrclock::now();
    
    // the maxmin point must be one of the voronoi radii
    dmax = rvor.head(i).maxCoeff(&imax);   // picks max dist and its index
    imax = irvor(imax);
    auto tend = hrclock::now();
    tmax += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();
    //dmax = ldmin2.maxCoeff(&imax);   // picks max dist and its index
    
    
    ifps(i) = imax;
    dfps(i-1) = dmax;
    xmax = X.row(imax);
    
    double four_dmax = 4*dmax;
    
    int nactive = 0;
    iactive = 0;
    tstart = hrclock::now();
    double sdcombo;    
    for (ssize_t j=0; j<i; ++j) {
        
        // computes distances to previously selected points and uses triangle inequality to find which 
        // voronoi sets might be affected by the newly selected point
        sdy(j) = lx2(imax) + lx2(ifps(j)) - 2*xmax.dot(X.row(ifps(j))); // these had already been computed....
        sdcombo = dmax + sdy(j) - rvor(j);
        if (four_dmax*sdy(j) > sdcombo*sdcombo) {
          iactive(j) = 1;
          nactive += 1;
        }        
    }
    tend = hrclock::now();
    tactive += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();    
    npoint_skip += i-nactive;
    
    tstart = hrclock::now();
    std::vector<int> newvor(0);
    
    double rvori{0}; int irvori{0};
    for (ssize_t k=0; k<i; ++k) { // loops over active points list
     if (not iactive[k]) continue;     
     std::vector<int> updvor(0); updvor.reserve(lvor[k].size());
     double rvork{0}; int irvork{0}; 
     for (int j : lvor[k]) { 
      //std::cout<<j<<" "; 
      //double combo = dmax +sdy(ivorj) -ldmin2(j);
      //if (iactive(ivorj) == 1) 
      { //  and  (four_dmax*sdy(ivorj) > combo*combo) ) {//ldmin2(j) > sdy(ivor(j)) ) {
        double dj2 = lx2(imax) + lx2(j) - 2*xmax.dot(X.row(j));
        ndist_eval ++;
        if ( dj2 < ldmin2(j) ) { // must assign to the new voronoi
          ivor(j) = i;
          newvor.push_back(j);
          ldmin2(j) = dj2;
          if (dj2 > rvori) {
            rvori = dj2;
            irvori = j;
          }
        }
        else
        {
          updvor.push_back(j);
          if ( ldmin2(j) > rvork ) {
            rvork = ldmin2(j);
            irvork = j;           
          }
        }
      }
     } 
     lvor[k] = updvor;     
     rvor(k) = rvork;
     irvor(k) = irvork;
    }
    lvor[i] = newvor;
    rvor(i) = rvori;
    irvor(i) = irvori;
    
    tend = hrclock::now();
    tloop += std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();    
    
  }
  dfps(nfps-1) = 0;
  auto gtend = hrclock::now();
  
  std::cout<<"Skipped "<<npoint_skip<<" FPS centers of "<< nfps*(nfps-1)/2<< " - "
    <<npoint_skip*100./(nfps*(nfps-1)/2)<<"%\n";  
  std::cout<<"Computed "<<ndist_eval<<" distances rather than "<< n*nfps << " - "
    <<ndist_eval*100./(n*nfps)<<" %\n";
  std::cout<<"Time total "<<std::chrono::duration_cast<std::chrono::nanoseconds>(gtend-gtstart).count()*1e-9<<"\n";
  std::cout<<"Time looking for max "<<tmax*1e-9<<"\n";
  std::cout<<"Time looking for active "<<tactive*1e-9<<"\n";
  std::cout<<"Time general loop "<<tloop*1e-9<<"\n";
  
  return std::make_tuple(ifps,dfps,ivor,rvor);
}
*/

}

PYBIND11_MODULE(libfps, m) {
    m.def("fps", &fps, "Selects points from a NxD dimensional feature matrix by farthest point sampling",
          py::arg("X"), py::arg("nfps"), py::arg("ifirst") );
    m.def("fps_voronoi", &fps_voronoi, 
        "Selects points from a NxD dimensional feature matrix by farthest point sampling, using a voronoi-cell algorithm to compute fewer distances",
          py::arg("X"), py::arg("nfps"), py::arg("ifirst") );
    m.def("fps_list", &fps_list, 
        "Selects points from a NxD dimensional feature matrix by farthest point sampling, using a voronoi-cell algorithm to compute fewer distances",
          py::arg("X").noconvert(), py::arg("nfps"), py::arg("ifirst") );
}

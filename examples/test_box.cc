#include <../src/neighbourhood_managers/neighbourhood_manager_base.hh>
#include <../src/neighbourhood_managers/neighbourhood_manager_cell.hh>
//#include <../src/neighbourhood_managers/neighbourhood_manager_lammps.hh>
#include <iostream>
#include <../src/basic_types.h>
#include <Eigen/StdVector>
#include <cmath>
#include <lattice.hh>
using namespace std;
using namespace proteus;

constexpr static int Natom{8};
constexpr static int dim{3};
using ptr_t = double**;
using vVector3d = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >;
using Dim_t = int;
using Ccoord_t = std::array<Dim_t,3>;
using Manager_t = proteus::NeighbourhoodManagerCell;
using Veci = typename Eigen::Array<int, dim, 1>;

const Ccoord_t lin2mult(const Dim_t& index, const Ccoord_t& shape)  {
  Ccoord_t retval{{0}};
  Dim_t factor{1};
  for (Dim_t i{0}; i < dim; ++i) {
    retval[i] = index/factor%shape[i];
    if (i != dim-1 ) {
      factor *= shape[i];
    }
  }
  return retval;
}

const Dim_t mult2lin( const Ccoord_t& coord, const Ccoord_t& shape)  {
    Dim_t index{0};
    Dim_t factor{1};
    for (Dim_t i = 0; i < dim; ++i) {
      index += coord[i]*factor;
      if (i != dim-1 ) {
        factor *= shape[i];
      }
    }
    return index;
}

  
template<class ManagerImplementation>
class Boxes {
  using Manager_t = proteus::NeighbourhoodManagerBase<ManagerImplementation>;
  using Vecd = typename Eigen::Array<double, dim, 1>;
  //using AtomRef_t = typename Manager_t::AtomRef;
  public:
    //! constructor
    Boxes(Manager_t & manager, Lattice& lattice, 
          const std::array<bool,3> & pbc, const double& cutoff_max)
      :sizes{sizes},   pbc{pbc}, manager{manager},lattice{lattice}{
        this->bin_size = cutoff_max;

        Vec3_t reciprocal_lenghts = this->lattice.get_reciprocal_lenghts();
        double face_distance{1};//! distance between cell faces -> height between 2 parallel planes belonging to the parallelepiped
        for (size_t i = 0; i < dim; ++i) {
          if (reciprocal_lenghts[i] > 0){
            face_distance = 1/reciprocal_lenghts[i];
          }
          else {
            face_distance = 1;
          }
          this->nbins[i] = std::max(static_cast<int>(face_distance / this->bin_size), 1 );
          this->nbins_d[i] = this->nbins[i];
          this->neigh_search[i] = static_cast<int>(std::ceil(this->bin_size*this->nbins[i]/face_distance));
          this->nneighbour *= (this->neigh_search[i]*2+1);
          this->nbin *= this->nbins[i];
        }
        //this->innerBoxes.resize(this->nbin);
        
      };
    //! copy constructor
    Boxes(const Boxes & other) = default;
    //! assignment operator
    Boxes & operator=(const Boxes & other) = default;
    virtual ~Boxes() = default;

    /**
     * iterators over `Boxes` dereferences to cell coordinates
     */
    class iterator
    {
    public:
      using value_type = Ccoord_t; //!< stl conformance
      using const_value_type = const value_type; //!< stl conformance
      using difference_type = std::ptrdiff_t; //!< stl conformance
      using iterator_category = std::forward_iterator_tag;//!<stl conformance
      using reference = value_type; //!< stl conformance

      //! constructor
      iterator(const Boxes & boxes, bool begin=true)
      :boxes{boxes}, index{begin? 0: boxes.nbin}{};
      //! dereferencing
      virtual ~iterator() = default;
      
      //! this is the object returned by the iteration
      inline value_type operator*() const{
        return lin2mult(this->index,this->boxes.nbins);
      };
      //! pre-increment
      inline iterator & operator++(){
        ++this->index;
        return *this;
      };
      //! inequality
      inline bool operator!=(const iterator & other) const{
        return (this->index != other.index) || (&this->boxes != &other.boxes);
      };
      //! equality
      inline bool operator==(const iterator & other) const{
        return !(*this!= other);
      };
    


    protected:
      const Boxes& boxes; //!< ref to pixels in cell
      size_t index; //!< index of currect pointed-to pixel
    };
    //! stl conformance
    inline iterator begin() const {return iterator(*this);}
    //! stl conformance
    inline iterator end() const {return iterator(*this, false);}
    //! stl conformance
    inline size_t size() const {return this->nbin;}

    inline void bin_centers(){
      this->binedCenters.resize(this->manager.size());
      Ccoord_t center_bin_coord;
      for (auto center: this->manager ){
        for (int ii{0}; ii<dim;++ii){
          center_bin_coord[ii] = static_cast<int>(std::floor(center.get_position()[ii]/this->nbins[ii]));
        }
        this->binedCenters[mult2lin(center_bin_coord,this->nbins)].push_back(center.get_index());
      }
    }  
    
  protected:
    //NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),count)
    Manager_t & manager;
    Lattice & lattice;
    Ccoord_t nbins{{1,1,1}}; //! number of boxes in each directions.
    Vecd nbins_d; //! to use in bin_centers()
    Dim_t nbin{1}; //! 
    Dim_t nneighbour{1}; //! Max number of neighbour boxes. counts also the center one.
    double bin_size;
    std::array<double,3> sizes; // size of the boxes along the 3 directions
    std::array<bool,3> pbc;
    Ccoord_t neigh_search{{1,1,1}};
    vector<vector<int>> binedCenters;
    //std::vector<BoxInner> innerBoxes;


};

template<class ManagerImplementation>
class Box {
  public:
    //! constructor
    Box(Boxes<ManagerImplementation> & boxes,const Ccoord_t & coord,const std::array<Ccoord_t,2> neigh_bounds)
    :boxes{boxes}
    { 
      
      //! warning: works only with negative a if |a| < b
      /*
      auto branchless_div_mod (int a, int b) {
        return std::array<int,2> d = {(a+b)/b,(a+b)%b};
      }
      */

      for (int dx{neigh_bounds[0][0]}; dx <= neigh_bounds[1][0]; ++dx){
        for (int dy{neigh_bounds[0][1]}; dy <= neigh_bounds[1][1]; ++dy){
          for (int dz{neigh_bounds[0][2]}; dz <= neigh_bounds[1][2]; ++dz){
            Ccoord_t shift{{dx,dy,dz}};
            std::array<Ccoord_t,2> neigh_bin;
            for (int ii{0};ii<dim;++ii){
              /*
              auto aa = branchless_div_mod(coord[ii]+shift[ii],this->boxes.nbins[ii]);
              neigh_bin[0][ii] = aa[0];
              neigh_bin[1][ii] = aa[1];
              */
            }

            for (auto icenter:this->boxes.binedCenters[mult2lin(neigh_bin[1],this->boxes.nbins)]){
              this->neighbour_ids.push_back(icenter);
              this->neighbour_shift.push_back(neigh_bin[0]);
            }
            //this->neighbour_boxes.push_back(neigh_coord);
          }
        }
      }
    };
    //! copy constructor
    Box(const Box & other) = default;
    //! assignment operator
    Box & operator=(const Box & other) = default;
    virtual ~Box() = default;
    /*
    inline const std::vector<Ccoord_t> get_neighbour_boxes(){
      return neighbour_boxes;
    }*/

  protected:
    Boxes<ManagerImplementation> & boxes;
    std::vector<Ccoord_t> neighbour_shift;
    std::vector<Dim_t> neighbour_ids;
};


int main()
{ 
  
  Ccoord_t  nbins{{3,2,2}};
  const std::array<double,3> & sizes{{5,3,3}};
  const std::array<bool,3> & pbc{{true,true,true}};
  double rc{2};
  Manager_t manager;
  

  Eigen::MatrixXd positions_test(3,22); // 3,22
  positions_test << 3.689540159937393,5.123016813620886,1.994119731169116,6.818437242389163,2.630056617829216,6.182500355729062,2.114977334498767,6.697579639059512,1.392155450018263,7.420401523540017,2.432242071439904,6.380314902118375,1.112656394115962,7.699900579442317,3.569715877854675,5.242841095703604,3.122826344932127,5.689730628626151,3.248684682453303,5.563872291104976,2.608353462112637,6.204203511445642,
                  5.035681855581504,2.134827911489532,0.946910011088814,6.223599755982222,4.168634519120968,3.001875247950068,1.980327734683430,5.190182032387606,2.943861424421339,4.226648342649697,5.457161501166098,1.713348265904937,1.501663178733906,5.668846588337130,5.208365510425203,1.962144256645833,2.728127406527150,4.442382360543885,2.839975217222644,4.330534549848392,0.744216089807768,6.426293677263268,
                  4.643695520786083,2.662204050783991,1.250682335857938,6.055217235712136,0.860905287815103,6.444994283754972,4.536108843695142,2.769790727874932,5.609177455068640,1.696722116501434,6.703053268421970,0.602846303148105,3.487609972580834,3.818289598989240,1.436734374347541,5.869165197222533,1.054504320562138,6.251395251007936,3.998423858825871,3.307475712744203,5.323662899811682,1.982236671758393;
  
  Lattice lattice;
  Cell_t cell;
  cell << 6.19,2.41,0.21,
          0.00,6.15,1.02,
          0.00,0.00,7.31;
  lattice.set_cell(cell);

  Boxes<Manager_t> boxes( manager,lattice, pbc,rc);
  cout << boxes.size()  << endl;

  //Eigen::Matrix<double,3,22> positions_sc;//get_scaled2cartesian
  //lattice.get_cartesian2scaled(positions_test,positions_sc);
  /*
  Vec3_t positions_sc;//get_scaled2cartesian
  

  
  Eigen::MatrixXd positions_sc_true(3,22);
  positions_sc_true <<  0.296723741750796,0.703639876645053,0.267449366982580,0.732914251413269,0.164593885207063,0.835769733188786,0.235325758326898,0.765037860068951,0.062053748440047,0.938309869955802,0.075557451185986,0.924806167209864,0.099306843325001,0.901056775070848,0.252988672836491,0.747374945559359,0.336207030045438,0.664156588350411,0.361801281872678,0.638562336523171,0.396587390963031,0.603776227432818,
                        0.713451112366376,0.286724809564553,0.125592877525700,0.874583044405229,0.658294016242431,0.341881905688497,0.219086555224869,0.781089366706060,0.351412276497280,0.648763645433649,0.735260445980754,0.264915475950175,0.165043890527786,0.835132031403143,0.814291210823212,0.185884711107717,0.419672726629966,0.580503195300962,0.371065952685266,0.629109969245663,0.000224293676929,0.999951628253999,
                        0.635252465223814,0.364186600654445,0.171091974809567,0.828347091068692,0.117770901205896,0.881668164672363,0.620534725539691,0.378904340338568,0.767329337218692,0.232109728659567,0.916970351357315,0.082468714520945,0.477101227439239,0.522337838439020,0.196543690061223,0.802895375817036,0.144255037012604,0.855184028865655,0.546980008047315,0.452459057830944,0.728271258524170,0.271167807354089;
  double aa{0};
  for (int jj{0};jj<22;++jj) {
    lattice.get_cartesian2scaled(positions_test.col(jj),positions_sc);
    for (int ii{0};ii<3;++ii) {
      aa = positions_sc_true(ii,jj) - positions_sc[ii];
      
    }
    }
  */
  
  for (auto box : boxes){
    
    for (int i{0};i<dim;++i){
      cout<< box[i]<<" ";
    }
    cout << endl;
    cout << mult2lin(box,nbins)<<endl;
    
  }

  return(0);
}
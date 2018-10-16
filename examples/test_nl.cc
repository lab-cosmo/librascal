#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "basic_types.hh"
#include <iostream>
#include <basic_types.hh>
#include <Eigen/StdVector>
#include <cmath>
using namespace std;

using Manager_t = rascal::StructureManagerCenters;
constexpr static int dim{3};
using Vector_t = Eigen::Matrix<double, dim, 1>;

int main()
{   

    // Manager_t manager{};
    // Eigen::MatrixXd positions(3,8);
    // Eigen::Matrix<int, Eigen::Dynamic, 1> numbers(8);
    // Eigen::MatrixXd cell(3, 3);
    // std::array<int, 3> pbc{{true,true,true}};
    // bool verbose{false};
    // // double cutoff{1.9};
    
    // cell <<
    //     2., 0., 0.,
    //     0., 2., 0.,
    //     0., 0., 2.;

    // positions <<
    //     0.4, 1.4, 0.4, 1.4, 0.4, 1.4, 0.4, 1.4,
    //     0.4, 0.4, 1.4, 1.4, 0.4, 0.4, 1.4, 1.4,
    //     0.4, 0.4, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4;

    // numbers << 1, 2, 3, 4, 5, 6, 7, 8;

    Manager_t manager{};
    Eigen::MatrixXd positions(22, 3);
    Eigen::Matrix<int, Eigen::Dynamic, 1> numbers(22);
    Eigen::MatrixXd cell(3, 3);
    std::array<int, 3> pbc{{true,true,true}};
    bool verbose{false};
    // double cutoff{1.9};
         
      cell <<
        6.19, 2.41, 0.21,
        0.00, 6.15, 1.02,
        0.00, 0.00, 7.31;

      positions <<
        3.689540159937393, 5.123016813620886, 1.994119731169116,
        6.818437242389163, 2.630056617829216, 6.182500355729062,
        2.114977334498767, 6.697579639059512, 1.392155450018263,
        7.420401523540017, 2.432242071439904, 6.380314902118375,
        1.112656394115962, 7.699900579442317, 3.569715877854675,
        5.242841095703604, 3.122826344932127, 5.689730628626151,
        3.248684682453303, 5.563872291104976, 2.608353462112637,
        6.204203511445642, 5.035681855581504, 2.134827911489532,
        0.946910011088814, 6.223599755982222, 4.168634519120968,
        3.001875247950068, 1.980327734683430, 5.190182032387606,
        2.943861424421339, 4.226648342649697, 5.457161501166098,
        1.713348265904937, 1.501663178733906, 5.668846588337130,
        5.208365510425203, 1.962144256645833, 2.728127406527150,
        4.442382360543885, 2.839975217222644, 4.330534549848392,
        0.744216089807768, 6.426293677263268, 4.643695520786083,
        2.662204050783991, 1.250682335857938, 6.055217235712136,
        0.860905287815103, 6.444994283754972, 4.536108843695142,
        2.769790727874932, 5.609177455068640, 1.696722116501434,
        6.703053268421970, 0.602846303148105, 3.487609972580834,
        3.818289598989240, 1.436734374347541, 5.869165197222533,
        1.054504320562138, 6.251395251007936, 3.998423858825871,
        3.307475712744203, 5.323662899811682, 1.982236671758393;

      positions.transposeInPlace();
      numbers << 20, 20, 24, 24, 15, 15, 15, 15, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8;

    manager.update(positions, numbers, cell,
                     Eigen::Map<Eigen::Matrix<int, 3, 1>>{pbc.data()});

    double cut_off{2.};
    
    int mult{10};
    double rc_max{mult*0.5 + cut_off};
    rascal::AdaptorNeighbourList<rascal::StructureManagerCenters> pair_manager{manager, rc_max };
    pair_manager.update();

    for (auto i{0}; i < mult; ++i) {
      auto cutoff_tmp = i*0.5 + cut_off;
      std::vector<std::vector<int>> neigh_ids;
      std::vector<std::vector<double>> neigh_dist;

      std::vector<std::vector<int>> neigh_ids_strict;
      std::vector<std::vector<double>> neigh_dist_strict;
      
      // rascal::AdaptorNeighbourList<rascal::StructureManagerCenters> pair_manager{manager, cutoff_tmp};
      // pair_manager.update();
      if (verbose) {}
      
      std::cout << "Setting up strict manager with rc="<<cutoff_tmp << std::endl;
      rascal::AdaptorStrict<rascal::AdaptorNeighbourList<rascal::StructureManagerCenters>> 
                                        adaptor_strict{pair_manager, cutoff_tmp};
      adaptor_strict.update();


      if (verbose) std::cout << "Setting get adaptor_strict info" << std::endl;
      for (auto center : adaptor_strict) {
        auto icenter{center.get_index()};
        std::vector<int> indices_{};
        std::vector<double> distances_{};

        if (verbose) {
          std::cout << "strict atom out " << center.get_index(); // get_index returns iteration index
          std::cout << " " << center.get_atom_index() << " " ; // get_atom_index returns index from        
          for (int ii{0};ii<3;++ii){
            std::cout << center.get_position()[ii] << " ";
          }
          std::cout << " " << center.get_atom_type() << std::endl;
        }
        int Nneigh{0};
        for (auto neigh : center) {
          double distance{(center.get_position()
                          - neigh.get_position()).norm()};
          Nneigh += 1;        
          indices_.push_back(neigh.get_atom_index());
          distances_.push_back(adaptor_strict.get_distance(neigh));
          
          if (verbose) {
              std::cout << "strict neigh out " << neigh.get_index();
              std::cout << " " << neigh.get_atom_index() << "\t " ;
                
              for (int ii{0};ii<3;++ii){
                std::cout << neigh.get_position()[ii] << ", ";
              }
              std::cout << "\t dist=" << distance;
              std::cout << "\t " << neigh.get_atom_type() << std::endl;
            }
          
        }

        std::cout << "Number of Neighbourg: " << Nneigh << std::endl;
        neigh_ids_strict.push_back(indices_);
        neigh_dist_strict.push_back(distances_);
        if (icenter > 1) break;
      }
    }
    return(0);
}

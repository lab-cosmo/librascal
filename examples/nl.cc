#include <../src/neighbourhood_managers/neighbourhood_manager_cell.hh>
#include <../src/neighbourhood_managers/neighbourhood_manager_lammps.hh>
#include <iostream>
#include <../src/basic_types.h>
#include <Eigen/StdVector>
using namespace std;

using Manager_t = proteus::NeighbourhoodManagerCell;
constexpr static int Natom{8};
constexpr static int dim{3};
using ptr_t = double**;
using vVector3d = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >;
int main()
{

  cout << Natom << " "<< dim << " test " << endl;
  Eigen::MatrixXd pos(3,8);
  pos << 0.00,2.68,1.79,0.00,1.79,0.89,0.89,2.68,
            0.00,2.68,0.00,1.79,1.79,2.68,0.89,0.89,
            1.79,2.68,0.00,0.00,1.79,0.89,2.68,0.89;
  cout << "Now the array pos is:" << endl << pos << endl;
  Eigen::MatrixXd cell(3,3);
  cell << 3.57,0.00,0.00,
            0.00,3.57,0.00,
            0.00,0.00,3.57;
  Eigen::VectorXi num(8);
  num << 6, 6, 6, 6, 6, 6, 6, 6;
  std::array<bool,3> pbc = {1,1,1};

  //std::array<int,3> sss = {-1,1,2};
  Eigen::Vector3d sss;
  sss << -1,1,2;
  Eigen::Vector3d eee;
  eee << 3.2,2,1;
  Eigen::Vector3d ddd;
  ddd = sss.array()*eee.array();

  std::vector<vVector3d> aa;
  


    //int nb_pairs;
  Manager_t manager;
  
  double rc_max{2};
  manager.build(pos,cell,pbc,rc_max);
  //manager.set_positions(pos);
  
  
  for (auto center:manager){
      //center.get_index(); type_name()

      cout << "Center id: " << center.get_atom_index() << endl;
      cout << "Center pos: " << center.get_position().transpose() << endl;
      
      cout << "Neighbour ids: " <<  endl;
      for (auto neigh : center){
          //neigh.get_index();
          cout << "Neigh idx: "<< neigh.get_atom_index() <<  endl ;
          cout << "Neigh pos: "<< neigh.get_position().transpose() <<  endl <<  endl;
          //auto atom = neigh.get_atoms()[1];
          //cout << "Neigh idx :"<< atom.get_index() <<  endl <<  endl;
      }
      cout <<  endl;
  }
  
  /*
  for (auto atom: manager) {
    for (auto pair: atom) {
        cout << " distance between atom " << atom.get_index() << " and atom " << pair.get_index() << "  " << (atom.get_position() - pair.get_position()).norm()  << endl;
    }
  }
  
  cout << manager.get_nb_clusters(1) << " test " << endl;
  */
  //cout <<  nb_pairs << " test " << endl;


    return(0);
}
#include <../src/neighbourhood_managers/neighbourhood_manager_simple.hh>
#include <../src/neighbourhood_managers/neighbourhood_manager_lammps.hh>
#include <iostream>

using namespace std;

using Manager_t = proteus::NeighbourhoodManagerLammps;
constexpr static int Natom{8};
constexpr static int dim{3};
using ptr_t = double**;

int main()
{

  cout << Natomnb << " "<< dim << " test " << endl;

  MatrixXdR pos << 0.00,0.00,1.79;
                  2.68,2.68,2.68;
                  1.79,0.00,0.00;
                  0.00,1.79,0.00;
                  1.79,1.79,1.79;
                  0.89,2.68,0.89;
                  0.89,0.89,2.68;
                  2.68,0.89,0.89;

  MatrixXdR cell << 3.57,0.00,0.00;
                  0.00,3.57,0.00;
                  0.00,0.00,3.57;

  MatrixXdI num << 6, 6, 6, 6, 6, 6, 6, 6;
  MatrixXdB pbc << 1,1,1;
    //int nb_pairs;
  Manager_t manager;

  manager.set_positions(pos);
  
  /*
  for (auto atom: manager) {
    for (auto pair: atom) {
        cout << " distance between atom " << atom.get_index() << " and atom " << pair.get_index() << "  " << (atom.get_x() - pair.get_x()).norm()  << endl;
    }
  }
  
  cout << manager.get_nb_clusters(1) << " test " << endl;
  */
  //cout <<  nb_pairs << " test " << endl;


    return(0);
}
#include <neighbourhood_managers/neighbourhood_manager_cell.hh>
#include <neighbourhood_managers/neighbourhood_manager_lammps.hh>
#include <neighbourhood_managers/field.hh>
#include <iostream>
#include <basic_types.h>
#include <Eigen/StdVector>
#include <cmath>
using namespace std;

using Manager_t = rascal::NeighbourhoodManagerCell;
constexpr static int Natom{8};
constexpr static int dim{3};
using Vector_t = Eigen::Matrix<double, dim, 1>;

int main()
{   
    Manager_t manager;
    Eigen::MatrixXd cell(3,3);
    cell << 6.19,2.41,0.21,
            0.00,6.15,1.02,
            0.00,0.00,7.31;
    Eigen::MatrixXd positions(3,22); // 3,22
    positions << 3.689540159937393,5.123016813620886,1.994119731169116,6.818437242389163,2.630056617829216,6.182500355729062,2.114977334498767,6.697579639059512,1.392155450018263,7.420401523540017,2.432242071439904,6.380314902118375,1.112656394115962,7.699900579442317,3.569715877854675,5.242841095703604,3.122826344932127,5.689730628626151,3.248684682453303,5.563872291104976,2.608353462112637,6.204203511445642,
                5.035681855581504,2.134827911489532,0.946910011088814,6.223599755982222,4.168634519120968,3.001875247950068,1.980327734683430,5.190182032387606,2.943861424421339,4.226648342649697,5.457161501166098,1.713348265904937,1.501663178733906,5.668846588337130,5.208365510425203,1.962144256645833,2.728127406527150,4.442382360543885,2.839975217222644,4.330534549848392,0.744216089807768,6.426293677263268,
                4.643695520786083,2.662204050783991,1.250682335857938,6.055217235712136,0.860905287815103,6.444994283754972,4.536108843695142,2.769790727874932,5.609177455068640,1.696722116501434,6.703053268421970,0.602846303148105,3.487609972580834,3.818289598989240,1.436734374347541,5.869165197222533,1.054504320562138,6.251395251007936,3.998423858825871,3.307475712744203,5.323662899811682,1.982236671758393;
    
    rascal::VecXi numbers(22);
    numbers << 20, 20, 24, 24, 15, 15, 15, 15,  8,  8,  8,  8,  8,  8,  8,  8,  8, 8,  8,  8,  8,  8;
    std::array<bool,3> pbc{{true,true,true}};
    double cutoff_max{3}; 
    rascal::VecXi center_ids(22);
    center_ids << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21;
    manager.build(positions,numbers,center_ids,cell,pbc,cutoff_max);
    
    Eigen::MatrixXd positions_test(3,22); // 3,22
    positions_test << 3.689540159937393,5.123016813620886,1.994119731169116,6.818437242389163,2.630056617829216,6.182500355729062,2.114977334498767,6.697579639059512,1.392155450018263,7.420401523540017,2.432242071439904,6.380314902118375,1.112656394115962,7.699900579442317,3.569715877854675,5.242841095703604,3.122826344932127,5.689730628626151,3.248684682453303,5.563872291104976,2.608353462112637,6.204203511445642,
                5.035681855581504,2.134827911489532,0.946910011088814,6.223599755982222,4.168634519120968,3.001875247950068,1.980327734683430,5.190182032387606,2.943861424421339,4.226648342649697,5.457161501166098,1.713348265904937,1.501663178733906,5.668846588337130,5.208365510425203,1.962144256645833,2.728127406527150,4.442382360543885,2.839975217222644,4.330534549848392,0.744216089807768,6.426293677263268,
                4.643695520786083,2.662204050783991,1.250682335857938,6.055217235712136,0.860905287815103,6.444994283754972,4.536108843695142,2.769790727874932,5.609177455068640,1.696722116501434,6.703053268421970,0.602846303148105,3.487609972580834,3.818289598989240,1.436734374347541,5.869165197222533,1.054504320562138,6.251395251007936,3.998423858825871,3.307475712744203,5.323662899811682,1.982236671758393;
        
    rascal::VecXi neighlist;
    //double rc2{cutoff_max*cutoff_max};
    double rc{cutoff_max};
    int atom_counter{0};
    //int pair_counter{0};

    /*
    rascal::Vec3i_t nbins_c,neigh_search,coord;
    neigh_search << 1,1,1;
    std::array<std::array<int, 3>,2> neigh_bounds{{ {-neigh_search[0],-neigh_search[1],-neigh_search[2]},
                                                        { neigh_search[0], neigh_search[1], neigh_search[2]} }};
    coord << 0,0,0;
    nbins_c << 1,2,2;
    rascal::Vec3i_t shift,neighbour_bin_idx_c;
    Vector_t neighbour_bin_shift;
    std::array<int,2> div_mod;
    int bin_id{0};
    for (int dx{neigh_bounds[0][0]}; dx <= neigh_bounds[1][0]; ++dx){
      for (int dy{neigh_bounds[0][1]}; dy <= neigh_bounds[1][1]; ++dy){
        for (int dz{neigh_bounds[0][2]}; dz <= neigh_bounds[1][2]; ++dz){
          shift << coord(0)+dx,coord(1)+dy,coord(2)+dz;
          
          for (int ii{0};ii<3;++ii){
            //branchless_div_mod(coord(ii)+shift(ii),nbins_c(ii),div_mod);
            rascal::internal::div_mod(shift(ii),nbins_c(ii),div_mod);
            neighbour_bin_shift[ii] = static_cast<double>(div_mod[0]);
            neighbour_bin_idx_c[ii] = div_mod[1];
          }

          bin_id = rascal::internal::mult2lin(neighbour_bin_idx_c,nbins_c);
          int trois{3};
          trois = bin_id*trois;
        }
      }
    std::vector<std::vector<AtomRef_t>> neighbour_bin_id;
    std::vector<size_t> number_of_neighbours_stride;
    std::vector<std::vector<AtomRef_t>> neighbour_atom_index;
    this->neighbour_atom_index[box_id].size()
    get_box_nb
    }*/
    Vector_t shift;
    /*
    Manager_t::Box box{manager.get_box(0)};
    for (size_t ii{0}; ii < box.get_number_of_neighbour_box(); ++ii){
        cout<<box.get_neighbour_bin_index(ii)<<": " << ii << ": ";
        shift = box.get_neighbour_bin_shift(ii);
        cout<< shift.transpose() <<endl;
    }
    cout << endl;

    for (size_t bin_id{0}; bin_id < manager.get_box_nb(); ++bin_id){
        if (bin_id != 0){continue;}
        Manager_t::Box box{manager.get_box(bin_id)};
        for(size_t neigh_part_id{0}; neigh_part_id < manager.neighbour_atom_index[bin_id].size(); ++neigh_part_id ){
            int part_index{manager.neighbour_atom_index[bin_id][neigh_part_id].get_index()};
            int shift_index{manager.neighbour_bin_id[bin_id][neigh_part_id].get_index()};
            shift = box.get_neighbour_bin_shift(shift_index);
            cout << bin_id << ": " << part_index << ": " << neigh_part_id <<  ": " << shift.transpose()  << endl;
        }
    }
    cout << endl;
    */
    for (auto center: manager) {
        if (atom_counter - center.get_atom_index() != 0 ){
            cout << "index " << center.get_atom_index() << endl;
        }
        ++atom_counter;
        
        
        for (int ii{3};ii<3;++ii){
            if (positions_test(ii,center.get_atom_index()) - center.get_position()[ii] != 0 ){
                cout << "position " << center.get_atom_index() << endl;
            }
        }


        if (center.get_atom_index() != 1){continue;}
        cout << "Center atom index: " << center.get_atom_index() << endl;
        cout << "Neighbour indices: "<< endl;
        for (auto neigh : center){
            //int part_index{neigh.get_atom_index()};
            
            //Vector_t r,rcent,shift;
            Vector_t r,rcent;
            shift = neigh.get_atom_shift();
            
            
            r = neigh.get_position() +  cell * shift;
            rcent = r -  center.get_position();
            double d{rcent.norm()};
            /*
            cout  << shift.transpose()  << endl;
            //cout << part_index << ": " << neigh.get_index() <<  ": " << shift.transpose()  << endl;
            cout  << "Neigh: " << part_index << endl;
            cout  << "Shift: " << shift.transpose() << endl;
            cout  << "pos: " << neigh.get_position().transpose() << endl;
            cout  << "pos shifted: " << r.transpose() << endl;
            cout  << "center pos: " << center.get_position().transpose() << endl;
            cout  << "distance: " << d << endl;
            cout << endl;
            */
            //double d2{(neigh.get_position() +  cell.transpose() * neigh.get_atom_shift() - center.get_position()).squaredNorm()};
            
            if (d < rc ){
                cout  << neigh.get_atom_index() << ":"<<neigh.get_atom_shift().transpose() << "; ";
            }
            
            
        }

        cout << endl;
        //break;
    
        
    }

    return(0);
}
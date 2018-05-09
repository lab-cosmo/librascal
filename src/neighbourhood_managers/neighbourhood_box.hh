/**
 * file   neighbourhood_box.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager for lammps neighbourhood lists
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef NEIGHBOURHOOD_BOX_H
#define NEIGHBOURHOOD_BOX_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/field.hh"
#include "lattice.hh"
#include "basic_types.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>


namespace rascal {

  template <class ManagerImplementation>
  class Box {
  public:
    using AtomRef_t = typename NeighbourhoodManagerBase<ManagerImplementation>::AtomRef;
    using Vector_ref = typename NeighbourhoodManagerBase<ManagerImplementation>::Vector_ref;
    using Vector_t = typename NeighbourhoodManagerBase<ManagerImplementation>::Vector_t;
    //! Default constructor
    Box() = default;

    //! constructor
    Box(NeighbourhoodManagerBase<ManagerImplementation> & manager, const Vec3i_t& coord, 
                const std::array<std::array<Dim_t, 3>,2>& neigh_bounds, const Vec3i_t& nbins_c)
                :manager{manager},centers{},neighbour_bin_shift{},neighbour_bin_index{},neighbour_ids{}
    { 
    
      //std::function<void (int,int,std::array<int,2>)> branchless_div_mod = [](int a, int b,std::array<int,2> d) {d = {a/b-1,(a+b)%b};};
      
      this->coordinates = coord;
      Vec3i_t shift,neighbour_bin_idx_c;
      Vector_t neighbour_bin_shift;
      std::array<int,2> div_mod;
      int bin_id{0};
      for (int dx{neigh_bounds[0][0]}; dx <= neigh_bounds[1][0]; ++dx){
        for (int dy{neigh_bounds[0][1]}; dy <= neigh_bounds[1][1]; ++dy){
          for (int dz{neigh_bounds[0][2]}; dz <= neigh_bounds[1][2]; ++dz){
            shift << dx,dy,dz;
            
            for (int ii{0};ii<3;++ii){
              //branchless_div_mod(coord(ii)+shift(ii),nbins_c(ii),div_mod);
              internal::div_mod(coord(ii)+shift(ii),nbins_c(ii),div_mod);
              neighbour_bin_shift[ii] = static_cast<double>(div_mod[0]);
              neighbour_bin_idx_c[ii] = div_mod[1];
            }

            bin_id = internal::mult2lin(neighbour_bin_idx_c,nbins_c);
            this->neighbour_bin_index.push_back(bin_id);
            this->neighbour_bin_shift.push_back(neighbour_bin_shift);
            
          }
        }
      }
    };
    //! copy constructor
    Box(const Box & other) = default;
    //! assignment operator
    Box & operator=(const Box & other) = default;


    virtual ~Box() = default;

    inline void push_center_back(const int& id){
      this->centers.push_back(AtomRef_t(this->manager,id));
    }
    inline void push_neighbour_back(const int& id){
      this->neighbour_ids.push_back(AtomRef_t(this->manager,id));
    }
    /*
    inline size_t get_number_of_centers(){
      return this->centers.size();
    }

    inline size_t get_number_of_neighbour(){
      size_t size{this->neighbour_ids.size()};
      return size;
    }

    inline std::vector<AtomRef_t> get_neighbour_ids(){
      std::vector<AtomRef_t> out(this->neighbour_ids);
      return out;
    }

    inline  int get_neighbour_index(int j_index){
      return this->neighbour_ids[j_index].get_index();
    }

    inline  std::vector<AtomRef_t> get_centers(){
      std::vector<AtomRef_t> out{this->centers};
      return out;
    }
    inline  std::vector<Dim_t> get_neighbour_bin_ids(){
      std::vector<Dim_t> out(this->neighbour_bin_index);
      return out;
    }*/

    void set_number_of_neighbour(const size_t neighbour_number,const std::vector<size_t>& neighbour_bin_numbers){
      this->neighbour_number = neighbour_number;
      int box_id{0};
      for (auto nb : neighbour_bin_numbers){
        this->neighbour_bin_numbers.push_back(nb);
        for (auto center_id{0};center_id<nb;++center_id){
          std::array<int,2> p1{{box_id,center_id}};
          this->neighbour_box_atom_ids.push_back(p1);
        }
        ++box_id;
      }
      
    }

    inline Vector_ref get_neighbour_bin_shift(const int& box_id){
      auto * xval{this->neighbour_box_atom_ids[box_id].data()};
      return Vector_ref(xval);
    }

    std::vector<AtomRef_t> centers;
    std::vector<Vector_t,Eigen::aligned_allocator<Vector_t>> neighbour_bin_shift; //stores double for the dot product with the Cell vectors
    std::vector<int> neighbour_bin_index;
    size_t neighbour_number;
    std::vector<size_t> neighbour_bin_numbers;
    std::vector<std::array<int,2> > neighbour_box_atom_ids; 
    // TODO replace neighbour_ids by an iterator that iterates over the centers in the neighbouring boxes. this should be done in the manager
    std::vector<AtomRef_t> neighbour_ids; 
    Vec3i_t coordinates;
    NeighbourhoodManagerBase<ManagerImplementation> & manager;

    template <class ManagerImplementation>
    class Neighbour_Box;
    /**
     * iterators over the centers of the neighbour boxes 
     */
    template <class ManagerImplementation>
    class iterator;
    using Iterator_t = iterator<ManagerImplementation>;

    //! stl conformance
    inline Iterator_t begin() const {return Iterator_t(*this);}
    //! stl conformance
    inline Iterator_t end() const {return Iterator_t(*this, false);}
    //! stl conformance
    inline size_t size() const {return this->neighbour_number;}

  protected:
    
  };

  template <class ManagerImplementation>
  class Box<ManagerImplementation>::Neighbour_Box {
      public:
         //! Default constructor
        Neighbour_Box() = delete;

        //! constructor from iterator
        //AtomRef(Manager_t & manager, int id): manager{manager}, index{id}{}
        Neighbour_Box(Box<ManagerImplementation> & box, int id): box{box}, index{id} {}
        //! Copy constructor
        Neighbour_Box(const Neighbour_Box &other) = default;

        //! Move constructor
        Neighbour_Box(Neighbour_Box &&other) = default;

        //! Destructor
        virtual ~Neighbour_Box() = default;

        //! Copy assignment operator
        Neighbour_Box& operator=(const Neighbour_Box &other) = delete;

        //! Move assignment operator
        Neighbour_Box& operator=(Neighbour_Box &&other) = default;

        //! return index
        inline int get_index() const {return this->index;}

        inline int get_neighbour_bin_index() {
          return this->box.neighbour_box_atom_ids[this->index][0];
        }
        inline int get_neighbour_atom_index() {
          return this->box.neighbour_box_atom_ids[this->index][1];
        }
        //! return position vector
        inline Vector_ref get_position_shift() {
          auto && bin_id{this->get_neighbour_bin_index()};
          return this->box.get_neighbour_bin_shift(bin_id);
        }
        //! return position vector


      protected:
        int index;
        Box<ManagerImplementation>& box;
    };
  
  template <class ManagerImplementation>
  class Box<ManagerImplementation>::iterator
    {
    public:
      using value_type = neighbour_box; //!< stl conformance
      using const_value_type = const value_type; //!< stl conformance
      using difference_type = std::ptrdiff_t; //!< stl conformance
      using iterator_category = std::forward_iterator_tag;//!<stl conformance
      using reference = value_type; //!< stl conformance

      //! constructor
      iterator(const Box& box, bool begin=true)
      :box{box}, index{begin? 0: box.neighbour_number}{};
      //! dereferencing
      virtual ~iterator() = default;
      
      //! this is the object returned by the iteration
      inline value_type operator*() const{
        return Neighbour_Box<>;
      };
      //! pre-increment
      inline iterator & operator++(){
        ++this->index;
        return *this;
      };
      //! inequality
      inline bool operator!=(const iterator & other) const{
        return (this->index != other.index) || (&this->box != &other.box);
      };
      //! equality
      inline bool operator==(const iterator & other) const{
        return !(*this!= other);
      };
    
    protected:
      const Box<ManagerImplementation>& box; //!< ref to pixels in cell
      size_t index; //!< index of currect pointed-to pixel
    };
  //----------------------------------------------------------------------------//
}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_CELL_H */

How to write a Representation Manager
-------------------------------



###Write a RepresentationManager


A representation manager is an object that builds a representation of the atomic structure contained in *a* structure manager. This class:

- inherits publicly from <rascal::RepresentationManagerBase> to follow its interface and use some of the common utilies shared by such class.

- is templated by a structure manager to be able to build the representation efficiently.

- uses the <rascal::Property> class to store the representation's features.

The behaviour of a representation manager is defined at construction by a structure manager and a structure containing a set of parameters. Within these parameters some
  - define behaviours or options for the class, e.g. choosing a particular implementation from a set of conceptually equivalent methods.
  - set the hyperparameters of the representation.

Note that there is one representation manager per structure manager.

To illustrate the basic structure that a new representation that would be implemented in ``representation_manager_custom.hh`` should follow, let's take the example of the  <rascal::RepresentationManagerSortedCoulomb>. 


The representation starts with the definition of some useful short hands

~~~~~~~~~~~~~{.hh}
/**
 * Implementation of the Environmental Coulomb Matrix
 */
template <class StructureManager>
class RepresentationManagerSortedCoulomb : public RepresentationManagerBase {
 public:
  using Manager_t = StructureManager;
  using ManagerPtr_t = std::shared_ptr<Manager_t>;
  using Parent = RepresentationManagerBase;
  // type of the hyperparameters
  using Hypers_t = typename Parent::Hypers_t;
  // numeric type for the representation features
  using Precision_t = typename Parent::Precision_t;
  // type of the data structure for the representation feaures
  using Property_t =
      Property<Precision_t, 1, 1, Manager_t, Eigen::Dynamic, 1>;
  template <size_t Order>
  // short hand type to help the iteration over the structure manager
  using ClusterRef_t = typename Manager_t::template ClusterRef<Order>;
  // type of the datastructure used to register the list of valid
  // hyperparameters
  using ReferenceHypers_t = Parent::ReferenceHypers_t;
~~~~~~~~~~~~~
followed by the definition of its constructors and destructor
~~~~~~~~~~~~~{.hh}
 //! Constructor
  RepresentationManagerSortedCoulomb(ManagerPtr_t sm, const Hypers_t & hyper)
      : structure_manager{std::move(sm)}, central_decay{},
        interaction_cutoff{}, interaction_decay{}, coulomb_matrices{*sm} {
    this->check_hyperparameters(this->reference_hypers, hyper);
    // Extract the options and hyperparameters
    this->set_hyperparameters(hyper);
    // additional checks specific to the coulomb matrix representation
    this->check_size_compatibility();
  }

  //! Copy constructor
  RepresentationManagerSortedCoulomb(
      const RepresentationManagerSortedCoulomb & other) = delete;

  //! Move constructor
  RepresentationManagerSortedCoulomb(
      RepresentationManagerSortedCoulomb && other) = default;

  //! Destructor
  virtual ~RepresentationManagerSortedCoulomb() = default;

  //! Copy assignment operator
  RepresentationManagerSortedCoulomb &
  operator=(const RepresentationManagerSortedCoulomb & other) = delete;

  //! Move assignment operator
  RepresentationManagerSortedCoulomb &
  operator=(RepresentationManagerSortedCoulomb && other) = default;
~~~~~~~~~~~~~
the declaration of the concrete implementation of the RepresentationManager interface
~~~~~~~~~~~~~{.hh}
//! compute representation
  void compute();

  //! set hypers
  void set_hyperparameters(const Hypers_t &);

  //! getter for the representation
  Eigen::Map<const Eigen::MatrixXd> get_representation_full() {
    auto nb_centers{this->structure_manager->size()};
    auto nb_features{this->get_n_feature()};
    auto & raw_data{this->coulomb_matrices.get_raw_data()};
    Eigen::Map<const Eigen::MatrixXd> representation(raw_data.data(),
                                                     nb_features, nb_centers);
    return representation;
  }

  //! get the raw data of the representation
  std::vector<Precision_t> & get_representation_raw_data() {
    return this->coulomb_matrices.get_raw_data();
  }

  Data_t & get_representation_sparse_raw_data() { return this->dummy; }

  //! get the size of a feature vector
  size_t get_feature_size() { return this->coulomb_matrices.get_nb_comp(); }

  //! get the number of centers for the representation
  size_t get_center_size() { return this->coulomb_matrices.get_nb_item(); }
~~~~~~~~~~~~~
and the declaration of some functions for internal use in the protected section.

The end of the class contains the different internal variables needed by the class
~~~~~~~~~~~~~{.hh}
  // Reference to the structure manager
  ManagerPtr_t structure_manager;
  // list of hyperparameters specific to the coulomb matrix
  // spherical cutoff for the atomic environment
  double central_cutoff{};

  double central_decay{};
  double interaction_cutoff{};
  double interaction_decay{};
  // at least equal to the largest number of neighours
  size_t size{};

  Data_t dummy{};

  Property_t coulomb_matrices;

  //! reference the requiered hypers
  ReferenceHypers_t reference_hypers{
      {"central_decay", {}},
      {"interaction_cutoff", {}},
      {"interaction_decay", {}},
      {"size", {}},
      {"sorting_algorithm", {"distance", "row_norm"}},
  };
~~~~~~~~~~~~~

##Write the python bindings of the new representation


We use ``pybind11`` to handle the generation of python bindings. The first step to register the new representation is to include the file in ``bindings/bind_include.hh``. Then, you need to explicitly register your representation manager for every possible structure manager stack that you want to make available to the python side in the <add_representation_managers>. Here is an example on how it is done for the sorted coulomb representation
~~~~~~~~~~~~~{.py}
// Defines a particular structure manager type
using Manager_t =
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>;
// Defines the representation manager type for the particular structure
// manager
using Representation1_t = RepresentationManagerSortedCoulomb<Manager_t>;
// Bind the interface of this representation manager
auto rep_sorted_coulomb =
    add_representation_manager<Representation1_t>(mod, m_throwaway);
~~~~~~~~~~~~~
The last step is to write a python class in ``bindings/rascal/representations/`` to simplify the use of the representation from the python side. You can use :class:`SortedCoulombMatrix` as a template



Efficient Implementation of Options


The switch between several implementations of conceptually equivalent parts of a representation can be implemented through several mechanisms such a virtual inheritance. We detail here how to implement such switch efficiently using the rascal::RepresentationManagerSortedCoulomb as an example.

The implementation of these two behaviour is encapsulated in the <rascal::internal::SortCoulomMatrix> class and the choice between them is done with a template parameter using template specialization. Note that a in this particular case templated functions could be sufficient but to underline how to implement the most general case a class is used.
~~~~~~~~~~~~~{.hh}
// Enum class defining the several possible sorting options of the Coulomb
// Matrix
enum class CMSortAlgorithm {
  Distance,  // sort according to the distance from the central atom
  RowNorm,   // sort according to the norm of the coulomb matrix rows
};

// Empty general template class implementing the determination of the
// sorting order for the coulomb matrix. It should never be instantiated.
template <CMSortAlgorithm Method>
struct SortCoulomMatrix {};
~~~~~~~~~~~~~
The specific implementation of the two options is done in with template specialization
~~~~~~~~~~~~~{.hh}
/**
 * Sort the coulomb matrix using the distance to the central atom
 * as reference order.
 *
 * @params distance_mat distance matrix between all the atoms in the
 *                      neighbourhood
 */
template <>
struct SortCoulomMatrix<CMSortAlgorithm::Distance> {
  static decltype(auto) get_coulomb_matrix_sorting_order(
      const Eigen::Ref<const Eigen::MatrixXd> & distance_mat,
      const Eigen::Ref<const Eigen::MatrixXd> &) {
    // initialize the distances to be sorted. the center is always first
    std::vector<double> distances_to_sort{0};
    distances_to_sort.reserve(distance_mat.cols());

    for (auto idx_i{1}; idx_i < distance_mat.cols(); ++idx_i) {
      distances_to_sort.push_back(distance_mat(idx_i, 0));
    }

    // find the sorting order
    std::vector<std::pair<size_t, distiter>> order_coulomb(
        distances_to_sort.size());
    size_t nn{0};
    for (distiter it{distances_to_sort.begin()};
         it != distances_to_sort.end(); ++it, ++nn) {
      order_coulomb[nn] = make_pair(nn, it);
    }

    // use stable sort
    std::stable_sort(order_coulomb.begin(), order_coulomb.end(),
                     ordering::ascending);
    return order_coulomb;
  }
};

/**
 * Sort the coulomb matrix using the distance to the central atom
 * as reference order.
 *
 * @params coulomb_mat coulomb matris between all the atoms in the
 *                      neighbourhood
 */
template <>
struct SortCoulomMatrix<CMSortAlgorithm::RowNorm> {
  static decltype(auto) get_coulomb_matrix_sorting_order(
      const Eigen::Ref<const Eigen::MatrixXd> &,
      const Eigen::Ref<const Eigen::MatrixXd> & coulomb_mat) {
    // initialize the distances to be sorted. the center is always first
    std::vector<double> distances_to_sort{};
    distances_to_sort.reserve(coulomb_mat.cols());

    auto row_norms = coulomb_mat.colwise().squaredNorm().eval();
    row_norms(0) = 1e200;
    for (auto idx_i{0}; idx_i < coulomb_mat.cols(); ++idx_i) {
      distances_to_sort.push_back(row_norms(idx_i));
    }

    std::vector<std::pair<size_t, distiter>> order_coulomb(
        distances_to_sort.size());
    size_t nn{0};
    for (distiter it{distances_to_sort.begin()};
         it != distances_to_sort.end(); ++it, ++nn) {
      order_coulomb[nn] = make_pair(nn, it);
    }

    // use stable sort
    std::stable_sort(order_coulomb.begin(), order_coulomb.end(),
                     ordering::descending);

    return order_coulomb;
  }
};
~~~~~~~~~~~~~
Finally, the switch between the two behaviours is done in the <rascal::RepresentationManagerSortedCoulomb::compute()> by calling the templated function <rascal::RepresentationManagerSortedCoulomb::compute_helper()> where the computation of the representation is actually implemented
~~~~~~~~~~~~~{.hh}
template <class Mngr>
void RepresentationManagerSortedCoulomb<Mngr>::compute() {
  auto option{this->options["sorting_algorithm"]};

  if (option == "distance") {
    compute_helper<internal::CMSortAlgorithm::Distance>();
  } else if (option == "row_norm") {
    compute_helper<internal::CMSortAlgorithm::RowNorm>();
  } else {
    auto error_message{std::string("Option '") + option +
                       std::string("' is not implemented.")};
    throw std::invalid_argument(error_message.c_str());
  }
}
~~~~~~~~~~~~~

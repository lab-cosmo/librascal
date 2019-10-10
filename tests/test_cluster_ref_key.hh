#ifndef TESTS_TEST_CLUSTER_REF_KEY_HH_
#define TESTS_TEST_CLUSTER_REF_KEY_HH_

#include "tests.hh"
#include <tuple>
#include "structure_managers/cluster_ref_key.hh"

/**
 * Prints the index sequence of an std::index_sequence type by recursively
 * printing the head element.
 *
 * @param _ the std::index_sequnce to print
 */
namespace rascal {
  struct LayerFixture {
    // {1}
    using order_one_seq = std::index_sequence<1>;
    // {1,1}
    using order_two_seq = std::index_sequence<1, 1>;
    // {1,1,0}
    using order_three_seq = std::index_sequence<1, 1, 0>;
  };

  /**
   * Prints the index sequence defined in the std::index_sequence type of the
   * parameter `sequence` by recursively printing the head element.
   *
   * @param sequence sequence to be printed
   */
  template <size_t Head>
  void print_index_sequence(std::index_sequence<Head> /*sequence*/) {
    std::cout << Head << std::endl;
  }

  template <size_t Head, size_t... Is>
  void print_index_sequence(std::index_sequence<Head, Is...> /*sequence*/) {
    std::cout << Head << ", ";
    print_index_sequence(std::index_sequence<Is...>{});
  }

  /**
   * Prints the index sequence defined in the std::index_sequence type of the
   * parameter `sequence` up to template parameter `Length`
   *
   * @tparam Length
   * @param sequence sequence to be printed
   */
  template <size_t Length, size_t... Ints>
  static void print_index_sequence(std::index_sequence<Ints...> /*sequence*/) {
    std::cout << std::array<size_t, Length>{Ints...} << std::endl;
  }
}  // namespace rascal
#endif // TESTS_TEST_CLUSTER_REF_KEY_HH_

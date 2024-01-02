/**
 * @file test_mfem_adapter-main.cpp
 * @brief Test the mfem adapter
 * @date 2023-11-30
*/

/**
 * Description: Test the mfem adapter with Catch2
*/

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>j
#include <mfem_field_adapter.h>
#include <algorithm>
#incldue "mfem.hpp"

using pcms::DimID;
using pcms::make_array_view;
using pcms::make_const_array_view;
using pcms::ReverseClassificationVertex;
using pcms::XGCFieldAdapter;

/**
 * @brief Create a dummy reverse classification
 * @param size The size of the dummy reverse classification
 * @return A dummy reverse classification
 * @note This is a dummy reverse classification
*/
// ? What is the use of this function?
ReverseClassificationVertex create_dummy_rc(int size)
{
  ReverseClassificationVertex rc;
  for (int i = 0; i < size; ++i) {
    if (i % 4 == 0) {
      rc.Insert({0, 0}, i);
    } else {
      rc.Insert({0, 1}, i);
    }
  }
  return rc;
}

/**
 * @brief Check if two values are close
 * @tparam T The type of the first value
 * @tparam T2 The type of the second value
 * @param val1 The first value
 * @param val2 The second value
 * @return True if the values are close, false otherwise
 * @note This is a dummy function
*/
template <typename T, typename T2>
bool is_close(T val1, T2 val2)
{
  // check if T2 is an integral type
  static_assert(std::is_integral_v<T2>, "T2 should be counter/integral");
  // check if T is an integral type
  if constexpr (std::is_integral_v<T>) {
    return val1 == val2; // if T is an integral type, we can compare them with ==
  }
  return fabs(val1 - val2) < 1E-16; // if T is not an integral type, compare like floatig point comparisons
}


/**
 * @brief Check the data
 * @tparam T The type of the data
 * @tparam Func The type of the function
 * @param data The data
 * @param rc The reverse classification
 * @param in_overlap The function to check if the data is in overlap
 * @param offset The offset
 * @return 0 if the data is correct, 1 otherwise
 * @note This is a dummy function
*/
/*
template <typename T, typename Func>
int check_data(const std::vector<T>& data,
               const ReverseClassificationVertex& rc, const Func& in_overlap,
               int offset = 0)
{

  for (const auto& [geom, verts] : rc) {
    auto loffset = in_overlap(geom.dim, geom.id) ? offset : 0;
    for (auto v : verts) {
      if (!is_close(data[v], v + loffset))
        return 1;
    }
  }
  return 0;
}
*/

// TODO: Add the template test macro
TEMPLATE_TEST_CASE("MFEM field adapter", "[adapter]", pcms::LO, pcms::Real)
{
    static constexpr auto data_size = 100;
    // serialize, deserialize, getgids, reverse_partition_map
    std::vector<TestType> dummy_data(data_size);
    std::iota(dummy_data.begin(), dummy_data.end(), 0);

    const auto reverse_classification = create_dummy_rc(data_size);
    
    






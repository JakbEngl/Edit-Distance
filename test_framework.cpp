#include <parlay/delayed_sequence.h>
#include <parlay/internal/get_time.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/utilities.h>
#include <sys/resource.h>
#include <unistd.h>

#include <fstream>

#include "dac_mm.h"
#include "dac_mm_k.h"
#include "edit_distance_block_hashing.h"
#include "edit_distance_dp.h"
#include "edit_distance_hashing.h"
#include "edit_distance_parallel.h"
#include "edit_distance_Bidirectional_Hash.h"
#include "edit_distance_SA_Pruning.h"
#include "minimum_edit_distance.h"
#include "range_min.h"
#include "suffix_array_parallel.h"
#include "minimum_edit_distance2.h"

//A*
#include "A_editDistance_Block_hashing.h"
#include "A_editDistance_Hashing.h"
#include "A_editDistance_Parallel.h"
//SecondBidirectional
#include "Bidirectional2_edit_distance_block_hash.h"
#include "Bidirectional2_edit_distance_hash.h"
#include "Bidirectional2_edit_distance_SA.h"
//SequentialBidirectional
#include "edit_distance_Bidirectional_Hash.h"
#include "edit_distance_Bidirectional_SA.h"
#include "edit_distance_Bidirerctional_block_hash.h"

constexpr size_t NUM_TESTS = 13;
size_t num_rounds = 3;

template <typename T>
auto generate_strings(size_t n, size_t k, size_t alpha, size_t seed = 0)
{
  //printf("Generating test case... (n: %zu, k: %zu, alpha: %zu)\n", n, k, alpha);
  parlay::sequence<T> A(n), B(n), A2(n), B2(n);
  /*parlay::parallel_for(0, n, [&](size_t i) {
    A[i] = B[i] = (parlay::hash32(i + seed) % 100 < 99) ? 0 : 1;
  });*/

  parlay::parallel_for(
      0, n, [&](size_t i)
      { A[i] = B[i] = parlay::hash32(i + seed) % alpha; });

  size_t _k = k / 3;

  parlay::parallel_for(0, _k, [&](size_t i) {
    size_t idx = parlay::hash32(i + seed) % n;
    B[idx] = parlay::hash32(idx + n + seed) % alpha;
  });

   // insertions and deletions
  auto pred1 = parlay::delayed_seq<bool>(
      n, [&](size_t i)
      { return parlay::hash32_2(i + seed) % n >= _k; });
  auto pred2 = parlay::delayed_seq<bool>(
      n, [&](size_t i)
      { return parlay::hash32_2(i + n + seed) % n >= _k; });
  A = pack(A, pred1);
  B = pack(B, pred2);
  return std::make_tuple(A, B);
}

/*auto generate_strings(size_t n, size_t k, size_t alpha, size_t seed = 0) {
  printf("Generating test case... (n: %zu, k: %zu, alpha: %zu)\n", n, k, alpha);
  parlay::sequence<T> A(n), B(n);
  // parlay::parallel_for(0, n, [&](size_t i) {
  //   A[i] = B[i] = (parlay::hash32(i + seed) % 100 < 99) ? 0 : 1;
  // });
  parlay::parallel_for(
      0, n, [&](size_t i) { A[i] = B[i] = parlay::hash32(i + seed) % alpha;
      });
  // Fills both sequences with a character based on the alphabet from i + seed -> 
  // THis is what leads to the not precise k
  // substitions, insertions, and deletions and roughly equally distributed
  size_t _k = k / 3;

  // substitutions
  parlay::parallel_for(0, _k, [&](size_t i) {
    size_t idx = parlay::hash32(i + seed) % n;

    //std::lock_guard<std::mutex> lock(mutex);
    B[idx] = parlay::hash32(i + n + seed) % alpha;
  });

  // insertions and deletions
  auto pred1 = parlay::delayed_seq<bool>(
      n, [&](size_t i) { return parlay::hash32_2(i + seed) % n >= _k; });
  auto pred2 = parlay::delayed_seq<bool>(
      n, [&](size_t i) { return parlay::hash32_2(i + n + seed) % n >= _k; });
  A = pack(A, pred1);
  B = pack(B, pred2);
  return std::make_tuple(A, B);
}*/

std::string test_name(int id)
{
  switch (id)
  {
  case 0:
    return "BFS-Hash";
    break;
  case 1:
    return "BFS-B-Hash";
    break;
  case 2:
    return "BFS-SA";
    break;
  case 3:
    return "SequentialDirection-Hash";
    break;
  case 4:
    return "SequentialDirection-B-Hash";
    break;
  case 5:
    return "SequentialDirection-SA";
    break;
  case 6:
    return "ParallelDirections-Hash";
    break;
  case 7:
    return "ParallelDirections-B-Hash";
    break;
  case 8:
    return "ParallelDirections-SA";
    break;
  case 9:
    return "A*-Hash";
    break;
  case 10:
    return "A*-B-Hash";
    break;
  case 11:
    return "A*-SA";
    break;
  case 12:
    return "DP - Windowed";
    break;
  default:
    abort();
  }
}

template <typename T>
double test(const parlay::sequence<T> &A, const parlay::sequence<T> &B,
            int id)
{
  //std::cout << "\nTest name: " << test_name(id) << std::endl;
  double total_time = 0;
  double building_time_total = 0;
  for (size_t i = 0; i <= num_rounds; i++)
  {
    parlay::internal::timer t;
    double b_time;
    size_t num_edits;
    switch (id)
    {
    case 0:
      num_edits = EditDistanceHashParallel(A, B, &b_time);
      break;
    case 1:
      num_edits = EditDistanceBlockHashParallel(A, B, &b_time);
      break;
    case 2:
      num_edits = EditDistanceSA(A, B, &b_time);
      break;
    case 3:
      num_edits = EditDistanceHashSeqBidirectional(A, B, &b_time);
      break;
    case 4:
      num_edits = EditDistanceBlockHashBidirectional(A, B, &b_time);
      // num_edits = DAC_MM<sequence<uint32_t>>(A, B).solve();
      break;
    case 5:
      num_edits = EditDistanceSA_Bidirectional1(A, B, &b_time, false);
      // num_edits = EditDistanceDP(A, B);
      break;
    case 6:
      num_edits = bidirectional2_edit_distance_hash(A, B, &b_time);
      // num_edits = minimum_edit_distance(A, B);
      break;
    case 7:
      num_edits = bidirectional2_edit_distance_block_hash(A, B, &b_time);
      break;
    case 8:
      num_edits = bidirectional2_edit_distance_SA(A,B, &b_time, false);
      // num_edits = EditDistanceBidirectional(A, B, &b_time);
      break;
    case 9:
      num_edits = EditDistanceHashParallel_A(A, B, &b_time);
      // num_edits = EditDistancePruning(A, B, &b_time);
      break;
    case 10:
      num_edits = EditDistanceBlockHashParallel_A(A, B, &b_time);
      break;
  case 11:
      num_edits = EditDistanceSA_A(A, B, &b_time, false);
      break;
  case 12:
      num_edits = EditDistanceDP(A, B);
      break;
    default:
      assert(0);
    }
    t.stop();
    if (i == 0)
    {
      printf("#edits: %zu\n", num_edits);
      printf("Warmup round: %f\n", t.total_time());
    }
    else
    {
      printf("Round %zu: %f\n", i, t.total_time());
      total_time += t.total_time();
      building_time_total += b_time;
    }
  }
  double average_time = total_time / num_rounds;
  printf("Average time: %f\n", total_time / num_rounds);
  printf("Average Building time: %f\n", building_time_total / num_rounds);
  //printf("%f\n", total_time / num_rounds);
  //printf("%f\n", building_time_total / num_rounds);

  return average_time;
}

template <typename T>
void run_all(const parlay::sequence<T> &A, const parlay::sequence<T> &B,
             int id = -1)
{
  std::vector<double> times;
  if (id == -1)
  {
    for (size_t i = 0; i < NUM_TESTS; i++)
    {
      times.push_back(test(A, B, i));
    }
  }
  else
  {
    times.push_back(test(A, B, id));
  }
  std::ofstream ofs("edit_distance.tsv", std::ios_base::app);
  for (auto t : times)
  {
    ofs << t << '\t';
  }
  ofs << '\n';
  ofs.close();
}

int main(int argc, char *argv[])
{
  int id = -1;
  size_t n = 1000000;
  size_t k = 1000;
  size_t alpha = n;
  if (argc == 1)
  {
    printf(
        "Usage: ./edit_distance <id> <n> <k> <alpha> <rounds>\n"
        "id: id of the algorithm\n"
        "n: length of strings\n"
        "k: estimated number of edits\n"
        "alpha: alphabet size\n"
        "rounds: number of rounds");
    exit(0);
  }
  if (argc >= 2)
  {
    id = atoi(argv[1]);
  }
  if (argc >= 3)
  {
    n = atoi(argv[2]);
  }
  if (argc >= 4)
  {
    k = atoi(argv[3]);
  }
  if (argc >= 5)
  {
    alpha = atoi(argv[4]);
  }
  if (argc >= 6)
  {
    num_rounds = atoi(argv[5]);
  }
  using Type = uint32_t;
  parlay::sequence<Type> A, B;
  std::tie(A, B) = generate_strings<Type>(n, k, alpha);

  // size_t mem_usage = 0;
  // struct rusage usage;
  // getrusage(RUSAGE_SELF, &usage);
  // mem_usage = usage.ru_maxrss;
  // std::cout << "Memory usage before allocation: " << (mem_usage >> 20) <<
  // std::endl;

  run_all(A, B, id);
  /* BSD, Linux, and OSX -------------------------------------- */
  // struct rusage prusage;
  // getrusage(RUSAGE_SELF, &prusage);
  // return (size_t)(prusage.ru_maxrss * 1024L);
  // cout << "peak mem usage: " << (prusage.ru_maxrss * 1024L >> 30) <<
  // std::endl; get memory usage info for (size_t i = 1; i <= 3000; i++) {
  //   for (size_t j = 3 * i; j >= 1; j -= 3) {
  //     printf("i: %zu, j: %zu\n", i, j);
  //     parlay::sequence<Type> A, B;
  //     std::tie(A, B) = generate_strings<Type>(i, j, 3 * i);
  //     // printf("A.size(): %zu, B.size(): %zu\n", A.size(), B.size());
  //     // printf("A: ");
  //     // for (size_t k = 0; k < A.size(); k++) {
  //     // printf("%u ", A[k]);
  //     //}
  //     // puts("");
  //     // printf("B: ");
  //     // for (size_t k = 0; k < B.size(); k++) {
  //     // printf("%u ", B[k]);
  //     //}
  //     // puts("");
  //     double b_time;
  //     size_t v1 = EditDistanceHashParallel(A, B, &b_time);
  //     size_t v2 = EditDistanceRollingHash(A, B, &b_time);
  //     if (v1 != v2) {
  //       printf("v1: %zu, v2: %zu\n", v1, v2);
  //       printf("wrong answer\n");
  //       return 0;
  //     }
  //   }
  // }

  // parlay::sequence<int> a(1000), b(1000);
  // for (int i = 0; i < 1000; i++) {
  //   a[i] = rand() % 1000;
  //   b[i] = rand() % 1000;
  // }
  // auto n = a.size(), m = b.size();

  // auto c = parlay::sequence<uint32_t>(n + m + 1);
  // parlay::parallel_for(0, n, [&](int i) { c[i] = a[i]; });
  // parlay::parallel_for(0, m, [&](int i) { c[i + n + 1] = b[i]; });
  // c[n] = std::numeric_limits<uint32_t>::max();
  // auto rank = parlay::sequence<unsigned int>();
  // auto sa = parlay::sequence<unsigned int>();
  // auto lcp = parlay::sequence<unsigned int>();
  // std::tie(rank, sa, lcp) = suffix_array_large_alphabet(c);
  // auto rmq = range_min(lcp);
  // auto GetLcp = [&](int i, int j) -> int {
  //   // std::cout << "GetLcp " << i << ' ' << j << '\n';
  //   if (i == n || j == m) return 0;
  //   assert(0 <= i && i < n && 0 <= j && j < m);
  //   for (int k = 0; k < 8; k++) {
  //     if (i + k >= n || j + k >= m || a[i + k] != b[j + k]) {
  //       return k;
  //     }
  //   }
  //   int l = rank[i], r = rank[j + n + 1];
  //   if (l > r) std::swap(l, r);
  //   assert(l < r);
  //   int id = rmq.query(l + 1, r);
  //   return lcp[id];
  // };

  // parlay::sequence<parlay::sequence<int>> table_A;
  // parlay::sequence<parlay::sequence<int>> table_B;
  // parlay::sequence<std::pair<int, int>> S_A;
  // parlay::sequence<std::pair<int, int>> S_B;
  // construct_table(a, b, table_A, table_B, S_A, S_B, std::min(n, m));
  // auto res = block_query_lcp(0, 0, a, b, table_A, table_B, S_A, S_B);
  // // cout << "res: " << res << endl;

  // for (int i = 0; i < 100; i++) {
  //   int p = rand() % n;
  //   int q = rand() % m;
  //   int res1 = GetLcp(p, q);
  //   int res2 = block_query_lcp(p, q, a, b, table_A, table_B, S_A, S_B);
  //   if (res1 != res2) {
  //     cout << "error: " << p << ' ' << q << endl;
  //     abort();
  //   }
  // }

  // cout << "good" << endl;

  return 0;
}

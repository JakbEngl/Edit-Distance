#include "rolling_hashing.h"
<<<<<<< HEAD
=======
// using hash_r_T = int32_t;
>>>>>>> cf8b1816a2cd43c5c6de7a3764204cd560128de0

// Returns an 32-bit integer of edit distance for two sequences A and B
// using bfs and hashing table query for lcp, in parallel
template <typename Seq>
int EditDistanceRollingHash(const Seq &a, const Seq &b, double *building_tm) {
  // build rolling hashing table with value s[0..i]

  /**
   * Record of building time
   */
  Timer tmr;
  int n = a.size();
  int m = b.size();
  if (n == 0) {
    return m;
  }
  if (m == 0) {
    return n;
  }
  parlay::sequence<hash_r_T> table_s1;
  parlay::sequence<hash_r_T> table_s2;
  // size_t mem_usage = 0;
  // struct rusage usage;
  // getrusage(RUSAGE_SELF, &usage);
  // mem_usage = usage.ru_maxrss;
  // std::cout << "Memory usage before allocation: " << mem_usage << std::endl;
  build_rolling(a, b, table_s1, table_s2);

  // getrusage(RUSAGE_SELF, &usage);
  // mem_usage = usage.ru_maxrss;
  // std::cout << "Memory usage after allocation: " << mem_usage << std::endl;
  auto Diag = [&](int i, int j) { return i - j + m; };
  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);

  *building_tm = tmr.elapsed();
  max_row[Diag(0, 0)] = query_rolling(a, b, table_s1, table_s2, 0, 0);
  int k = 0;
  for (;;) {
    if (max_row[Diag(n, m)] == n) break;  // find path
    k++;
    int l = Diag(0, std::min(k, int(m)));
    int r = Diag(std::min(k, int(n)), 0);
    parlay::parallel_for(l, r + 1, [&](int id) {
      int t = -1;
      if (max_row[id] != -1) {
        int i = max_row[id];
        int j = i + m - id;
        if (i == n || j == m) {
          t = i;
        } else {
          int get_lcp = query_rolling(a, b, table_s1, table_s2, i + 1, j + 1);
          t = i + 1 + get_lcp;
        }
      }
      if (id > 0 && max_row[id - 1] != -1) {
        int i = max_row[id - 1];
        int j = i + m - id + 1;
        if (i == n) {
          t = n;
        } else {
          t = std::max(
              t, i + 1 + query_rolling(a, b, table_s1, table_s2, i + 1, j));
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          t = std::max(t,
                       i + query_rolling(a, b, table_s1, table_s2, i, j + 1));
        }
      }
      // assert(t <= n);
      temp[id] = t;
    });
    parlay::parallel_for(l, r + 1,
                         [&](int id) { max_row[id] = std::min(temp[id], id); });
  }
  return k;
}

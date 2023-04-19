#include "rolling_hashing.h"
using hash_r_T = int32_t;

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
  // cout << "seq" << endl;
  // for (int i = 0; i < 20; i++) {
  //   cout << a[i] << " ";
  // }
  // cout << endl;
  build_rolling(a, b, table_s1, table_s2);
  // cout << "rolling table 1: " << endl;
  // for (int i = 1; i < 50; i++) {
  //   // hash_r_T res = table_s1[i * 4];
  //   // for (int j = i * 4 + 1; j < (i + 1) * 4; j++) {
  //   //   res = res * PRIME + table_s1[j];
  //   // }
  //   cout << table_s1[i * 4 - 1] << " ";
  // }
  // cout << endl;
  // cout << "rolling table 2: " << endl;
  // for (int i = 1; i < 50; i++) {
  //   // hash_r_T res = table_s2[i * 4];
  //   // for (int j = i * 4 + 1; j < (i + 1) * 4; j++) {
  //   //   res = res * PRIME + table_s2[j];
  //   // }
  //   // cout << res << " ";
  //   cout << table_s2[i * 4 - 1] << " ";
  // }
  // cout << endl;

  auto Diag = [&](int i, int j) { return i - j + m; };
  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);

  *building_tm = tmr.elapsed();
  max_row[Diag(0, 0)] = query_rolling(a, b, table_s1, table_s2, 0, 0);
  // std::cout << "rolling lcp: "
  //           << query_rolling(a, b, table_s1, table_s2, 1000, 1000) << std::endl;

  // assert(query_rolling(a, b, table_s1, table_s2, 0, 0) == test_lcp(a, b,
  // 0, 0)); bfs for path
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
#include "hash_block_parallel.h"
#include "rolling_hashing_blk.h"
#include <chrono>
#include <iostream>
#include <queue>
#include <mutex>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>
#include "utils.h"
#include <queue>
#include <unordered_map>

using hash_T = int64_t;
using hash_r_b_T = uint64_t;
// Returns an 32-bit integer of edit distance for two sequences A and B
// using bfs and hashing table query for lcp, in parallel
template <typename Seq>
int EditDistanceBlockHashParallel_A(const Seq &a, const Seq &b,
                                  double *building_tm) {
  parlay::internal::timer tt;
  // build sparse table

  /**
   * rocord of building time
   */
  Timer tmr;

  struct Cell
  {
    int i, j;  // Position in the DAG
    int edits; // number of edits
  };

  int n = a.size();
  int m = b.size();
  if (n == 0) return m;
  if (m == 0) return n;
  parlay::sequence<hash_r_b_T> table_A;
  parlay::sequence<hash_r_b_T> table_B;
  //parlay::sequence<parlay::sequence<hash_T>> table_A;
  //parlay::sequence<parlay::sequence<hash_T>> table_B;
  //parlay::sequence<parlay::sequence<HashString>> table_A;
  //parlay::sequence<parlay::sequence<HashString>> table_B;
  //parlay::sequence<std::pair<hash_T, hash_T>> S_A;
  //parlay::sequence<std::pair<hash_T, hash_T>> S_B;
  //construct_table(a, b, table_A, table_B, S_A, S_B, std::min(n, m));
  build_rolling_blk(a,b,table_A, table_B);
  std::vector<int> furthestProgress[3] = {std::vector<int>(n + m + 1, -1), std::vector<int>(n + m + 1, -1),
                                          std::vector<int>(n + m + 1, -1)};
  std::queue<Cell> queues[3];

  auto Diag = [&](int i, int j) { return i - j + m; };
  auto Heuristic = [&](int i, int j)
  { return std::abs(Diag(n, m) - Diag(i, j)); };
  //parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);

  *building_tm = tmr.elapsed();
  //tt.next("building");
  // std::cout << building_tm << ", ";
  //durat
  queues[0].push({0, 0, 0}); 
  size_t result = std::numeric_limits<size_t>::max();

  int activeQueue = 0;

  auto updateCellAndQueue = [&](int i, int j, int edits, int cur_heur_value, int operation)
  {
    // 0 => substitution, 1 => insertion, 2 => deletion
    int new_i = i + (operation != 1 ? 1 : 0); // Øger I medmindre insertion
    int new_j = j + (operation != 2 ? 1 : 0); // Øger j medmindre deletion
    int new_edits = edits + 1;
    int olddiag = Diag(i,j);
    int diag = Diag(new_i, new_j);
    int heuristic_cost = (operation == 0) ? 0 : Heuristic(new_i, new_j) - cur_heur_value;
    //int heuristic_cost = (operation == 0) ? 0 : diag_heur[diag] - diag_heur[olddiag];
    int newQueueIndex = (activeQueue + heuristic_cost + 1) % 3;

    if (furthestProgress[newQueueIndex][diag] < new_i)
    {
      furthestProgress[newQueueIndex][diag] = new_i;
      queues[newQueueIndex].push({new_i, new_j, new_edits});
    }
  };

  while (!queues[0].empty() || !queues[1].empty() || !queues[2].empty())
  {
    if (queues[activeQueue].empty())
    {
      // Hvis køen bliver tom, så gå videre (MAKE IT WRAP AROUND)
      activeQueue = (activeQueue + 1) % 3;
      continue;
    }

    Cell current = queues[activeQueue].front();
    queues[activeQueue].pop();

    //int lcp = block_query_lcp(current.i, current.j, a, b, table_A, table_B, S_A, S_B);
    int lcp = query_rolling_blk(a,b, table_A, table_B, current.i, current.j);
    if (lcp > 0)
    {
      int skipped_i = current.i + lcp;
      int skipped_j = current.j + lcp;
      // Ny current
      Cell skipped_cell = {skipped_i, skipped_j, current.edits};
      current = skipped_cell; // Opdateret current med LCP
    }

    // Slutcondition
    if (current.i == n && current.j == m)
    {
      result = current.edits;
      break; // Success er opnået!
    }

    int cur_heur_value = Heuristic(current.i, current.j);
    int newQueueIndex = activeQueue + 1;

    if (current.i < n && current.j < m)
    { // Substitution (0)
      updateCellAndQueue(current.i, current.j, current.edits, cur_heur_value, 0);
    }
    if (current.j < m)
    { // Insertion (1)
      updateCellAndQueue(current.i, current.j, current.edits, cur_heur_value, 1);
    }
    if (current.i < n)
    { // Deletion (2)
      updateCellAndQueue(current.i, current.j, current.edits, cur_heur_value, 2);
    }
  }
  return result;
}
  
#include "edit_distance_parallel.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include "dc3.h"
#include "parlay/sequence.h"
#include "range_min.h"
#include "sparse_table_sequential.h"
#include "suffix_array_parallel.h"
#include "utils.h"
#include <queue>
#include <unordered_map>

// #define COMPUTE_AVERAGE_LCP

size_t EditDistanceSA_A(const parlay::sequence<uint32_t> &a,
                      const parlay::sequence<uint32_t> &b, double *building_tm,
                      bool use_DC3)
{
  Timer tmr;
#ifdef COMPUTE_AVERAGE_LCP
  std::atomic_uint64_t lcp_total = 0;
  std::atomic_uint64_t lcp_cnt = 0;
#endif
  struct Cell
  {
    int i, j;  // Position in the DAG
    int edits; // number of edits
  };
  int n = a.size(), m = b.size();
  auto c = parlay::sequence<uint32_t>(n + m);
  parlay::parallel_for(0, n, [&](int i)
                       { c[i] = a[i]; });
  parlay::parallel_for(0, m, [&](int i)
                       { c[i + n] = b[i]; });
  auto rank = parlay::sequence<unsigned int>();
  auto sa = parlay::sequence<unsigned int>();
  auto lcp = parlay::sequence<unsigned int>();
  if (use_DC3)
  {
    std::tie(rank, sa, lcp) = DC3(c);
  }
  else
  {
    std::tie(rank, sa, lcp) = suffix_array_large_alphabet(c);
  }
  auto rmq = range_min(lcp);
  auto GetLcp = [&](int i, int j) -> int
  {
    if (i == n || j == m)
      return 0;
    assert(0 <= i && i < n && 0 <= j && j < m);
    for (int k = 0; k < 8; k++)
    {
      if (i + k >= n || j + k >= m || a[i + k] != b[j + k])
      {
        return k;
      }
    }
    int l = rank[i], r = rank[j + n];
    if (l > r)
      std::swap(l, r);
    int id = rmq.query(l + 1, r);
    return std::min(lcp[id], (unsigned int)n - i);
  };

  *building_tm = tmr.elapsed();
  // std::cout << " building time of SA: " << building_tm << std::endl;
  // A* tomfoolery
  // std::unordered_map<int, int> furthestProgress[3];
  std::vector<int> furthestProgress[3] = {std::vector<int>(n + m + 1, -1), std::vector<int>(n + m + 1, -1),
                                          std::vector<int>(n + m + 1, -1)};
  std::queue<Cell> queues[3];

  auto Diag = [&](int i, int j)
  { return i - j + m; };
  auto Heuristic = [&](int i, int j)
  { return std::abs(Diag(n, m) - Diag(i, j)); };
  /*
  //Sparse table for Heuristic værdier
  int endDiag = Diag(n, m);
  std::vector<int> diag_heur = std::vector<int>(n + m + 1, -1);
  parlay::parallel_for(0, diag_heur.size(), [&](int i) {
    diag_heur[i] = std::abs(i - endDiag);
});*/

  // Start med 0 edits
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

    int lcp = GetLcp(current.i, current.j);
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
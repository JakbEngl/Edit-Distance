#include "edit_distance_parallel.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>
#include <chrono>
#include <queue>
#include <mutex>
#include <cmath>

#include "dc3.h"
#include "parlay/sequence.h"
#include "range_min.h"
#include "sparse_table_sequential.h"
#include "suffix_array_parallel.h"
#include "utils.h"

// #define COMPUTE_AVERAGE_LCP

size_t EditDistanceSA(const parlay::sequence<uint32_t> &a,
                      const parlay::sequence<uint32_t> &b, double *building_tm,
                      bool use_DC3)
{
  Timer tmr;
#ifdef COMPUTE_AVERAGE_LCP
  std::atomic_uint64_t lcp_total = 0;
  std::atomic_uint64_t lcp_cnt = 0;
#endif
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

  auto Diag = [&](int i, int j)
  { return i - j + m; };

  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);
  max_row[Diag(0, 0)] = GetLcp(0, 0);
#ifdef COMPUTE_AVERAGE_LCP
  lcp_total += GetLcp(0, 0);
  lcp_cnt++;
#endif
  int k = 0;
  
  /*std::atomic<int> Count0(0);
  std::atomic<int> Count4(0);
  std::atomic<int> Count8(0);
  std::atomic<int> Count16(0);
  std::atomic<int> Count32(0);
  std::atomic<int> Count64(0);
  std::atomic<int> Count128(0);
  std::atomic<int> Count256(0);
  std::atomic<int> CountBig(0);
  std::atomic <long long> durationZero(0);
  std::atomic <long long> duration4 = (0);
  std::atomic <long long> duration8 = (0);
  std::atomic <long long> duration16 = (0);
  std::atomic <long long> duration32 = (0);
  std::atomic <long long> duration64 = (0);
  std::atomic <long long> duration128 = (0);
  std::atomic <long long> duration256 = (0);
  std::atomic <long long> durationBig = (0);*/
  /*
  std::atomic<int> max(0);
  std::priority_queue<int, std::vector<int>, std::greater<int>> minHeap;
  std::mutex heapMutex;
  const size_t MAX_TOP_VALUES = 10;
  std::atomic<int> min(INT_MAX);*/
  //std::atomic<int> lcp_count(0);
  //std::atomic<int> lcp_HIGH_count(0);
  //std::atomic <long long> timingsLow(0);
  //std::atomic <long long> timingsAll(0);
  //std::atomic <long long> stdLow(0);
  //std::atomic <long long> stdAll(0);

  /*std::atomic<int> count50(0);
  std::atomic<int> count60(0);
  std::atomic<int> count70(0);
  std::atomic<int> count80(0);
  std::atomic<int> count140(0);
  std::atomic<int> count260(0);
  std::atomic<int> count500(0);
  std::atomic<int> count1000(0);
  std::atomic<int> countOutliers(0);*/


  for (;;)
  {
    if (max_row[Diag(n, m)] == n)
      break; // find path
    k++;
    int l = Diag(0, std::min(k, m)), r = Diag(std::min(k, n), 0);
    parlay::parallel_for(l, r + 1, [&](int id)
                         {
      int t = -1;
      if (max_row[id] != -1) {
        int i = max_row[id];
        int j = i + m - id;
        if (i == n || j == m) {
          t = i;
        } else {
          //auto start = std::chrono::high_resolution_clock::now();
          int _lcp = GetLcp(i + 1, j + 1);
          //int _lcp_x = GetLcp(i + 1, j + 1);
          //auto stop = std::chrono::high_resolution_clock::now();
          //auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

          /*switch (_lcp){
            case 0:
              Count0++;
              durationZero += duration;
              break;
            case 1 ... 4:
              Count4++;
              duration4 += duration;
              break;
            case 5 ... 8:
              Count8++;
              duration8 += duration;
              break;
            case 9 ... 16:
              Count16++;
              duration16 += duration;
              break;
            case 17 ... 32:
              Count32++;
              duration32 += duration;
              break;
            case 33 ... 64:
              Count64++;
              duration64 += duration;
              break;
            case 65 ... 128:
              Count128++;
              duration128 += duration;
              break;
            case 129 ... 256:
              Count256++;
              duration256 += duration;
              break;
            default:
              CountBig++;
              durationBig += duration;
          }*/

          /* switch (duration){
            case 0 ... 50: 
              count50++;
              break;
            case 51 ... 60: 
              count60++;
              break;
            case 61 ... 70: 
              count70++;
              break;
            case 71 ... 80: 
              count80++;
              break;
            case 81 ... 140: 
              count140++;
              break;
            case 141 ... 260: 
              count260++;
              break;
            case 261 ... 500: 
              count500++;
              break;
            case 501 ... 1000: 
              count1000++;
              break;
            default: 
              countOutliers++;
          }*/

          /*timingsAll += (duration);
          lcp_count++;
          auto deviation = duration - 61;
          auto squaredDeviation = deviation * deviation;
          if (duration < 35 * 61){
            lcp_HIGH_count++;
            timingsLow += duration;
            stdLow += squaredDeviation;
          } 
          stdAll += squaredDeviation;*/

          /*
          std::lock_guard<std::mutex> lock(heapMutex);
          if (minHeap.size() < MAX_TOP_VALUES) {
              minHeap.push(duration);
          } else if (duration > minHeap.top()) {
              minHeap.pop();
              minHeap.push(duration);
            }
  
          int currentMin = min.load();
            while (duration < currentMin) {
              if (min.compare_exchange_weak(currentMin, duration)) {
                break;
                }
              }*/
          
          //lcp_count++;
          //summedDuration += duration;
          t = i + 1 + _lcp;
#ifdef COMPUTE_AVERAGE_LCP
          lcp_total += _lcp;
          lcp_cnt++;
#endif
        }
      }
      if (id > 0 && max_row[id - 1] != -1) {
        int i = max_row[id - 1];
        int j = i + m - id + 1;
        if (i == n) {
          t = n;
        } else {
          //auto start = std::chrono::high_resolution_clock::now();
          int _lcp_2 = GetLcp(i + 1, j);
          //int _lcp_2_x = GetLcp(i + 1, j);
          //auto stop = std::chrono::high_resolution_clock::now();
          //auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

          /*switch (_lcp_2){
            case 0:
              Count0++;
              durationZero += duration;
              break;
            case 1 ... 4:
              Count4++;
              duration4 += duration;
              break;
            case 5 ... 8:
              Count8++;
              duration8 += duration;
              break;
            case 9 ... 16:
              Count16++;
              duration16 += duration;
              break;
            case 17 ... 32:
              Count32++;
              duration32 += duration;
              break;
            case 33 ... 64:
              Count64++;
              duration64 += duration;
              break;
            case 65 ... 128:
              Count128++;
              duration128 += duration;
              break;
            case 129 ... 256:
              Count256++;
              duration256 += duration;
              break;
            default:
              CountBig++;
              durationBig += duration;
          }*/

          /* switch (duration){
            case 0 ... 50: 
              count50++;
              break;
            case 51 ... 60: 
              count60++;
              break;
            case 61 ... 70: 
              count70++;
              break;
            case 71 ... 80: 
              count80++;
              break;
            case 81 ... 140: 
              count140++;
              break;
            case 141 ... 260: 
              count260++;
              break;
            case 261 ... 500: 
              count500++;
              break;
            case 501 ... 1000: 
              count1000++;
              break;
            default: 
              countOutliers++;
          }*/

          /*timingsAll += duration;
          lcp_count++;
          auto deviation = duration - 61;
          auto squaredDeviation = deviation * deviation;
          if (duration < 35 * 61){
            lcp_HIGH_count++;
            timingsLow += duration;
            stdLow += squaredDeviation;
          } 
          stdAll += squaredDeviation;*/
          //auto deviation = duration.count() - 61;
          //auto squaredDeviation = deviation * deviation;
          //summedSquaredDeviation += squaredDeviation;
          /*
          std::lock_guard<std::mutex> lock(heapMutex);
          if (minHeap.size() < MAX_TOP_VALUES) {
              minHeap.push(duration);
          } else if (duration > minHeap.top()) {
              minHeap.pop();
              minHeap.push(duration);
            }

          int currentMin = min.load();
            while (duration < currentMin) {
              if (min.compare_exchange_weak(currentMin, duration)) {
                break;
                }
              }*/
          
          //lcp_count++;
          //summedDuration += duration;
          t = std::max(t, i + 1 + _lcp_2);
#ifdef COMPUTE_AVERAGE_LCP
          lcp_total += _lcp_2;
          lcp_cnt++;
#endif
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          //auto start = std::chrono::high_resolution_clock::now();
          int _lcp_3 = GetLcp(i, j + 1);
          //int _lcp_3_x = GetLcp(i, j + 1);
          //auto stop = std::chrono::high_resolution_clock::now();
          //auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();

          /*switch (_lcp_3){
            case 0:
              Count0++;
              durationZero += duration;
              break;
            case 1 ... 4:
              Count4++;
              duration4 += duration;
              break;
            case 5 ... 8:
              Count8++;
              duration8 += duration;
              break;
            case 9 ... 16:
              Count16++;
              duration16 += duration;
              break;
            case 17 ... 32:
              Count32++;
              duration32 += duration;
              break;
            case 33 ... 64:
              Count64++;
              duration64 += duration;
              break;
            case 65 ... 128:
              Count128++;
              duration128 += duration;
              break;
            case 129 ... 256:
              Count256++;
              duration256 += duration;
              break;
            default:
              CountBig++;
              durationBig += duration;
          }*/

          /*switch (duration){
            case 0 ... 50: 
              count50++;
              break;
            case 51 ... 60: 
              count60++;
              break;
            case 61 ... 70: 
              count70++;
              break;
            case 71 ... 80: 
              count80++;
              break;
            case 81 ... 140: 
              count140++;
              break;
            case 141 ... 260: 
              count260++;
              break;
            case 261 ... 500: 
              count500++;
              break;
            case 501 ... 1000: 
              count1000++;
              break;
            default: 
              countOutliers++;
          }
          */
          /*timingsAll += duration;
          lcp_count++;
          auto deviation = duration - 61;
          auto squaredDeviation = deviation * deviation;
          if (duration < 35 * 61){
            lcp_HIGH_count++;
            timingsLow += duration;
            stdLow += squaredDeviation;
          } 
          stdAll += squaredDeviation;*/
          //auto deviation = duration.count() - 61;
          //auto squaredDeviation = deviation * deviation;
          //summedSquaredDeviation += squaredDeviation;
          /*
          std::lock_guard<std::mutex> lock(heapMutex);
          if (minHeap.size() < MAX_TOP_VALUES) {
              minHeap.push(duration);
          } else if (duration > minHeap.top()) {
              minHeap.pop();
              minHeap.push(duration);
            }

          int currentMin = min.load();
            while (duration < currentMin) {
              if (min.compare_exchange_weak(currentMin, duration)) {
                break;
                }
              }*/
          
        

          t = std::max(t, i + _lcp_3);
#ifdef COMPUTE_AVERAGE_LCP
          lcp_total += _lcp_3;
          lcp_cnt++;
#endif
        }
      }
      assert(t <= n);
      temp[id] = t; });
    parlay::parallel_for(l, r + 1,
                         [&](int id)
                         { max_row[id] = std::min(temp[id], id); });
  }
#ifdef COMPUTE_AVERAGE_LCP
  std::cout << "Lcp total: " << lcp_total << std::endl;
  std::cout << "lcp_cnt: " << lcp_cnt << std::endl;
  std::cout << "average: " << (double)(lcp_total) / lcp_cnt << std::endl;
#endif
  /*std::cout << "Duration for zeroes: " << durationZero.count()/count << " Count: " << count << std::endl;
  std::cout << "Duration for small: " << durationSmall.count()/Smallcount << " Count small: " << Smallcount << std::endl;
  std::cout << "Duration for mid: " << durationMid.count()/Midcount << " Count mid: " << Midcount << std::endl;
  std::cout << "Duration for large: " << durationLarge.count()/Largecount << " Count large: " << Largecount << std::endl;*/
  //std::cout << "Duration for LCP average: " << summedDuration.count()/lcp_count << " Count LCP: " << lcp_count << std::endl;

  /*std::cout << "Min: " << min << std::endl;
  std::priority_queue<int, std::vector<int>, std::greater<int>> copy = minHeap;
  std::cout << "Top durations are: ";
  while (!copy.empty()) {
        std::cout << copy.top() << " ns, ";
        copy.pop();
    }
    std::cout << std::endl;*/
  /*double averageNoOutliers = (timingsLow) / lcp_HIGH_count;
  double averageAll = (timingsAll) / lcp_count;
  double varianceAll = stdAll / (lcp_count - 1);
  double varianceLow = stdLow / (lcp_HIGH_count - 1);
  double standardDeviationAll = sqrt(varianceAll);
  double standardDeviationLow = sqrt(varianceLow);

  std::cout << "Count by all: " << lcp_count << " Nanoseconds" << std::endl;
  std::cout << "Count without outliers: " << lcp_HIGH_count << " Nanoseconds" << std::endl;
  std::cout << "Total time taken by all: " << timingsAll << " Nanoseconds" << std::endl;
  std::cout << "Total time taken without outliers: " << timingsLow << " Nanoseconds" << std::endl;
  std::cout << "Average time taken by all: " << averageAll << " Nanoseconds" << std::endl;
  std::cout << "Average time taken without outliers: " << averageNoOutliers << " Nanoseconds" << std::endl;
  std::cout << "The standard deviation for all: " << standardDeviationAll << " Nanoseconds" << std::endl;
  std::cout << "The standard deviation without outliers: " << standardDeviationLow << " Nanoseconds" << std::endl;
  */

  /*double average0;
  double average4;
  double average8;
  double average16;
  double average32;
  double average64;
  double average128;
  double average256;
  double averageBig;
  if (Count0 == 0){
    average0 = 0.0;
  }
  else {average0 = durationZero/Count0;}

  if (Count4 == 0){
    average4 = 0.0;
  }
  else {average4 = duration4/Count4;}

  if (Count8 == 0){
    average8 = 0.0;
  }
  else {average8 = duration8/Count8;}

  if (Count16 == 0){
    average16 = 0.0;
  }
  else {average16 = duration16/Count16;}

  if (Count32 == 0){
    average32 = 0.0;
  }
  else {average32 = duration32/Count32;}

  if (Count64 == 0){
    average64 = 0.0;
  }
  else {average64 = duration64/Count64;}

  if (Count128 == 0){
    average128 = 0.0;
  }
  else {average128 = duration128/Count128;}

  if (Count256 == 0){
    average256 = 0.0;
  }
  else {average256 = duration256/Count256;}

  if (CountBig == 0){
    averageBig = 0.0;
  }
  else {averageBig = durationBig/CountBig;}

  std::cout << Count0 << std::endl;
  std::cout << Count4 << std::endl;
  std::cout << Count8 << std::endl;
  std::cout << Count16 << std::endl;
  std::cout << Count32 << std::endl;
  std::cout << Count64 << std::endl;
  std::cout << Count128 << std::endl;
  std::cout << Count256 << std::endl;
  std::cout << CountBig << std::endl;

  std::cout << average0 << std::endl;
  std::cout << average4 << std::endl;
  std::cout << average8 << std::endl;
  std::cout << average16 << std::endl;
  std::cout << average32 << std::endl;
  std::cout << average64 << std::endl;
  std::cout << average128 << std::endl;
  std::cout << average256 << std::endl;
  std::cout << averageBig << std::endl;*/

  /*std::cout << count50 << std::endl;
  std::cout << count60 << std::endl;
  std::cout << count70 << std::endl;
  std::cout << count80 << std::endl;
  std::cout << count140 << std::endl;
  std::cout << count260 << std::endl;
  std::cout << count500 << std::endl;
  std::cout << count1000 << std::endl;
  std::cout << countOutliers << std::endl;*/
  return k;
}

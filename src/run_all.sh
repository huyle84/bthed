#!/bin/bash
# Stage-1 experiment batch. Output: ../results/*.csv
set -e
R=../results
( ./bthed_sim2 stats 1000000 3 11 > $R/stats_mldsa44_be3.csv
  for b in 1 2 7 8; do ./bthed_sim2 stats 200000 $b 12 > $R/stats_mldsa44_be$b.csv; done
  for L in 8 20; do ./bthed_sim2 session 20000 $L 3 13 > $R/session_mldsa44_L$L.csv; done
  ./bthed_sim2 bench 2001 3 > $R/bench_mldsa44.csv ) &
( ./bthed_sim3 stats 300000 3 21 > $R/stats_mldsa65_be3.csv
  ./bthed_sim5 stats 300000 3 31 > $R/stats_mldsa87_be3.csv
  ./bthed_sim3 stats 100000 8 22 > $R/stats_mldsa65_be8.csv
  ./bthed_sim5 stats 100000 8 32 > $R/stats_mldsa87_be8.csv
  for L in 8 20; do ./bthed_sim3 session 5000 $L 3 23 > $R/session_mldsa65_L$L.csv;
                    ./bthed_sim5 session 5000 $L 3 33 > $R/session_mldsa87_L$L.csv; done
  ./bthed_sim3 bench 1001 3 > $R/bench_mldsa65.csv
  ./bthed_sim5 bench 1001 3 > $R/bench_mldsa87.csv ) &
wait
for nt in "3 2" "5 3" "7 4" "10 6" "16 9" "20 11"; do ./bthed_sim2 rss $nt; done > $R/rss.csv
echo DONE

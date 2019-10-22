# TAC_FSRPD

1. Fig3_FSR_and_conf_sets_triangle.m
    - FSR analysis (95% confidence region and support) of
      unicycle dynamics with fixed turning rate sequence and
      velocity as a triangular distributed rv (about 10 with
      end points at 9 and 11).
    - This computation took about 20 seconds
2. Fig6_alg1_vs_alg2.m 
    - Comparison between Algorithms 1 and 2. Use figure
      number (1, 2, or 3) within this script to get other
      plots
3. pursuit_problem.m
    - Target pursuit problem
4. Fig10_avoid_set.m
    - Avoid set corresponding to obstacle as a Dubins
      vehicle --- construct a finite union of convex and
      compact sets

## Requires

- SReachTools
    - MPT3
    - CVX
- CharFunTools
- EllipsoidalToolbox v1.1.3 (older version)

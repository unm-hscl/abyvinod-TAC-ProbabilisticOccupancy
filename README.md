# abyvinod-TAC-ProbabilisticOccupancy

| Title      | Probabilistic Occupancy via Forward Stochastic Reachability for Markov Jump Affine Systems |
|------------|--------------------------------------------------------------------------------------------|
| Authors    | Abraham P. Vinod & Meeko M. K. Oishi                                                       |
| Journal    | IEEE Transactions on Automatic Control                                                     |

To generate the figures associated with the computations, run the following
MATLAB scripts:

1. `Fig3_FSR_and_conf_sets_triangle.m`
    - FSR analysis (95% confidence region and support) of unicycle dynamics with
      fixed turning rate sequence and velocity as a triangular distributed rv
      (about 10 with end points at 9 and 11).
    - This computation took about 20 seconds
2. `Fig6_runall.m` (which internally runs `Fig6_alg1_vs_alg2.m`)
    - Comparison between Algorithms 1 and 2. Use figure number (1, 2, or 3)
      within this script to get other plots
    - This computation took about 2 minutes
3. `Fig789_pursuit_problem.m`
    - Target pursuit problem discussed in Section 6.A
    - This computation took about 8 minutes
4. `Fig10_avoid_set.m`
    - Keep-out set construction discussed in Section 6.B
    - Avoid set corresponding to obstacle as a Dubins vehicle --- construct a
      finite union of convex and compact sets
    - Change the `figure_number` to get different figures

## Requires

- SReachTools https://sreachtools.github.io
    - MPT3 https://www.mpt3.org/
    - CVX http://cvxr.com/cvx/
- CharFunTools https://github.com/witkovsky/CharFunTool
    - All codes were tested on the bleeding edge:
      https://github.com/witkovsky/CharFunTool/commit/959c9ac6612dd521bbe87385913bc871d84a199f

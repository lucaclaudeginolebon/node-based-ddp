# Node-based DDP

This project contains the code used to generate the results in the following papers:

1) L. C. G. Lebon and C. Altafini, Geometric Control Theory Over Networks: Minimal Node Cardinality Disturbance Decoupling Problems, accepted in IEEE Transactions on Automatic Control, 2026, https://doi.org/10.1109/TAC.2026.3701288.
2) L. C. G. Lebon, J. de Vlas, J. K. Fichte, and C. Altafini, A Polynomial-Time Algorithm for the Minimal Input/Output Node Cardinality Disturbance Decoupling Problem via Output Feedback, submitted to ECC 2027, 2026, https://doi.org/10.48550/arXiv.2610.22735.

From 1): In particular, Zmax and Smin from Propositions 4 and 5, the solvability conditions of Theorems 2 and 4, Algorithms 1--3, the computation of the friend B'*A from Theorem 9, the output injection A*C' of Theorem 10, the output feedback B'*A*C' of Theorem 11 and the full and reduced order compensator of Theorem 12, can all be computed using the MATLAB functions provided in this project. Moreover, functions identifying the minimal cardinality input nodes for the DDPSF, the suboptimal input/output allocations for the DDPOF, the minimal cardinality input/output solutions for the DDPDF and the corresponding (C,A,B)-pairs of node sets are also provided.

From 2): function min_card_ddpof.m is the optimal polynomial-time solver for the minimal input/output DDPOF.

The script example_ddp.m contains a demonstration of the main functions implementing 1) and 2).

For additional information on the usage of the functions, use the `help' command in the Command Window of MATLAB.
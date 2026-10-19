% We provide a worked-out example using the main MATLAB functions of the project.
% ========================================================================
% Lebon Luca Claude Gino, October 2026 (Matlab R2022b)
% ========================================================================

clear;
close all;
clc;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%% Setup - select testcase

a = 1; % testcase 1
% a = 2; % testcase 2
% a = 3; % testcase 3 (from Ex.1 of "A Polynomial-Time Algorithm for the Minimal Input/Output Node Cardinality Disturbance Decoupling Problem via Output Feedback")

switch a

    case 1

        n = 21;
        i = [3 4 5 6 7 7 7 7 8 9 10 11 12 12 13 13 14 14 15 16 17 17 18 19 20 20 20 21 11];
        j = [1 1 1 2 3 4 5 6 7 7 7 7 8 9 10 11 12 13 13 14 14 15 15 15 16 17 18 19 5];
        v = ones(size(i));
        A = sparse(i, j, v, n, n); % define adjacency (can be sparse or full)
        G = digraph(A'); % define graph (can be weighted)

        figure;
        plot(G,'b');
        title('Open-loop network')

        D = [2,1]; T = [20,21]; % assign disturbances and targets

    case 2

        n = 19;
        A = zeros(n,n);
        A(3,1) = 1; A(3,2) = 1;
        A(4,3) = 1; A(5,3) = 1;
        A(6,4) = 1; A(6,5) = 1;
        A(7,6) = 1; A(8,6) = 1;
        A(9,6) = 1; A(10,7) = 1;
        A(11,8) = 1; A(11,9) = 1;
        A(12,10) = 1; A(12,11) = 1;
        A(13,12) = 1; A(14,12) = 1;
        A(6,15) = 1; A(15,3) = 1;
        A(16,6) = 1; A(12,16) = 1;
        A(17,6) = 1; A(12,17) = 1;
        A(18,12) = 1; A(3,19) = 1;
        A = magic(n).*A; % add weights
        G = digraph(A'); % define graph

        figure;
        plot(G,'b');
        title('Open-loop network')

        D = [1 2 19]; T = [13 14 18]; % assign disturbances and targets

    case 3

        n = 10;
        A = zeros(n,n);
        A(4,1) = 1; A(4,2) = 1;
        A(5,3) = 1; A(3,7) = 1;
        A(6,4) = 1; A(6,5) = 1;
        A(7,4) = 1; A(7,5) = 1;
        A(8,6) = 1; A(9,7) = 1; A(10,7) = 1;
        A = magic(n).*A; % add weights
        G = digraph(A'); % define graph

        figure;
        plot(G,'b');
        title('Open-loop network')

        D = [1 2 3]; T = [8 9 10]; % assign disturbances and targets

end

%%% Minimal input cardinality DDPSF over networks

V_in = min_card_ddpsf(G, D, T, 'V_in'); % min cardinality input nodes allocation to DDPSF over networks
V_in_ext = min_card_ddpsf_ext(G, D, T, V_in, 'V_in', 'extensive'); % more optimal solutions
V_out = min_card_ddpsf(G, D, T, 'V_out'); % min cardinality output nodes allocation to DDPSF over networks
V_out_ext = min_card_ddpsf_ext(G, D, T, V_out, 'V_out', 'extensive'); % more optimal solutions

% Allocate output measurements and control inputs

V_in1 = V_in_ext{1};
V_out1 = V_in_ext{3};
V_in2 = V_in_ext{3};
V_out2 = V_in_ext{1};

% Check solvability of DDPSF and DDPDF via graphical condition

[Zmax1,Smin1,d1] = ddp_condition(A, D, T, V_in1, 0, V_out1); % d1 = [1 0]
% nodal DDPSF solvable via V_in = V_in1, nodal DDPDF not solvable via V_in = V_in1 and V_out = V_out1
[~,~,d2] = ddp_condition(A, D, T, V_in2, 0, V_out2); % d2 = [1 1]
% DDPSF solvable via V_in = V_in2 and DDPDF solvable via V_in = V_in2 and V_out = V_out2 i.e., (Zmax2,Smin2) is a (C,A,B)-pair of nodes

% Compute state feedback B'*A (F_q = 0)

[Fn, ~, An, Bn, ~, ~, ~, TT] = static_friend(A, V_in1, Zmax1, D, T, 'V_in');
F = Fn*inv(TT);

figure;
plot(digraph((TT*(An-Bn*Fn)*inv(TT))'),'b');
title('Closed-loop network, state feedback')

% Compute output injection A*C' (H_q = 0)

[~, Hn, An, ~, Cn, ~, ~, TT2] = static_friend(A, V_in1, Zmax1, D, T, 'V_out');
H = TT2*Hn;

figure;
plot(digraph((TT2*(An-Hn*Cn)*inv(TT2))'),'b');
title('Closed-loop network, output injection')


%%% Minimal input/output DDPOF over networks

[V_in_optimal, V_out_optimal] = min_card_ddpof(G,D,T); % optimal DDPOF solution w.r.t. |B|+|C|
V_in_optimal = {V_in_optimal}; V_out_optimal = {V_out_optimal}; 

for x = 1:numel(V_in_optimal)
    figure;
    p = plot(G,'b');
    title('Minimal input/output DDPOF over networks')
    nodeColors = 1 * ones(n, 1);
    nodeColors(V_in_optimal{x}) = 3;
    nodeColors(V_out_optimal{x}) = 2;
    p.NodeCData = nodeColors;
    colormap(jet);
    p.MarkerSize = 8;
    hold on;
    legendEntries = {'Input nodes','Output nodes','Other nodes'};
    hRed = scatter(nan, nan, 100, 'r', 'filled');
    hGreen = scatter(nan, nan, 100, 'g', 'filled');
    hBlue = scatter(nan, nan, 100, 'b', 'filled');
    legend([hRed, hGreen, hBlue], legendEntries, 'Location', 'best');
    hold off;
end

% Minimal edge-cut DDPOF

[V_in_opt_e, V_out_opt_e] = min_edge_cut_ddpof(G,D,T); % optimal DDPOF solution w.r.t. edges
V_in_opt_e = {V_in_opt_e}; V_out_opt_e = {V_out_opt_e};

for x = 1:numel(V_in_opt_e)
    figure;
    p = plot(G,'b');
    title('Minimal edge-cut DDPOF over networks')
    nodeColors = 1 * ones(n, 1);
    nodeColors(V_in_opt_e{x}) = 3;
    nodeColors(V_out_opt_e{x}) = 2;
    p.NodeCData = nodeColors;
    colormap(jet);
    p.MarkerSize = 8;
    hold on;
    legendEntries = {'Input nodes','Output nodes','Other nodes'};
    hRed = scatter(nan, nan, 100, 'r', 'filled');
    hGreen = scatter(nan, nan, 100, 'g', 'filled');
    hBlue = scatter(nan, nan, 100, 'b', 'filled');
    legend([hRed, hGreen, hBlue], legendEntries, 'Location', 'best');
    hold off;
end

% Other heuristics for the solutions to the minimal input/output DDPOF over networks

[V_in_opt, V_out, C1] = constrained_optimal_solution(G, D, T, V_out_ext, 1, 'V_in');
[V_in, V_out_opt, C2] = constrained_optimal_solution(G, D, T, V_in_ext, 1, 'V_out');
[V_in_best, V_out_best, ~, ~] = global_constrained_optimal_solution(V_in_opt, V_out, C1, V_in, V_out_opt, C2, 1); % suboptimal solution

if ~isempty(V_in_best)

    for x = 1:numel(V_in_best)
        figure;
        p = plot(G,'b');
        title('$\mathcal{G}$')
        nodeColors = 1 * ones(n, 1);
        nodeColors(V_in_best{x}) = 3;
        nodeColors(V_out_best{x}) = 2;
        p.NodeCData = nodeColors;
        colormap(jet);
        p.MarkerSize = 8;
        hold on;
        legendEntries = {'Input nodes','Output nodes','Other nodes'};
        hRed = scatter(nan, nan, 100, 'r', 'filled');
        hGreen = scatter(nan, nan, 100, 'g', 'filled');
        hBlue = scatter(nan, nan, 100, 'b', 'filled');
        legend([hRed, hGreen, hBlue], legendEntries, 'Location', 'best');
        hold off;
    end

end

[Zmax_OF, Smin_OF, ~] = ddp_condition(A, D, T, V_in_best{1}, 0, V_out_best{1}); % compute S_min(V_out_best{1}) and Z_max(V_in_best{1})
[~, COMr] = compensator_ddp(G, D, T, V_out_best{1}, V_in_best{1}, Smin_OF, Zmax_OF); % reduced order compensator from (C,A,B)-pair (S_min(V_out_best{1}),Z_max(V_in_best{1}))

Gg = COMr.GAMMA.N; % static output feedback
B = eye(n);
B = B(:, V_in_best{1});
C = eye(n);
C = C(V_out_best{1}, :);

figure;
plot(digraph((A-B*Gg*C)'),'b');
title('Closed-loop network, output feedback')

%%% Minimal input/output DDPDF via (C,A,B)-pairs of node sets

[V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair(G, D, T, 1); % Compute (C,A,B)-pairs of node sets
k = 7; % select an output set
V_in_OF = place_controls(A, V_out{k}, T, S_Min{k}); % Input nodes allocation for DDPOF via V_in_OF and V_out{k}
k = 6; % select an input set
V_out_OF = place_sensors(A, V_in{k}, D, Z_Max{k}); % Output nodes allocation for DDPOF via V_in{k} and V_out_OF

k = 9; % select a (C,A,B)-pair
[~, COMr] = compensator_ddp(G, D, T, V_out{k}, V_in{k}, S_Min{k}, Z_Max{k}); % dynamical feedback synthesis
% (order = |Z_Max - S_Min|, F_q = H_q = K_q = 0)

plot_block_colored_digraph(COMr.A_e, n);

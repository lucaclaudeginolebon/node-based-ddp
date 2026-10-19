%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           Compute V_in - V_out           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computing V_in or V_out for DDP by exploiting the mincut/maxflow algorithm,
% This function minimizes the edges to be cut from D to T (it does not guarantee that the solutions minimizes card(V_in)+card(V_out)).
% The function returns when possible a non-degenerate optimal solution.
% INPUTS: G digraph, D disturbance nodes set, T target nodes set
% OUTPUTS: V_in - V_out is an output feedback solution s.t. Im(D) \subseteq W(V_in,V_out) \subseteq Ker(T), W(V_in,V_out) controlled and conditioned
% invariant set of nodes obtained by finding the min cut from D to T in the original graph G
% EXAMPLE : [V_in, V_out] = min_edge_cut_ddpof(G,D,T);
% ========================================================================
% Lebon Luca Claude Gino, May 2026 (Matlab R2022b)
% ========================================================================

function [V_in, V_out] = min_edge_cut_ddpof(G,D,T)

    [Cost_edges_0, V_in, V_out] = min_edge_cut_ddpof_main(G,D,T);
    D = D(:).'; T = T(:).'; V_in = V_in(:).'; V_out = V_out(:).';
    I_d = intersect(D,V_out);
    I_t = intersect(T,V_in);
    Index_0 = numel(union(I_d,I_t));
    A = adjacency(G);
    mask_IT = ~any(A(I_d, T), 2);
    
    if any(mask_IT)
    
        I_keep = I_d(mask_IT);
        mask_receivers = any(A(I_keep, V_in), 1);
        D_new = [D, V_in(mask_receivers)];
        [Cost_edges, V_in_new, V_out_new] = min_edge_cut_ddpof_main(G,D_new,T);
        Index = numel(union(intersect(D,V_out_new),intersect(T,V_in_new)));
    
        if Cost_edges <= Cost_edges_0 && Index <= Index_0
            V_in = V_in_new;
            V_out = V_out_new;
        end
    
    end

end


function [Cost_edges, V_in, V_out] = min_edge_cut_ddpof_main(G,D,T)

    n = numnodes(G);
    
    if any(mod(D,1) ~= 0) || any(mod(T,1) ~= 0) || any(D < 1) || any(T < 1) || any(D > n) || any(T > n)
        error("Error: Target or Disturbance set contains invalid node indices")
    end
    
    k = numedges(G);
    G.Edges.Weight = ones(k,1);
    n_dist = length(D);
    n_targ = length(T);
    
    if ~isempty(T) && ~isempty(D)
    
        for i = 1:n_dist
            G = addedge(G,n+1,D(i),Inf);
        end
    
        for i = 1:n_targ
            G = addedge(G,T(i),n+2,Inf);
        end
    
        [Cost_edges,~,cs,ct] = maxflow(G,n+1,n+2);
    
        if ~ismember(n+1,cs)
            ck = cs;
            cs = ct;
            ct = ck;
        end
      
        [s,t] = findedge(G);
        mask = ismember(s,cs) & ismember(t,ct);
        css = s(mask);
        ctt = t(mask);
        V_out = unique(css);
        V_in = unique(ctt);
    
    else
        error("Error: Target or Disturbance set is empty")
    end

end
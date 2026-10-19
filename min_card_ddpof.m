%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           Compute V_in - V_out           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computing V_in or V_out for DDP by exploiting the mincut/maxflow algorithm,
% This function minimizes the cost function card(V_in)+card(V_out).
% The function returns when possible a non-degenerate optimal solution.
% INPUTS: G digraph, D disturbance nodes set, T target nodes set
% OUTPUTS: V_in - V_out is an output feedback solution s.t. Im(D) \subseteq W(V_in,V_out) \subseteq Ker(T), W(V_in,V_out) controlled and conditioned
% invariant set of nodes obtained by finding the min cut from D to T in the extended graph G_bar
% EXAMPLE : [V_in, V_out] = min_card_ddpof(G,D,T);
% ========================================================================
% Lebon Luca Claude Gino, June 2026 (Matlab R2022b)
% ========================================================================
% Acknowledgement: the sketch of G_bar results from a discussion between L.C.G. Lebon and J. de Vlas in May 2026
% ========================================================================

function [V_in, V_out] = min_card_ddpof(G,D,T)

    [Cost_nodes_0, V_in, V_out] = min_card_ddpof_main(G,D,T);
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
        [Cost_nodes, V_in_new, V_out_new] = min_card_ddpof_main(G,D_new,T);
        Index = numel(union(intersect(D,V_out_new),intersect(T,V_in_new)));
    
        if Cost_nodes <= Cost_nodes_0 && Index <= Index_0
            V_in = V_in_new;
            V_out = V_out_new;
        end
    
    end

end


function [Cost_nodes, V_in, V_out] = min_card_ddpof_main(G,D,T)

    n = numnodes(G);
    
    if any(mod(D,1) ~= 0) || any(mod(T,1) ~= 0) || any(D < 1) || any(T < 1) || any(D > n) || any(T > n)
        error("Error: Target or Disturbance set contains invalid node indices")
    end
    
    k = numedges(G);
    G.Edges.Weight = Inf*ones(k,1);
    n_dist = length(D);
    n_targ = length(T);
    G0 = G;
    
    if ~isempty(T) && ~isempty(D)
    
        for i = 1:n
            G = addedge(G,i+2*n,i,1);
            G = addedge(G,i,i+n,1);
            [~,nid] = outedges(G0,i);
            for j = 1:length(nid)
                G = rmedge(G,i,nid(j));
                G = addedge(G,i,nid(j)+2*n,Inf);
                G = addedge(G,i+n,nid(j),Inf);
            end
        end
    
        for i = 1:n_dist
            G = addedge(G,3*n+1,D(i),Inf);
        end
    
        for i = 1:n_targ
            G = addedge(G,T(i),3*n+2,Inf);
        end
    
        [Cost_nodes,~,cs,ct] = maxflow(G,3*n+1,3*n+2);
    
        if ~ismember(3*n+1,cs)
            ck = cs;
            cs = ct;
            ct = ck;
        end
    
        s = G.Edges.EndNodes(:,1);
        t = G.Edges.EndNodes(:,2);
        mask = ismember(s,cs) & ismember(t,ct);
        sCut = s(mask);
        tCut = t(mask);
        css = sCut(tCut >= n+1 & tCut <= 2*n);
        ctt = tCut(sCut >= 2*n+1 & sCut <= 3*n);
        V_out = unique(css);
        V_in = unique(ctt);
    
    else
    
    error("Error: Target or Disturbance set is empty")
    
    end

end
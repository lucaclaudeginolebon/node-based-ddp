%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute V_out from V_in (suboptimal) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finds V_out for DDPOF once V_in is fixed (e.g. via min_card_ddpsf.m)
% INPUTS: Adjacency matrix A, Control nodes V_in (computed e.g. via min_card_ddpsf.m), D disturbance nodes set, Z_max max controlled invariant set
% OUTPUTS: V_out output nodes for output feedback s.t. Im(D) \subseteq W \subseteq Ker(T), W controlled and conditioned invariant set of nodes
% EXAMPLE : V_out = place_sensors(A, V_in, D, Z_max);
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================

function V_out = place_sensors(A, V_in, D, Z_max)
    
    G = digraph(A');
    m = length(V_in);
    d = length(D);
    V_out = [];
    
    for i = 1:m
        [sout, ~] = findedge(G,inedges(G,V_in(i))); 
        V_out_tmp = intersect(sout, Z_max);
        sout = setdiff(V_out_tmp,V_out);
        l_tmp = length(sout);
        for j = 1:l_tmp
            for k = 1:d
                if ~isempty(shortestpath(G,D(k),sout(j)))
                    V_out = [V_out; sout(j)];
                else
                end
            end
        end
        V_out = unique(V_out);
    end
    
    V_out = setdiff(V_out, V_in);
    V_out = V_out';

end
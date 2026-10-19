%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute V_in from V_out (suboptimal) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finds V_in for DDPOF once V_out is given (e.g. via min_card_ddpsf.m)
% INPUTS: Adjacency matrix A, Control nodes V_in (computed e.g. via min_card_ddpsf.m), T target nodes set, S_min min conditioned invariant node set
% OUTPUTS: V_in input nodes for feedback s.t.  Im(D) \subseteq W \subseteq Ker(T), W controlled and conditioned invariant set of nodes
% EXAMPLE : V_in = place_controls(A, V_out, T, S_min);
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================

 function V_in = place_controls(A, V_out, T, S_min) 

     G = digraph(A');
     n = size(A,1);
     p = length(V_out);
     t = length(T);
     V_in = [];
    
     for i = 1:p
         [~,s_in] = findedge(G,outedges(G,V_out(i)));  
         V_in_tmp = intersect(s_in, setdiff(1:n,S_min)); 
         s_in = setdiff(V_in_tmp,V_in); 
         l_tmp = length(s_in);
         for j = 1:l_tmp
             for k = 1:t
                 if ~isempty(shortestpath(G,s_in(j),T(k)))
                     V_in = [V_in; s_in(j)]; 
                 else
                 end
             end
         end
         V_in = unique(V_in); 
     end
    
     V_in = setdiff(V_in, V_out);
     V_in = V_in';
     
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   (C,A,B)-pairs construction at minimal cost   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computation of (C,A,B)-pairs of node sets at minimal cost for DDPDF (via compensator_ddp.m)
% INPUTS: G digraph, D disturbance nodes set, T target nodes set, show variable (if 1 display, otherwise not)
% OUTPUTS: V_out, V_in = cell list of (C,A,B)-pair generating input/output node sets, s.t.
% exists F|(A-BF)Z_max\subseteqZ_max and exists H|(A-HC)S_min\subseteqS_min,
% and exists N|(A-BNC)S_min\subseteqZ_max, TotalCost = card(V_in)+card(V_out) array, S_min, Z_max (C,A,B)-pair
% EXAMPLE : [V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair(G, D, T);
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================


function [V_out, V_in, TotalCost, S_Min, Z_Max] = cab_pair(G, D, T, show)
    
    if nargin == 3
        show = 0;
    end
    
    if nargin > 4 || nargin < 3    
        disp("Error: number of function inputs is incorrect, they must be 3 or 4");
        V_out = [];
        V_in = [];
        TotalCost = [];
        return    
    end
    
    Vout = min_card_ddpsf(G,D,T,'V_out');
    A = full(adjacency(G))';
    n = numnodes(G);
    n_dist = length(D);
    n_targ = length(T);
    V_out = [];
    V_in = [];
    TotalCost = [];
    Min_dim = []; % stores intersect(Z,complement(S)) = Z-S = dim reduced order compensator
    S_Min = [];
    Z_Max = [];
    
    if ~isempty(Vout)
    
        V_out_effective_tot = [];
        G_cell_effective_tot = [];
        S_min_effective_tot = [];
        S_min_out_effective_tot = [];
        V_in_effective_tot = [];
        Z_max_effective_tot = [];
        Vin_cost_tot = [];
        [V_out_cell, G_cell] = min_card_ddpsf_ext(G,D,T,Vout,'V_out','extensive'); % Provides V_out combinations for DDP at minimum cost
        Vout = [];
        cond1 = 1;
        L2 = 0;
    
        while cond1 && L2 <= max(length(D),length(T))
    
            xx = zeros(numel(V_out_cell),1);
            L1 = length(V_out_cell{1});
            for i = 1:numel(V_out_cell)
    
                Smin = ddp_condition(A, D, T, V_out_cell{i}, 0, 'S_min');
                S_min{i} = Smin;
                Smin_out = [];
    
                for ii = 1:length(Smin)
    
                    [~,s_in] = findedge(G,outedges(G,Smin(ii)));
                    Smin_out = [Smin_out; s_in];
    
                end
    
                Smin_out = unique(Smin_out)';
                S_min_out{i} = Smin_out;
                xx(i) = 1;
    
            end
        
            V_out_effective = V_out_cell(xx ~= 0);
            G_cell_effective = G_cell(xx ~= 0);
            S_min_effective = S_min(xx ~= 0);
            S_min_out_effective = S_min_out(xx ~= 0);    
            V_out_effective_tot = [V_out_effective_tot, V_out_effective];
            G_cell_effective_tot = [G_cell_effective_tot, G_cell_effective];
            S_min_effective_tot = [S_min_effective_tot, S_min_effective];
            S_min_out_effective_tot = [S_min_out_effective_tot, S_min_out_effective];  
            Vout = [Vout, V_out_cell];
            [V_out_cell, G_cell] = min_card_ddpsf_ext(G,D,T,Vout,'V_out','extensive',G_cell{end});
            L2 = length(V_out_cell{1});
            cond1 = L2-L1 > 0;
    
            if ~cond1
    
                result = ones(1, numel(V_out_cell));
    
                for iii = 1:numel(V_out_cell)
                    for jjj = 1:numel(Vout)
    
                        if length(Vout{jjj}) == length(V_out_cell{iii}) && all(ismember(Vout{jjj}, V_out_cell{iii})) && all(ismember(V_out_cell{iii}, Vout{jjj}))
    
                            result(iii) = 0;
                            break
    
                        end
    
                    end
    
                end
    
                if any(result ~= 0)
    
                    cond1 = 1;
                    nzInd = result ~= 0;
                    V_out_cell = V_out_cell(nzInd);
    
                end
    
            end
    
        end
    
        G = digraph(A); % Transpose: invert arrow direction
    
        for j = 1:numel(V_out_effective_tot)
    
            Smin = S_min_effective_tot{j};
            Smin_out = S_min_out_effective_tot{j};
    
            Vin = min_card_ddpsf(G,T,D,'V_out');
            [V_in_cell, G_cell_Vin] = min_card_ddpsf_ext(G,T,D,Vin,'V_out','extensive');
    
            cond2 = 1;
            V_in_all = [];
            LL2 = 0;
    
            while cond2 && LL2 <= max(length(D),length(T)) %+5 % second condition in order to faster the search
    
                LL1 = length(V_in_cell{1});
                V_in_all = [V_in_all, V_in_cell];
    
                for z = 1:numel(V_in_cell)
    
                    Zmax = ddp_condition(A, D, T, V_in_cell{z}, 0, 'Z_max');
                    isSubset1 = all(ismember(Smin, Zmax));
                    isSubset2 = all(ismember(Smin_out, Zmax));
                    isSubset2 = 1;
    
                    if isSubset1 && isSubset2
    
                        %%% (C,A,B) pairs
                        V_in_effective = V_in_cell(z);
                        V_out = [V_out; V_out_effective_tot(j)];
                        V_in = [V_in;  V_in_effective];
                        S_Min = [S_Min; {Smin}];
                        Z_Max = [Z_Max; {Zmax}];
                        Total_cost = length(V_out_effective_tot{j})+length(V_in_cell{z});
                        TotalCost = [TotalCost; Total_cost];
                        Min_dim = [Min_dim; length(setdiff(Zmax,Smin))];
    
                    end
    
                end
    
                Vin = [Vin, V_in_cell];
                [V_in_cell, G_cell_Vin] =  min_card_ddpsf_ext(G,T,D,Vin,'V_out','extensive',G_cell_Vin{end});
                LL2 = length(V_in_cell{1});
                cond2 = LL2-LL1 > 0;
    
                if ~cond2
    
                    result = ones(1, numel(V_in_cell));
    
                    for iii = 1:numel(V_in_cell)
    
                        for jjj = 1:numel(V_in_all)
    
                            if length(V_in_all{jjj}) == length(V_in_cell{iii}) && all(ismember(V_in_all{jjj}, V_in_cell{iii})) && all(ismember(V_in_cell{iii}, V_in_all{jjj}))
    
                                result(iii) = 0;
                                break
    
                            end
    
                        end
    
                    end
    
                    if any(result ~= 0)
    
                        cond2 = 1;
                        nzInd = result ~= 0;
                        V_in_cell = V_in_cell(nzInd);
    
                    end
    
                end
    
            end
    
        end
    
    else
    
        disp('Error: V_out = [], DDP already solved, consider to redefine T, D, or G');
        return
    
    end
    
    V_out = V_out(:); V_in = V_in(:); TotalCost = TotalCost(:); S_Min = S_Min(:); Z_Max = Z_Max(:);

    if show
        T = table((1:numel(V_out))', V_out, V_in, cellfun(@numel, V_out), cellfun(@numel, V_in),  TotalCost, Min_dim, 'VariableNames', {'(C,A,B)-pair number:', 'V_out_opt', 'V_in_opt', '# of Outputs', '# of Inputs', 'Minimal Cost C','Minimal reduced compensator dimension'});
        disp(T);
    end
    
end
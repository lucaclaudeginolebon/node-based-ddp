%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Best constrained optimal V_in-V_out pair when conditioning for V_in or for V_out %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given the list of input/output pairs V_in-V_out provided by constrained_optimal_solution.m find solution that minimize C = card(V_in)+card(V_out), conditioning over V_in or over V_out
% INPUTS: [V_in_opt cell, V_out cell, C1 cost] once constrained w.r.t. V_in, and [V_in cell, V_out_opt cell, C2 cost] once constrain w.r.t. V_out (the order is important), show variable (if 1 display, otherwise not)
% OUTPUTS: V_in cell, V_out cell, C cost, S additional text cell telling towards which among V_in or V_out I have conditioned (V_in and V_out are cell arrays since the optima may not be unique)
% e.g. S = {'V_in', 'V_in', 'V_out'} meaning that I have three V_in-V_out pairs, the first two obtained when conditioning w.r.t. V_in, the last when conditioning w.r.t. V_out
% EXAMPLE: [V_in_opt, V_out, C1] = constrained_optimal_solution(G,D,T,Vout,'V_in'); [V_in, V_out_opt, C2] = constrained_optimal_solution(G,D,T,Vout,'V_out');
% [V_in_best, V_out_best, C, S] = global_constrained_optimal_solution(V_in_opt, V_out, C1, V_in, V_out_opt, C2)
% Warning: cell arrays V_in_opt, V_out and V_in, V_out_opt must be of the same size
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================

function [V_in_best, V_out_best, C, S] = global_constrained_optimal_solution(V_in_opt, V_out, C1, V_in, V_out_opt, C2, show)

    if nargin == 6
        show = 0;
    end

    if ~isempty(C1) && ~isempty(C2)
    
        if C1 < C2
            V_in_best = V_in_opt;
            V_out_best = V_out;
            C = C1;
            Ss = ones(numel(V_in_opt),1);
            S = cell(size(Ss)); 
            S(Ss == 1) = {'V_in'}; 
        elseif C2 < C1
            V_in_best = V_in;
            V_out_best = V_out_opt;
            C = C2;
            Ss = ones(numel(V_in),1);
            S = cell(size(Ss)); 
            S(Ss == 1) = {'V_out'}; 
        elseif C1 == C2
            V_in_best = [V_in_opt; V_in];
            V_out_best = [V_out; V_out_opt];
            C = C1;
            x = cell(1,numel(V_in_best));
    
            for i = 1:numel(V_in_best)
                a = sort(V_in_best{i});
                b = sort(V_out_best{i});
                c = [a b];
                x{i} = c;
            end
    
            V_str = cellfun(@mat2str, x, 'UniformOutput', false);
            [~, idx] = unique(V_str, 'stable');
            allIdx = 1:numel(x);
            repeatedIdx = setdiff(allIdx, idx);
            [~, ~, ic] = unique(V_str, 'stable');
            idxx = find(histcounts(ic, 1:max(ic)+1) > 1);
            repeatedOccurrences = ismember(ic, idxx);
            Ss = [ones(numel(V_in_opt),1); zeros(length(V_in),1)];
            S = cell(size(Ss)); % initialize the cell array
            S(Ss == 1) = {'V_in'}; % fill 'V_in'
            S(Ss == 0) = {'V_out'}; % fill 'V_out'
            S(repeatedOccurrences) = {'V_in and V_out'};
            V_in_best(repeatedIdx) = [];
            V_out_best(repeatedIdx) = [];
            S(repeatedIdx) = [];
        end

        V_in_best = V_in_best(:); V_out_best = V_out_best(:); S = S(:);

        if show
            T = table(V_in_best, V_out_best, cellfun(@numel, V_in_best), cellfun(@numel, V_out_best),  C.*ones(numel(S),1), S, 'VariableNames', {'V_in_best', 'V_out_best','# of Inputs','# of Outputs', 'Minimal Cost C', 'Constraint w.r.t.'});
            disp(T);
        end
    
    else
    
        V_in_best = [];
        V_out_best = [];
        C = [];
        S =[];

        if show
            disp('Error: Target or Disturbance may be empty')
        end
    
    end

end
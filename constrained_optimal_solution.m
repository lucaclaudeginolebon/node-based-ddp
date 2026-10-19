%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimizes V_out resp. V_in conditioned to min(card(V_in)) resp. min(card(V_out)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given the list of optima V_in (resp. V_out), computes (V_in-V_out) pairs and returns the solution V_out (resp. V_in) that gives the minimal card(V_in) (resp. card(V_out))
% INPUTS: G connected directed unweighted graph, D disturbance nodes set, T target nodes set, V (V_in or V_out) nodes list of optima solutions computed via min_card_ddpsf_ext.m
% (extensive constrained optimization search if varargin{2} = 'extensive' in min_card_ddpsf_ext.m), show variable (if 1 display, otherwise not)
% ADDITIONAL INPUTS: varargin{1} = 'V_in' (given opt V_out finds best V_in), varargin{1} = 'V_out' (given opt V_in finds best V_out)
% OUTPUTS: V_in cell, V_out cell, C = cost = card(V_in)+card(V_out) (V_in and V_out are cell arrays since the optima may not be unique)
% EXAMPLE : [V_in, V_out, C] = constrained_optimal_solution(G,D,T,V_out,'V_in'), V_in{i}-V_out{i} is the i-th pair
% Warning: the solution is not globally optimal in general, since we start with the optimal V_in (resp. V_out) and from that we find V_out (resp. V_in) minimizing C 
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================

function [V_in, V_out, C] = constrained_optimal_solution(G, D, T, V, show, varargin)

    if nargin > 5
        if strcmp(varargin{1}, 'V_in')
            a = 1;
            if show
                disp('Search for min V_in, given V_out');
            end
        elseif strcmp(varargin{1}, 'V_out')
            a = 2;
            if show
                disp('Search for min V_out, given V_in');
            end
        else
            disp('Error: additional input typing error, choose ''V_in'' or ''V_out'' as additional input');
            return
        end
    else  
        disp('Error: not enough inputs, choose ''V_in'' or ''V_out'' as additional input');   
        return
    end

    A = adjacency(G)';

    switch a

        case 1

            if ~isempty(V)
                
                m = numel(V);
                C = zeros(m,1);
                Cc = C;
                V_in = cell(1, m);

                for i = 1:m
                    Vout = V{i};
                    S_min = ddp_condition(A, D, T, Vout, 0, 'S_min');
                    Vin = place_controls(A, Vout, T, S_min);
                    V_in{i} = Vin;
                    Cc(i) = length(Vin);
                end

                minValue = min(Cc);
                idx = find(Cc == minValue); 
                C = minValue + length(Vout);
                V_out = V(idx);
                V_in = V_in(idx);
                V_in = V_in(:); V_out = V_out(:);

            else
                if show
                    disp('Error: Target or Disturbance may be empty')
                end
                V_in = [];
                V_out = [];
                C = [];
                return
            end

        case 2

            if ~isempty(V)

                m = numel(V);
                C = zeros(m,1);
                Cc = C;
                V_out = cell(1, m);

                for i = 1:m
                    Vin = V{i};
                    Z_max = ddp_condition(A, D, T, Vin, 0, 'Z_max');
                    Vout = place_sensors(A, Vin, D, Z_max);
                    V_out{i} = Vout;
                    Cc(i) = length(Vout);
                end

                minValue = min(Cc);
                idx = find(Cc == minValue); 
                C = minValue + length(Vin);
                V_in = V(idx);
                V_out = V_out(idx); 
                V_in = V_in(:); V_out = V_out(:);
            else
                if show
                    disp('Error: Target or Disturbance may be empty')
                end
                V_in = [];
                V_out = [];
                C = [];
                return
            end

    end

end
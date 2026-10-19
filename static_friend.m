%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     State feedback and output injection friends    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes F = B'*A s.t. Z_max is (A-BF)-invariant in closed-loop (F_q = 0), and H = A*C' s.t. S_min is (A-HC)-invariant in closed-loop (H_q = 0)
% INPUTS: adjacency matrix A, control nodes V_in (computed e.g. by min_card_ddpsf.m), Z_max max control invariant node set (computed by ddp_condition.m), D disturbance node set, T target node set
% ADDITIONAL INPUT: If varargin{1} = 'V_in' compute Fn s.t. Z_max is An-Bn*Fn invariant via state feedback friend Fn;
% if varargin{1} = 'V_out' compute Hn s.t. S_min is An-Hn*Cn invariant via output injection friend Hn; if empty compute Fn (default)
% OUTPUTS: state feedback Fn, output injection Hn, An, Bn, Cn, Tn, Dn, Cn new state-space matrices after change of basis, TT change of basis matrix
% EXAMPLE : V_in = min_card_ddpsf(G,D,T,'V_in'); Z_max = ddp_condition(A, D, T, V_in, 0,'Z_max'); [Fn, ~, An, Bn, ~, Tn, Dn, TT] = static_friend(A, V_in, Z_max, D, T,'V_in');
% Warning: F in the original basis is F = Fn*inv(TT) = B'*A; H in the original basis is H = TT*Hn = A*C'
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================


function [Fn, Hn, An, Bn, Cn, Tn, Dn, TT] = static_friend(A, V, Z, D, T, varargin)

    if nargin > 5
        if strcmp(varargin{1}, 'V_in')
            V_in = V;
            Z_max = Z;
            a = 1;
            Hn = [];
            Cn = [];
        elseif strcmp(varargin{1}, 'V_out')
            V_out = V;
            S_min = Z;
            a = 2;
            Fn = [];
            Bn = [];
        else
            disp('Error: specify additional input: "V_in" for state feedback, "V_out" for output injection.')
            return
        end
    else
        V_in = V; % default
        Z_max = Z;
        a = 1;
        Hn = [];
        Cn = [];
    end
    
    switch a
    
        case 1
    
            n = size(A,1);
            s = length(Z_max);
            m = length(V_in);
            d = length(D);
            t = length(T);
            K = unique([Z_max V_in],'stable');
            c = length(K);
            TT = zeros(n,n);
    
            for i = 1:c
                TT(K(i),i) = 1;
            end
    
            P = setdiff(1:1:n,K);
    
            for i = 1:(n-c)
                TT(P(i),i+c) = 1;
            end
    
            Bi = zeros(n,m);
    
            for i = 1:m
                Bi(V_in(i),i) = 1;
            end
    
            Di = zeros(n,d);
    
            for i = 1:d
                Di(D(i),i) = 1;
            end
    
            Ti = zeros(t,n);
    
            for i = 1:t
                Ti(i,T(i)) = 1;
            end
    
            An = inv(TT)*A*TT;
            Bn = inv(TT)*Bi;
            Dn = inv(TT)*Di;
            Tn = Ti*TT;
            V = [eye(s,s);zeros(n-s,s)];
            Kn = [V Bn];
            AV = An*V;
            X = Kn'*AV;
            U = X(end-(m-1):end,:);
            Fn = U*V';
    
        case 2
    
            n = size(A,1);
            d = length(D);
            t = length(T);
            y = length(V_out);
            K = unique(S_min,'stable');
            c = length(K);
            TT = zeros(n,n);
    
            for i = 1:c
                TT(K(i),i) = 1;
            end
    
            P = setdiff(1:1:n,K);
    
            for i = 1:(n-c)
                TT(P(i),i+c) = 1;
            end
    
            Di = zeros(n,d);
    
            for i = 1:d
                Di(D(i),i) = 1;
            end
    
            Ti = zeros(t,n);
    
            for i = 1:t
                Ti(i,T(i)) = 1;
            end
    
            Ci = zeros(y,n);
    
            for i = 1:y
                Ci(i,V_out(i)) = 1;
            end
    
            An = inv(TT)*A*TT;
            Dn = inv(TT)*Di;
            Tn = Ti*TT;
            Cn = Ci*TT;
    
    
            [~, col] = find(Cn);
            Hn = zeros(n,y);
            G_new = digraph(An');
            children = [];
    
            for i = 1:length(col)
                children = [children; successors(G_new, col(i))];
            end
    
            children = unique(children);
            m = length(children);
    
            for i = 1:m
                for j = 1:y
                    [~, col_j] = find(Cn(j,:));
                    Hn(children(i),j) = An(children(i),col_j);
                end
            end
    
    end

end
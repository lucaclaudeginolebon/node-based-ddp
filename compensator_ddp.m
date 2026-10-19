%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compensator synthesis using (C,A,B)-pairs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Full/reduced-order compensator synthesis via minimal input-ouput (C,A,B)-pairs of node sets: s.t. d/dt[x;w] = A_e*[x;w]+D_e*d,
% t = T_e*[x;w] and T_e*exp(A_e*t)*D_e = 0, forall t, or equiv.T_e*inv(sI-A_e)*D_e = 0, forall s. A_e = [A+BNC BM; LC K], D_e = [D;0], T_e = [T 0]
% If exists F|(A-BF)Z_max\subseteqZ_max and exists H|(A-HC)S_min\subseteqS_min and exists N|(A-BNC)S_min\subseteqZ_max,
% then the full state compensator (K,L,M,N) exists, with K = A-BF-HC+BNC, L = H-BN, M = NC-F (F_q, H_q, K_q = 0)
% INPUTS: G, D, T, V_out, V_in, S_Min, Z_Max (V_out, V_in, S_Min, Z_Max can be computed via cab_pair.m)
% OUTPUTS: COM, COMr: COM = full-state compensator structure COM.A_e, COM.D_e, COM.T_e, COM.V_in,
% COM.V_out, COM.GAMMA (GAMMA = [K L; M N] compensator with F_q, H_q, K_q = 0) for the (C,A,B)-pair specified by V_out/V_in computed by cab_pair.m
% COMr = minimal-reduced order compensator structure COMr.A_e, COMr.D_e, COMr.T_e, COMr.GAMMA (gammar = [Kr Lr; Mr N] compensator)
% EXAMPLE : [COM, COMr] = compensator_ddp(G, D, T, V_out, V_in, S_Min, Z_Max)
% Warning: The algorithm needs compatible (C,A,B)-pairs of node sets (consider precompute via cab_pair.m)
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================


function [COM, COMr] = compensator_ddp(G, D, T, V_out, V_in, S_Min, Z_Max)

    COM = struct(); COMr = struct();
    
    if ~isempty(V_out)
    
        A = adjacency(G, 'weighted')';
        [Fn, ~, ~, Bn, ~, Tn, Dn, TT1] = static_friend(A, V_in, Z_Max, D, T);
        Bo = TT1*Bn;
        To = Tn*inv(TT1);
        Do = TT1*Dn;
        Fo = Fn*inv(TT1);
        [~, Hn, ~, ~, Cn, ~, ~, TT2] = static_friend(A, V_out, S_Min, D, T, 'V_out');
        N = zeros(length(V_in), length(V_out));
    
        for i = 1:length(V_in)
            for j = 1:length(V_out)
                if ~isempty(A(V_in(i),V_out(j)))
                    N(i,j) = A(V_in(i),V_out(j));
                end
            end
        end
    
        Co = Cn*inv(TT2);
        Ho = TT2*Hn;
        K = A-Bo*Fo-Ho*Co+Bo*N*Co;
        L = -Bo*N+Ho;
        M = -Fo+N*Co;
        gamma.K = K;
        gamma.L = L;
        gamma.M = M;
        gamma.N = N;
        Ae = [A-Bo*N*Co Bo*M; L*Co K];
        De = [Do; zeros(size(Do,1),size(Do,2))];
        Te = [To zeros(size(To,1),size(To,2))];
        COM.A_e = Ae; COM.D_e = De; COM.T_e = Te; COM.V_in = V_in; COM.V_out = V_out; COM.GAMMA = gamma;
        Q = eye(size(A,1));
        x_hat_set = setdiff(Z_Max,S_Min); % min order
        Q = Q(x_hat_set, :);
        r = numel(x_hat_set);
        Kr = Q*K*Q';
        Lr = Q*L;
        Mr = M*Q';
        gammar.Kr = Kr;
        gammar.Lr = Lr;
        gammar.Mr = Mr;
        gammar.N = N;
        Aer = [A-Bo*N*Co Bo*Mr; Lr*Co Kr];
        Der = [Do; zeros(r,size(Do,2))];
        Ter = [To zeros(size(To,1),r)];
        COMr.A_e = Aer; COMr.D_e = Der; COMr.T_e = Ter; COMr.GAMMA = gammar;
    
    else
        disp('Error: V_out = [], DDP already solved, consider to redefine T, D, or G')
        return
    end

end
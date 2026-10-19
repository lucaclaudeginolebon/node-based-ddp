%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Compute multiple min cardinality solutions V_in or V_out  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computing an extensive set of min cardinality solutions V_in or V_out for DDP exploiting the mincut/maxflow algorithm
% INPUTS: G digraph, D disturbance nodes set, T target nodes set, V_in/V_out nodes optimal solution computed by min_card_ddpsf.m
% ADDITIONAL INPUTS: varargin{1} = 'V_in' (finds minimal sets of control nodes), varargin{1} = 'V_out' (finds minimal sets of measurement nodes), default = 'V_in',
% varargin{2} = 'extensive' or 'partial' (finds the sets of solutions, 'extensive' is a (rich albeit not complete)
% combinatorial search of mincut solutions by recursively modifying the weights of G_extended , while 'partial' has a computational cost of approx. M*I*N,
% where M is the length of the solution, I is the maxflow-mincut polynomial time solution cost, and N is the number of optimal solutions), default = 'partial',
% varargin{3} = G_extended: the algorithm searches all solutions equi-cost by starting with the preset weighted graph G_extended (e.g. G_extended =
% G_cell{end}, and V_in_out = V{end} outputs of a previous run of min_card_ddpsf_ext.m)
% OUTPUTS: V_in/V_out min cardinality solutions, G_cell cell array of all weighted graphs for all optimal solutions
% EXAMPLE 1: V_in_ext = min_card_ddpsf_ext(G,D,T,V_in,'V_in','partial')
% EXAMPLE 2: [V_out_cell, G_cell] = min_card_ddpsf_ext(G,D,T,V_in,'V_in','extensive'); [V_out_cell2, G_cell2] = min_card_ddpsf_ext(G,D,T,V_out_cell,'V_out','extensive',G_cell{end});
% Warning: when the function has 7 inputs (i.e., varargin{3} = G_extended, then also the 4th input V_in_out must be a cell array compatible with the cell G_extended)
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================


function varargout = min_card_ddpsf_ext(G,D,T,V_in_out,varargin)

    if nargin > 4
        if strcmp(varargin{1}, 'V_in')
            a = 1; 
        elseif strcmp(varargin{1}, 'V_out')
            a = 2;
        else
            disp('Error: additional input typing error, choose ''V_in'' or ''V_out'' as first additional input and ''extensive'' or ''partial'' as second additional input');
            return
        end
        if nargin > 5
            if strcmp(varargin{2}, 'extensive')
                b = 2;
            elseif strcmp(varargin{2}, 'partial')
                b = 1;
            else
                disp('Error: additional input typing error, choose ''V_in'' or ''V_out'' as first additional input and ''extensive'' or ''partial'' as second additional input');
                return
            end
        else
            b = 1;
        end
        if nargin == 7
            a = 3; % use preset G_extended structure
            G_extended = varargin{3};
            b = 3; % search for suboptimal solutions starting by G_extended, maxflow value is increased
%            disp('Warning: the solutions are found by starting with a preset weighted G_extended graph. Such solutions may not be optimal in the sense of the cardinality')
        elseif nargin > 7
            disp('Error: there are too many inputs. Max nargin = 7')
            return
        end
    elseif nargin == 4
        a = 1; % set 'V_in' as default if varargin is empty
        b = 1; % set 'partial' as default if varargin is empty
    else
        disp('Error: not enough function inputs, choose ''V_in'' or ''V_out'' as first additional input and ''extensive'' or ''partial'' as second additional input');
        return
    end

    n = numnodes(G);
    k = numedges(G);
    G.Edges.Weight = Inf*ones(k,1);
    n_dist = length(D);
    n_targ = length(T);
    
    switch a
    
        case 1
    
            if ~isempty(T) && ~isempty(D)
    
                for i = 1:n
                    [~,nid] = outedges(G,i);
                    for j = 1:length(nid)
                        G = rmedge(G,i,nid(j));
                        G = addedge(G,i+n,nid(j),Inf);
                    end
                    if isempty(intersect(i,D))
                        G = addedge(G,i,i+n,1);
                    else
                        G = addedge(G,i,i+n,Inf);
                    end
                end
    
                for i = 1:n_dist
                    G = addedge(G,2*n+1,D(i),Inf);
                end
    
                for i = 1:n_targ
                    G = addedge(G,T(i)+n,2*n+2,Inf);
                end
    
            else
                V = {};
                disp('Error: Target or Disturbance set is empty')
                return
            end
    
        case 2
    
            if ~isempty(T) && ~isempty(D)
    
                for i = 1:n
                    [~,nid] = outedges(G,i);
                    for j = 1:length(nid)
                        G = rmedge(G,i,nid(j));
                        G = addedge(G,i+n,nid(j),Inf);
                    end
                    if isempty(intersect(i,T))
                        G = addedge(G,i,i+n,1);
                    else
                        G = addedge(G,i,i+n,Inf);
                    end
                end
    
                for i = 1:n_dist
                    G = addedge(G,2*n+1,D(i),Inf);
                end
    
                for i = 1:n_targ
                    G = addedge(G,T(i)+n,2*n+2,Inf);
                end
    
            else
                V = {};
                disp('Error: Target or Disturbance set is empty')
                return
            end
    
        case 3
    
            G = G_extended;
    
    end
    
    switch b
    
        case 1
    
            maxFlow = length(V_in_out);
            maxFlow_n = maxFlow;
            ii = 1;
            V{ii} = V_in_out;
            mu = 0.5; % parameter < 1 strictly
            p = mu/maxFlow;
    
            while maxFlow_n == maxFlow
    
                for i = 1:maxFlow
                    G.Edges.Weight(V_in_out(i)) = 1+p;
                end
    
                [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
    
                if ~ismember(2*n+1,cs)
                    ck = cs;
                    cs = ct;
                    ct = ck;
                end
    
                css = [];
                ctt = [];
    
                for i = 1:length(cs)
                    for j = 1:length(ct)
                        if findedge(G,cs(i),ct(j)) > 0
                            css = [css, cs(i)];
                            ctt = [ctt, ct(j)];
                        end
                    end
                end
    
                V_in_out_new = unique(css);
                V_in_out = V_in_out_new;
                ii = ii+1;
                V{ii} = V_in_out;
    
            end
    
        case 2
    
            G_cell{1} = G;
            maxFlow = length(V_in_out); 
            maxFlow_n = maxFlow;
            ii = 1;
            kk = ii;
            V_is_increasing = 1;
            V{ii} = V_in_out;
            mu = 0.5; % parameter < 1 strictly
            p = mu/maxFlow;
            MMaxFlow = [];
    
            while V_is_increasing
    
                V_in_out = V{kk};
    
                for s = 1:maxFlow
    
                    G.Edges.Weight(V_in_out(s)) = 1+p;
                    [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                    MMaxFlow = [MMaxFlow; maxFlow_n];
    
                    if maxFlow_n == maxFlow
                        if ~ismember(2*n+1,cs)
                            ck = cs;
                            cs = ct;
                            ct = ck;
                        end
    
                        css = [];
                        ctt = [];
    
                        for i = 1:length(cs)
                            for j = 1:length(ct)
                                if findedge(G,cs(i),ct(j)) > 0
                                    css = [css, cs(i)];
                                    ctt = [ctt, ct(j)];
                                end
                            end
                        end
    
                        V_in_out_new = unique(css);
                        ii = ii+1;
                        V{ii} = V_in_out_new;
                        G_cell{ii} = G;
    
                    end
    
                    G.Edges.Weight(V_in_out(s)) = 1; 
    
                end
    
                for i = 1:maxFlow
                    G.Edges.Weight(V_in_out(i)) = 1+p;
                end
    
                [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                MMaxFlow = [MMaxFlow; maxFlow_n];
    
                if maxFlow_n == maxFlow
                    if ~ismember(2*n+1,cs)
                        ck = cs;
                        cs = ct;
                        ct = ck;
                    end
    
                    css = [];
                    ctt = [];
    
                    for i = 1:length(cs)
                        for j = 1:length(ct)
                            if findedge(G,cs(i),ct(j)) > 0
                                css = [css, cs(i)];
                                ctt = [ctt, ct(j)];
                            end
                        end
                    end
    
                    V_in_out_new = unique(css);
                    ii = ii+1;
                    V{ii} = V_in_out_new;
                    G_cell{ii} = G;
    
                end
    
                V_sor = cellfun(@sort, V, 'UniformOutput', false);
                V_str = cellfun(@mat2str, V_sor, 'UniformOutput', false);
                [~, idx] = unique(V_str, 'stable');
                V = V(idx);
                V = V(~cellfun('isempty', V));
                G_cell = G_cell(idx);
                G_cell = G_cell(~cellfun('isempty', G_cell));
                kk = kk+1;
    
                if kk > numel(V) % no more solutions with maxFlow-cost are found
                    varargout{1} = V;
                    varargout{2} = G_cell;
                    return
                end
    
            end
    
        case 3

            maxFlow_in = length(V_in_out{end});
            mu_preset = 1.1;
            VV = [];

            for i = 1:numel(V_in_out)
                V_tmp = V_in_out{i};
                L = length(V_tmp);
                p_preset = mu_preset/L;
                for j = 1:L
                     if isempty(intersect(V_tmp(j),VV)) % updated weight added only once to each node in previous solutions sets
                        G.Edges.Weight(V_tmp(j)) = G.Edges.Weight(V_tmp(j))+p_preset;
                     end
                end
                VV = [V_tmp, VV];
            end

            [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
            mu_increase = 1; 
            condition = 1;
            tol = 0.09;
            k = 1;
            par = 2;

            while condition
                VV = []; % reset
                if abs(1+maxFlow_in-maxFlow_n) <= tol
                    condition = 0;
                    maxFlow = maxFlow_n;
                else
                    for i = 1:numel(V_in_out)
                        V_tmp = V_in_out{i};
                        L = length(V_tmp);
                        p_increase = mu_increase/L;
                        for j = 1:L
                            if isempty(intersect(V_tmp(j),VV)) % updated weight added only once to each node in previous solutions sets
                                G.Edges.Weight(V_tmp(j)) = G.Edges.Weight(V_tmp(j))+p_increase;
                            end
                        end
                        VV = [V_tmp, VV];
                    end
                    [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                    maxFlow_in = maxFlow_in + 1;
                end
                k = k + 1;
                if k == par*n
%                    disp('Warning: unable to compute new solutions; the solution might already exist')
                    maxFlow = maxFlow_n;
                    break
                end
            end

            if ~ismember(2*n+1,cs)
                ck = cs;
                cs = ct;
                ct = ck;
            end

            css = [];
            ctt = [];

            for i = 1:length(cs)
                for j = 1:length(ct)
                    if findedge(G,cs(i),ct(j)) > 0
                        css = [css, cs(i)];
                        ctt = [ctt, ct(j)];
                    end
                end
            end

            V_in_out_new = unique(css);
            V{1} = V_in_out_new;
            G_cell{1} = G;
            ii = 1;
            kk = 1;
            V_is_increasing = 1;
            mu = 0.5; % parameter < 1 strictly
            p = mu/maxFlow;
            MMaxFlow = [];
    
            while V_is_increasing
    
                V_in_out = V{kk};
    
                for s = 1:length(V_in_out)
    
                    G.Edges.Weight(V_in_out(s)) = G.Edges.Weight(V_in_out(s))+p;
                    [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                    MMaxFlow = [MMaxFlow; maxFlow_n];
    
                    if maxFlow_n == maxFlow
                        if ~ismember(2*n+1,cs)
                            ck = cs;
                            cs = ct;
                            ct = ck;
                        end
    
                        css = [];
                        ctt = [];
    
                        for i = 1:length(cs)
                            for j = 1:length(ct)
                                if findedge(G,cs(i),ct(j)) > 0
                                    css = [css, cs(i)];
                                    ctt = [ctt, ct(j)];
                                end
                            end
                        end
    
                        V_in_out_new = unique(css);
                        ii = ii+1;
                        V{ii} = V_in_out_new; 
                        G_cell{ii} = G;
    
                    end
    
                    G.Edges.Weight(V_in_out(s)) = G.Edges.Weight(V_in_out(s))-p;
    
                end

                for i = 1:length(V_in_out)
                    G.Edges.Weight(V_in_out(i)) = G.Edges.Weight(V_in_out(i))+p;
                end
    
                [maxFlow_n,~,cs,ct] = maxflow(G,2*n+1,2*n+2); % compute new solution
                MMaxFlow = [MMaxFlow; maxFlow_n];
    
                if maxFlow_n == maxFlow
                    if ~ismember(2*n+1,cs)
                        ck = cs;
                        cs = ct;
                        ct = ck;
                    end
    
                    css = [];
                    ctt = [];
    
                    for i = 1:length(cs)
                        for j = 1:length(ct)
                            if findedge(G,cs(i),ct(j)) > 0
                                css = [css, cs(i)];
                                ctt = [ctt, ct(j)];
                            end
                        end
                    end
    
                    V_in_out_new = unique(css);
                    ii = ii+1;
                    V{ii} = V_in_out_new;
                    G_cell{ii} = G;
    
                end
    
                V_sor = cellfun(@sort, V, 'UniformOutput', false);
                V_str = cellfun(@mat2str, V_sor, 'UniformOutput', false);
                [~, idx] = unique(V_str, 'stable');
                V = V(idx);
                V = V(~cellfun('isempty', V));
                G_cell = G_cell(idx);
                G_cell = G_cell(~cellfun('isempty', G_cell));
                kk = kk+1;
    
                if kk > numel(V)
                    varargout{1} = V;
                    varargout{2} = G_cell;
                    return
                end
    
            end
    
    end   
    
    nargoutchk(1, 2);
    
    switch nargout
    
        case 2
    
            if b == 2
                V = V(:); G_cell = G_cell(:);
                varargout{1} = V;
                varargout{2} = G_cell;
            else
                disp('Warning: the G_cell is displayed only if varargin{2} = ''extensive''')
                return;
            end
    
        case 1

            V = V(:);
            varargout{1} = V;
    
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           Compute V_in or V_out           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computing V_in or V_out for DDP by exploiting the mincut/maxflow algorithm
% INPUTS: G digraph, D disturbance nodes set, T target nodes set
% ADDITIONAL INPUT: varargin = 'V_in' (finds minimal control nodes), varargin = 'V_out' (finds minimal measured nodes), default = V_in
% OUTPUTS: V_in / V_out is the min cardinality solution s.t. Im(D)\subseteqZ_max(V_in) / S_min(V_out)\subseteqKer(T)
% EXAMPLE : V_in = min_card_ddpsf(G,D,T,'V_in');
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================


function V_in_out = min_card_ddpsf(G,D,T,varargin)

    if nargin > 3
        if strcmp(varargin{1}, 'V_in')
            a = 1;
        elseif strcmp(varargin{1}, 'V_out')
            a = 2;
        else
            disp('Error: additional input typing error, choose ''V_in'' or ''V_out''');
            return
        end
    else
        a = 1; % set V_in as default if varargin is empty
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
                
                [~,~,cs,ct] = maxflow(G,2*n+1,2*n+2);
                
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
                
                V_in_out = unique(css);
        
            else
                 V_in_out = [];
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
                
                [~,~,cs,ct] = maxflow(G,2*n+1,2*n+2);
                
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
                
                V_in_out = unique(css);
        
            else
                 V_in_out = [];
                 disp('Error: Target or Disturbance set is empty')
                 return
            end

    end    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphical solvability condition for DDPS/DF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sufficient condition for DDPSF: Construct Zmax = max controlled invariant node set in ker(T)
% Sufficient condition for DDPDF: If varargin = V_out: Construct Smin = min conditioned invariant node set containing D
% This function can also be exploited to compute Zmax or Smin only
% INPUTS: Adjacency matrix A, D disturbance nodes set, T target nodes set, Input nodes V_in, show variable (if 1 display, otherwise not)
% ADDITIONAL INPUTS: varargin{1} = Output nodes V_out (DDPDF), varargin{1} = 'Z_max' compute only Zmax, varargin{1} = 'S_min' compute only Smin
% OUTPUTS: Zmax max control invariant node subset, Smin min conditioned invariant node superset (if varargin not empty, otherwise Smin = []),
% d = Boolean variable if V_out is empty otherwise d = [d t] with d,t Boolean variables
% (d = 1 if nodal DDPSF is solvable, 0 otherwise, t = 1 if nodal DDPDF is solvable, 0 otherwise)
% The function can give the three outputs [Zmax, Smin, d] or only one of them Zmax or Smin if requested
% EXAMPLE 1: Zmax = ddp_condition(A, D, T, V_in, 0, 'Z_max') % no show, compute Zmax only;
% EXAMPLE 2: [Zmax, ~, d] = ddp_condition(A, D, T, V_in, 0) % no show, DDPSF (Smin = []);
% EXAMPLE 3: Smin = ddp_condition(A, D, T, V, 0, 'S_min') % no show, compute Smin, V = V_out;
% EXAMPLE 4:[Zmax, Smin, d] = ddp_condition(A, D, T, V_in, 0, V_out);
% ========================================================================
% Lebon Luca Claude Gino, March 2026 (Matlab R2022b) 
% ========================================================================

function varargout = ddp_condition(A, D, T, V, show, varargin)
        
    if nargin > 5
        if strcmp(varargin{1}, 'Z_max')
            V_in = V;
            a = 3;
        elseif strcmp(varargin{1}, 'S_min')
            V_out = V;
            a = 2;
        else
            V_in = V;
            V_out = varargin{1};
            a = 1;
        end
    else
        V_in = V;
        V_out = [];
        a = 1;
    end
    
    switch a
    
        case 1
    
            if isempty(V_out)
                Sast = [];
                disp('V_out is empty')
            end
    
            n = size(A,1); 
            Zo = setdiff(1:n,T);
            Zcurr = Zo;
            flag_done = 0;
    
            while ~flag_done
                Zelim = [];
                for i = 1:length(Zcurr)
                    ind = find(A(:,Zcurr(i)));
                    if ~isempty(ind)
                        flag = ismember(ind,union(Zcurr,V_in)); 
                        if sum(flag) < length(ind) 
                            Zelim = union(Zelim, Zcurr(i));
                        end
                    end
                end
                if isempty(Zelim)
                    flag_done = 1;
                    Zast = Zcurr;
                    if show
                        disp(['Size of Zmax = ',num2str(length(Zcurr))]);
                    end
                else
                    if show
                        disp(['Elimininate ',num2str(length(Zelim)),' out of ',num2str(length(Zcurr))]);
                    end
                    Zcurr = setdiff(Zcurr,Zelim);
                end
            end
    
            % check DDP by state feedback
    
            flag_D_in_Zast = ismember(D,Zast);
    
            if sum(flag_D_in_Zast) == length(D)
                d = 1;
                if show
                    disp('The DDPSF problem is solvable');
                end
            else
                d = 0;
                if show
                    disp(['The solvability condition for DDPSF fails: ', num2str(sum(flag_D_in_Zast)),' out of ',...
                        num2str(length(D)),' are decouplable with these ',num2str(length(V_in)),' controls']);
                end
            end
    
            if ~isempty(V_out)
                So = D;
                Scurr = So;
                flag_done = 0;
                while ~flag_done
                    Sadd = [];
                    Scurr_no_out = setdiff(Scurr,V_out);
                    if ~isempty(Scurr_no_out)
                        for i = 1:length(Scurr_no_out)
                            ind = find(A(:,Scurr_no_out(i))); 
                            flag = ismember(ind,Scurr); 
                            if sum(flag) < length(ind) 
                                Sadd = union(Sadd, setdiff(ind,Scurr));
                            end
                        end
                    end
                    if isempty(Sadd)
                        flag_done = 1;
                        Sast = Scurr;
                        if show
                            disp(['Size of Smin = ',num2str(length(Scurr))])
                        end
                    else
                        if show
                            disp(['Add ',num2str(length(Sadd)),' to Scurr of size ',num2str(length(Scurr))])
                        end
                        Scurr = union(Scurr,Sadd);
                    end
                end
    
                % check DDP by dynamical feedback
    
                flag_Sast_in_Zast = ismember(Sast,Zast);
                if sum(flag_Sast_in_Zast) == length(Sast)
                    t = 1;
                    if show
                        disp('The DDPDF problem is solvable by dynamical feedback')
                    end
                else
                    t = 0;
                    if show
                        disp('The solvability condition for DDPDF fails: ')
                        disp(['dim(Smin)= ',num2str(length(Sast)),', dim(Zmax)= ',num2str(sum(flag_Sast_in_Zast)),...
                            ', with ',num2str(sum(flag_Sast_in_Zast)),' nodes of Smin in Zmax, ',...
                            'with ',num2str(length(V_in)),' controls'])
                    end
                end
                d = [d t];
            end
    
        case 2
    
            if ~isempty(V_out)
                So = D;
                Scurr = So;
                flag_done = 0;
                while ~flag_done
                    Sadd = [];
                    Scurr_no_out = setdiff(Scurr,V_out);
                    if ~isempty(Scurr_no_out)
                        for i = 1:length(Scurr_no_out)
                            ind = find(A(:,Scurr_no_out(i)));
                            flag = ismember(ind,Scurr);
                            if sum(flag) < length(ind)
                                Sadd = union(Sadd, setdiff(ind,Scurr));
                            end
                        end
                    end
                    if isempty(Sadd)
                        flag_done = 1;
                        Sast = Scurr;
                        if show
                            disp(['Size of Smin = ',num2str(length(Scurr))])
                        end
                    else
                        if show
                            disp(['Add ',num2str(length(Sadd)),' to Scurr of size ',num2str(length(Scurr))])
                        end
                        Scurr = union(Scurr,Sadd);
                    end
                end
            else
                disp('Error: V_out is empty')
            end
    
        case 3
    
            if ~isempty(V_in)
                n = size(A,1); 
                Zo = setdiff(1:n,T);
                Zcurr = Zo;
                flag_done = 0;
    
                while ~flag_done
                    Zelim = [];
                    for i = 1:length(Zcurr)
                        ind = find(A(:,Zcurr(i)));
                        if ~isempty(ind)
                            flag = ismember(ind,union(Zcurr,V_in)); 
                            if sum(flag) < length(ind)
                                Zelim = union(Zelim, Zcurr(i));
                            end
                        end
                    end
                    if isempty(Zelim)
                        flag_done = 1;
                        Zast = Zcurr;
                        if show
                            disp(['Size of Zmax = ',num2str(length(Zcurr))]);
                        end
                    else
                        if show
                            disp(['Elimininate ',num2str(length(Zelim)),' out of ',num2str(length(Zcurr))]);
                        end
                        Zcurr = setdiff(Zcurr,Zelim);
                    end
                end
    
            else
                disp('Error: V_in is empty')
            end
    end
    
    nargoutchk(1, 3);
    
    switch nargout 

        case 3

            varargout{1} = Zast;
            varargout{2} = Sast';
            varargout{3} = d;

        case 1

            if nargin > 5
                if strcmp(varargin{1}, 'Z_max')
                    varargout{1} = Zast;
                elseif strcmp(varargin{1}, 'S_min')
                    varargout{1} = Sast';
                else
                    disp('Varargin typing error: choose ''Z_max'' or ''S_min'' as first additional input')
                    return
                end
            else
                disp('Error: additional input missing, assign ''Z_max'' or ''S_min'' in varargin{1}')
                return
            end

    end

end
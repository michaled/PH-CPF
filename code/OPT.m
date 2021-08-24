classdef OPT < handle
    
    properties (Access='public')
    end
    
    methods
    end
    
    methods (Static)
        function [ x ] = eigs_fun( CG, b )
            x = pcg(CG,b,[],500);
%             x = gmres(CG,b,[],[],500);
        end
        
        function [ x ] = fmincon( CE, dim )
            % add a constraint to avoid the trivial solution
            % A * ev <= b
            x0 = -ones(dim,1); A = ones(dim,1)'; b = -1;
            options = optimoptions('fmincon',...
                'Display','iter','GradObj','on','Hessian','lbfgs');
            
            x = fmincon(CE,x0,A,b,[],[],[],[],[],options);
        end
        
        function [ x, output ] = minfunc( CE, x0, opt )
            
            [x,~,~,output] = minFunc(CE,x0,opt);
        end

        function [ x,output ] = ipopt( CE, CG, x0, opt )                        
            auxdata = {} ;
            dim = size(x0,1);

            % Set up the auxiliary data.
            options.auxdata = auxdata ;
            
            % Set up the lower and upper bounds
            options.lb = -Inf*ones(dim,1);
            options.ub = Inf*ones(dim,1);

            % Set the IPOPT options.
            options.ipopt.jac_d_constant   = 'no';
            options.ipopt.hessian_constant = 'no';
            options.ipopt.mu_strategy      = 'adaptive';
            options.ipopt.max_iter         = 1000;
            options.ipopt.tol              = 1e-10;
  
            funcs.objective         = CE;
            funcs.gradient          = CG;
            options.ipopt.hessian_approximation      = 'limited-memory';
            options.ipopt.limited_memory_update_type = 'bfgs' ; 
    
            [x,output] = ipopt_auxdata(x0,funcs,options);
        end
        
        function [ x,output ] = lsqnonlin( CE, x0, opt_in ) 
            opt = optimoptions('lsqnonlin');
            
            opt.Algorithm                   = 'levenberg-marquardt';
%            opt.Algorithm                   = 'trust-region-reflective';
            opt.Display                     = 'off'; % 'iter-detailed';
            opt.ScaleProblem                = 'jacobian';
            opt.MaxIterations               = 2000;
            opt.SpecifyObjectiveGradient    = true;
            opt.FiniteDifferenceType        = 'central';
            opt.CheckGradients              = false;
            
            if nargin > 2                
                opt.FunctionTolerance           = opt_in.funcTol;
            end

            [x,~,~,~,output] = lsqnonlin(CE, x0, [], [], opt);         
        end
        
        function [ x, output ] = direct( CH, x0, opt )
            [B,f,l,u] = opt.directparams();
            
            opt = optimoptions('quadprog','Display','iter',...
                               'HessMult',CH,'TolFun',opt.TolFun,...
                               'Algorithm','trust-region-reflective');
                           
            [x,~,~,output] = quadprog(B,f,[],[],[],[],l,u,...
                                      x0,opt);
        end

        function [ x, output ] = iter( CE, CG, opttool, x0, opt )
            dim = size(x0,1);
            if strcmp(opttool,'lsqnonlin') == 1
                % 
                [x,output] = OPT.lsqnonlin(CE,x0,opt);
            elseif strcmp(opttool,'ipopt') == 1
                % COFMAPS IPOPT: fast, gives the best results
                [x,output] = OPT.ipopt(CE,CG,x0);
            elseif strcmp(opttool,'minfunc') == 1    
                % COFMAPS MINFUNC: fastest, but IPOPOT gives better results
                [x,output] = OPT.minfunc(CE,x0,opt);
            elseif strcmp(opttool,'eigs') == 1    
                % COFMAPS EIGS: slowest, but gives a good approximation
                [x,~] = eigs(@(b) OPT.eigs_fun(CG,b),dim,1,'SM');
            elseif strcmp(opttool,'fmincon') == 1    
                % COFMAPS FMINCON: slow, gives a poor result
                x = OPT.fmincon(CE,dim);
            end
            
            x = x(:,1);
        end
    end
    methods (Static)
        % Note that the epsilon and allowed error are fixed, 
        % and do not depend on the range of f. Check output manually
        % if function fails check
        function ae = dbg_grad_ls(fun, x0, sf, n, varargin)
          delta = 1e-3; 
          se = 0;
          
          fprintf(' Iter       j             err');
          fprintf('           g_est               g               f\n')

          if isempty(sf)
              sf = 1:numel(x0);
          end
          for i=1:n
            T = x0;
            j = randsample(sf,1);
            T0=T; T0(j) = T0(j)-delta;
            T1=T; T1(j) = T1(j)+delta;

            [f,g] = fun(T, varargin{:});
            f0 = fun(T0, varargin{:});
            f1 = fun(T1, varargin{:});

            g_est = (f1-f0) / (2*delta);
            e = norm(g(:,j) - g_est);

            fprintf('% 5d  % 6d % 15g % 15f % 15f % 15f\n', ...
                    i,j,e,full(g(1,j)),g_est(1),f(1));

            se = se + e;
            
            if e > delta || isnan(e)
                [full(g(:,j)), g_est]
                error('?')
            end
          end

          ae = se/n;        
        end
        
        % Generate n random invertible 2x2 matrices whose abs(determinant) is 
        % bigger then mind. Works by rejection
        % returned as a n x 4 vector
        function a = rand_invertible_mats2(n, mind)
            a = rand(n,4);
            s = a(:,1).*a(:,4) - a(:,2).*a(:,3);
            locs = find(abs(s) < mind);
            while ~isempty(locs)
                a(locs,:) = rand(length(locs),4);
                s = a(:,1).*a(:,4) - a(:,2).*a(:,3);
                locs = find(abs(s) < mind);
            end
        end
    end
end
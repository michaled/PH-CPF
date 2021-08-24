classdef CPF < handle

    properties (Access='public')
        p       % parameters
        
        name            % input mesh name 
        M               % input mesh
                
        siz             % which sizing to use for the vector fields, 'i'=id * gradsz, 
                        % 'k' = curvature metric corresponding to siz, 2nf x 2nf
        g               % Should be per face, and symmetric positive definite
        gd, gde         % determinant(per face, and edge) of g, for smoothness weighting
        gw              % weighting for smoothness energy
        mel2d           % 2d average edge length according to M and g          

        SO              % Shape operator, per face, symmetric 2nf x 2nf
        
        n               % n-direction fields
        lsm             % weight for smoothness energy 
        lc              % weight of continuity constraint 
        lsi             % weight of sizing constraint
        la              % weight for alignment constraint
        lin             % weight for injectivity constraint
        lcj             % weight for conjugacy constraint
        lo              % weight for orthogonality constraint
        k_eps           % Epsilon for curvature sizing, as a fraction of 
                        % 1/model's bounding box
        delta           % max distance from surface, as a fraction of 
                        % model's bounding box
        eta             % Area of max element, as a fraction of total 
                        % surface area    
        rhow            % 1, weigh the alignment with RHO squared, 2 - polynomial weighting
        rhowf, rhowe    % rho weighting on faces and edges
        smw             % diagonal matrix weighting for inner edges
        
        iter_opt        % if > 0, iterate the optimization additional 
                        % iter_opt times, while increasing the integrability 
                        % penalty, default 0
        iter_opt_mult   % multiplicity constant for iter_opt, default 2
        cont_n          % n used for continuity 
        
        inj_br_s        % parameter for the injectivity barrier function
        inj_br_eps      % parameter for the injectivity barrier function
        
        orth_obj        % Orthogonality constraint: 1 - unnormalized, 2 - normalized
        
        opt,opttool     % optimization tool used 
        abseps
        function_tol    % functionTolerance for optimization
        optinit         % @see opt_init
        optinitval      % @see opt_init
        optC            % global scale value for energy

        z0,U0,V0        % initial values
        cons2, cons4    % Alignment constraints from SHAPEOP
        w               % Initial u

        align_boundary  % 0 - don't align, n > 0 - align with n
        align_curvature % 0 - don't align, n > 0 - align damin with n
        align_guiding   % 0 - don't align, 
                        % 1 - align with guiding field with connected
                        % components supression of damin
                        % 2 - align according to guiding field (mostly
                        % damin, depending on curvature regions)
                        % 3 - align according to smooth 2-rosy with cons2
                        % and cons4 constraints
                        
        A               % Alignment info:
                        % A.f - faces with alignment constraints, na x 1
                        % A.d - the constraints vectors, 2na x 1
                        % A.n - the constraint degree, scalar or na x 1
                        % A.Ai, A.Bi - arrays for computing the constraints
                        
        remove_boundary_faces_from_alignment
        ell_cons_guiding
                        
        Bf              % nfi x nf, subset of identity matrix, with rows 
                        % only for interior faces

        z               % power vector of pulled back edges, nie x 2
        U,V             % final coordinate vector fields, nf x 3, each
        U_prev, V_prev  % Coordinate vector fields from the previous iteration, nf x 3, each
        Uc, Vc          % final combed coordinate vector fields 
                        % (U,V) -> (gradu, gradv) -> rawfield (depends on n) ->
                        % combed rawfield -> (graduc, gradvc) -> (Uc, Vc)

        gradu, gradv    % final candidate gradients, uncombed, nf x 3 each
        graduc, gradvc  % final candidate gradients, combed, nf x 3 each

        Fc,             % Combed rawfield, nf x 3*n
        mtch, mtchc,    % combed and uncombed matching, ne x 1
        pd              % parameterization data file name (between comb and param)
        
        fvals, grads    % partial values of objectives and gradients at the 
                        % current iteration
        Je              % Rotates by pi/2 a 2*nie vector
                        
        M_c_ni, M_c     % Cut meshes, with non integer (ni) and integer param
        cv2v, S         % map from cut vertices to original vertices, 
                        % singularities
        Mr              % Remeshed mesh
        Mrd             % dual remeshed mesh (if dualize == 1)
        Mp              % Planarized remeshed mesh (if planar == 1)
        
        uv,fuv          % parameterization (output of mesher)
        lr              % parametrization edge length for parameterization and meshing
        auto_lr         % 1 - compute lr automatically
        
        E,ES,EPSI       % energy values during minimization
        EPSI_C, EPSI_S
        EPSI_CJ, EPSI_CO
        EPSI_A
        EPSI_IN
        AR              % Total parameterized area during minimization
        mtime,miter     % computation time and iteration count
        epf             % output solution
        perr            % Poisson error before integers
        x2cornerMat     % convert from poisson output Ax\b to function values 
                        %  at vertices
        
        profile         % Run profiler on optimization
        
        dualize         % Output dual of remeshed mesh
        planar          % Output also planarization of output mesh
        mp              % Max planarity error, default 1
        psnap           % Snap to boundary when planarizing
        
        figs            % open figures while running
        
        cutterdir       % mesh cutter dir and bin
        cutterbin
        
        paramerdir      % parameterizer dir and bin
        paramerbin
 
        mesherdir       % mesher dir and bin
        mesherbin
        
        dualizerdir     % dualizer dir and bin
        dualizerbin

        planarizerdir   % planarizer dir and bin
        planarizerbin
        
        datadir         % meshes dir    
        outdir          % output dir 
        outdir_tmp    % additional files' output dir
        
        planarize_sing_cor
        
    end
    
    methods
        % Constructor 
        function [ cpf ] = CPF( name )
            
            % parameters
            cpf.p = inputParser;
                                   
            % cross fields
            addOptional(cpf.p,'n',6);
            addOptional(cpf.p,'cdv',[]);
            addOptional(cpf.p,'so',.1);
            
            % sizing
            addOptional(cpf.p,'siz','i');
            addOptional(cpf.p,'delta',.1);
            addOptional(cpf.p,'eta',.05);

            % injectivity
            addOptional(cpf.p, 'inj_br_s', 0.1);
            addOptional(cpf.p, 'inj_br_eps', 1e-3);
            
            % orthogonality
            addOptional(cpf.p, 'orth_obj', 1);
            
            % optimization parameters
            addOptional(cpf.p,'lsm',.1);
            addOptional(cpf.p,'lc',.9);
            addOptional(cpf.p,'lsi',.1);            
            addOptional(cpf.p,'la',.1);
			addOptional(cpf.p,'lin',.1);
            addOptional(cpf.p,'lcj',.1);
            addOptional(cpf.p,'lo',.1);
            
            % Entering alignment constraints directly
            addOptional(cpf.p,'A',[]);

            % Specific alignment types
            addOptional(cpf.p,'align_boundary',0);
            addOptional(cpf.p,'align_curvature',0);
            addOptional(cpf.p,'align_guiding',-1);
            addOptional(cpf.p,'rhow',-1);
            addOptional(cpf.p,'remove_boundary_faces_from_alignment',1);
            addOptional(cpf.p,'ell_cons_guiding',1);

            addOptional(cpf.p,'iter_opt',0);
            addOptional(cpf.p,'iter_opt_mult',2);
            
            addOptional(cpf.p,'cont_n',0);           
            
            addOptional(cpf.p,'opt','iter');
            addOptional(cpf.p,'opttool','lsqnonlin');
            addOptional(cpf.p,'abseps',1e-12);
            addOptional(cpf.p,'function_tol',1e-6);
            addOptional(cpf.p,'optinit','r'); 
            addOptional(cpf.p,'optinitval',[0,0,0]); 

            addOptional(cpf.p,'optC',1); 

            % cutter, parameterizer & mesher
            addOptional(cpf.p,'cutterdir','external/MeshCutter/');
            addOptional(cpf.p,'cutterbin','MeshCutter_bin');
            addOptional(cpf.p,'lr',-1);
            addOptional(cpf.p,'auto_lr',1);
            addOptional(cpf.p,'dualize',1);
            addOptional(cpf.p,'planar',1);
            addOptional(cpf.p,'mp',.1);            
            addOptional(cpf.p,'psnap',0);            
            addOptional(cpf.p,'paramerdir','external/Parameterization/');
            addOptional(cpf.p,'paramerbin','Parametrization_bin');
            addOptional(cpf.p,'mesherdir','external/Mesher/');
            addOptional(cpf.p,'mesherbin','Mesher_bin');
            addOptional(cpf.p,'dualizerdir','external/Dualizer/');
            addOptional(cpf.p,'dualizerbin','Dualizer_bin');
            addOptional(cpf.p,'planarizerdir','external/Planarizer/');
            addOptional(cpf.p,'planarizerbin','Planarize_bin');

            % Profiler
            addOptional(cpf.p,'profile',0);
            
            % Figures
            addOptional(cpf.p,'figs',0);
            
            % data directories
            addOptional(cpf.p,'datadir','../data/');
            addOptional(cpf.p,'outdir','results/');
            addOptional(cpf.p,'outdir_tmp','');
            
            addOptional(cpf.p,'planarize_sing_cor',0);
            
            % create MESH data structures
            cpf.name = name; cpf.M = MESH( cpf.name );
            % Check mesh is single component
            
            addOptional(cpf.p,'U0', rand(cpf.M.nf,3)); 
            
            
            sz = sparse(cpf.M.nie, cpf.M.nie);
            so = speye(cpf.M.nie, cpf.M.nie);            
            cpf.Je = [sz,-so;so,-sz];
            
            % interior faces
            cpf.Bf = speye(cpf.M.nf,cpf.M.nf);
        end
        
        % Problem parameters
        function set_param( cpf, varargin )
            
            varargin = varargin{1};
            
            % parse & store input parameters
            parse( cpf.p, varargin{:} );
            
            % figs
            cpf.figs = cpf.p.Results.figs;
            
            % cross fields
            cpf.n = cpf.p.Results.n;
            
            % interior faces
            cpf.remove_boundary_faces_from_alignment  = cpf.p.Results.remove_boundary_faces_from_alignment;
            if cpf.remove_boundary_faces_from_alignment 
                cpf.Bf = cpf.Bf(cpf.M.bf == 0,:);
            end
            
            
            % sizing
            cpf.siz = cpf.p.Results.siz;
            cpf.lr = cpf.p.Results.lr;
            cpf.auto_lr = cpf.p.Results.auto_lr;
            if cpf.lr>0
                cpf.auto_lr = 0;
            end

            % injectivity
            cpf.inj_br_s = cpf.p.Results.inj_br_s;
            cpf.inj_br_eps = cpf.p.Results.inj_br_eps;
            
            % orthogonality
            cpf.orth_obj = cpf.p.Results.orth_obj;
            
            % Scale by model's bounding box diagonal
            cpf.delta = cpf.p.Results.delta * cpf.M.bbd;
            cpf.eta = cpf.p.Results.eta * sum(cpf.M.ta);
            cpf.k_eps = cpf.delta / cpf.eta;
            
            % identity
            if cpf.siz == 'i'
                cpf.g = cpf.eta * speye(2*cpf.M.nf);  
            % Curvature
            elseif cpf.siz == 'k'
                e = cpf.k_eps; r = cpf.delta;
                [cpf.g,S] = SHAPEOP.curvature_sizing(cpf.M, e, r);
            % Hodge sizing
            elseif cpf.siz == 'h'
                [ey,~] = eigs(cpf.M.hodge,[],4,1e-5); 
                ey = ey(:,1);
                s = MESH.normv(ey,2);
                cpf.g = spdiags([s;s].^2,0,2*cpf.M.nf,2*cpf.M.nf);
            else
                error('unsupported sizing parameter');
            end
            
            cpf.SO = SHAPEOP.shapeop(cpf.M);
            

            % compute average 2d edge length according to g and edge
            % lengths
            [cpf.mel2d,el2] = cpf.comp_2dmel;
            cpf.gw = speye(cpf.M.nie, cpf.M.nie);
            
            % optimization weights
            cpf.lsm = cpf.p.Results.lsm;
            cpf.lc = cpf.p.Results.lc;
            cpf.lsi = cpf.p.Results.lsi;
            cpf.la = cpf.p.Results.la;
			cpf.lin = cpf.p.Results.lin;
            cpf.lcj = cpf.p.Results.lcj;
            cpf.lo = cpf.p.Results.lo;

            cpf.iter_opt = cpf.p.Results.iter_opt;
            cpf.iter_opt_mult = cpf.p.Results.iter_opt_mult;            
            cpf.cont_n = cpf.p.Results.cont_n;
            if cpf.cont_n == 0
                cpf.cont_n = cpf.n;
            end

            cpf.rhow = cpf.p.Results.rhow;
            if cpf.rhow==-1
                if cpf.n==4
                    cpf.rhow = 1;
                elseif cpf.n==6
                    cpf.rhow = 2;
                else 
                    error('rhow is undefined.');
                end
            end
            
            % normalize optimization parameters so they have the same units
            % according to n, g and mesh edge lengths            
            % 2 vectors per face, 2d mel^2    
            lsi_units = cpf.M.nf*2;

            % edge length ^ n ^ 2
            lc_units = cpf.mel2d^(2*cpf.n)*cpf.M.ne;
            lsm_units = lc_units;

            lin_units = cpf.M.nf;
            lcj_units = cpf.M.nf*2;
            lo_units = cpf.M.nf*2;          
            la_units = cpf.mel2d^2;
            
            % Make parameters unitless
            cpf.lsi = cpf.lsi / lsi_units;
            cpf.lc = cpf.lc / lc_units;
            cpf.lsm = cpf.lsm / lsm_units;
            cpf.lin = cpf.lin / lin_units;
            cpf.lcj = cpf.lcj / lcj_units;
            cpf.lo = cpf.lo / lo_units;
            cpf.la = cpf.la / la_units;
            
            % Alignment
            cpf.align_boundary = cpf.p.Results.align_boundary;
            cpf.align_curvature = cpf.p.Results.align_curvature;
            cpf.align_guiding = cpf.p.Results.align_guiding;

            if cpf.align_guiding==-1
                if cpf.n==4
                    cpf.align_guiding = 2;
                elseif cpf.n==6
                    cpf.align_guiding = 3;
                else 
                    error('align_guiding is undefined.');
                end
            end
            
            cpf.A.f = []; cpf.A.d = []; cpf.A.ns = [];
            cpf.add_bdry_align();
            cpf.add_curv_align();
            cpf.add_guide_align();

            la_units = la_units*length(cpf.A.f);
            if la_units > 0
                cpf.la = cpf.la / la_units;
            end
            
            % Visualize alignment
            cpf.A.do = cpf.A.d;
            cpf.A.no = cpf.A.ns;
            cpf.show_alignment(cpf.A);
            
            % Initialize alignment structure
            na = length(cpf.A.f); 
            nf = cpf.M.nf;
            nie = cpf.M.nie;
            II = [1:4*na];
            JJ = repmat(cpf.A.f',1,4) + [zeros(1, na), nf*ones(1,na), ...
                2*nf*ones(1,na), 3*nf*ones(1,na)];
            SS = ones(size(II));
            % Ai: Indexes which faces we need, 4na x 4nf
            cpf.A.Ai = sparse(II, JJ, SS, 4*na, 4*cpf.M.nf);
            % Bi: Indexes which outputs we need, na x 2na
            % Only imaginary values, so second half
            II = [1:na]; JJ = II + na; SS = ones(size(II));
%            % For aligning V with prescribed direction use:
%             if cpf.n==6
%                 JJ = II;
%             end
            cpf.A.Bi = sparse(II,JJ,SS, na, 2*na);            
            % weigh with RHO
            if cpf.align_guiding ~= 0 && cpf.rhow > 0
                w = cpf.A.RHO;
                if cpf.rhow == 2
                    % More agressive weighing near 1
                    w = w.^3 - 3*w.^2 + 3*w;
                end
                cpf.rhowf = w;
                cpf.A.Bi = spdiags(w(cpf.A.f),0,na,na) * cpf.A.Bi;
            else
                cpf.rhowf = zeros(nf,1);
                cpf.rhowf(cpf.A.f) = 1;
            end
            cpf.rhowe = mean(cpf.rhowf(abs(cpf.M.e2t(cpf.M.inner_edges,1:2))),2);
            
            % Choose diagonal weights for smoothness
            cons2 = SHAPEOP.guiding_field( cpf.M );            
            locs2 = find(MESH.normv(cons2) > 1e-5);
            ww = zeros(nf,1); 
%            ww(locs2) = 1;
            wwe = mean(ww(abs(cpf.M.e2t(cpf.M.inner_edges,1:2))),2);
            wwe(wwe > 0) = 1;
            cpf.smw = spdiags(1-[wwe;wwe],0,2*nie,2*nie);
            
            % Constraints: 2*na x 1
            % If needed project on basis
            locs2 = [cpf.A.f; cpf.A.f + cpf.M.nf];
            locs3 = [locs2; cpf.A.f + 2*cpf.M.nf];
            if size(cpf.A.d,2) == 3
                cpf.A.d = cpf.M.EB(locs2, locs3)*cpf.A.d(:);
            % If needed, vectorize
            elseif size(cpf.A.d,2) == 2
                cpf.A.d = cpf.A.d(:);
            end
            % Imaginary is 0 --> can be both positive and negative
            assert(sum(mod(cpf.A.ns,2))==0,'alignment n should be even');
            cpf.A.ns = cpf.A.ns / 2;

            % optimization parameters
            cpf.abseps = cpf.p.Results.abseps;
            cpf.function_tol = cpf.p.Results.function_tol;
            cpf.opt = cpf.p.Results.opt;
            cpf.opttool = cpf.p.Results.opttool;
            cpf.optinit = cpf.p.Results.optinit;            
            cpf.optinitval = cpf.p.Results.optinitval;
            cpf.U0 = cpf.p.Results.U0;
            cpf.optC = cpf.p.Results.optC;
            
            % cutter, parameterizer,  mesher, dualizer & planarizer
            cpf.cutterdir = cpf.p.Results.cutterdir;
            cpf.cutterbin = cpf.p.Results.cutterbin;
            cpf.lr = cpf.p.Results.lr;
            cpf.paramerdir = cpf.p.Results.paramerdir;
            cpf.paramerbin = cpf.p.Results.paramerbin;
            cpf.mesherdir = cpf.p.Results.mesherdir;
            cpf.mesherbin = cpf.p.Results.mesherbin;
            cpf.dualize = cpf.p.Results.dualize;
            cpf.dualizerdir = cpf.p.Results.dualizerdir;
            cpf.dualizerbin = cpf.p.Results.dualizerbin;
            cpf.planar = cpf.p.Results.planar;
            cpf.planarizerdir = cpf.p.Results.planarizerdir;
            cpf.planarizerbin = cpf.p.Results.planarizerbin;
            
            % Dualize
            if cpf.n==6
                cpf.dualize = 1;
            end
            
            % Planarity parameter 
            cpf.mp = cpf.p.Results.mp;
            
            % Snap to boundary when planarizing
            cpf.psnap = cpf.p.Results.psnap;

                       
            % Profiler
            cpf.profile = cpf.p.Results.profile;
            
            % data directories
            cpf.datadir = cpf.p.Results.datadir;
            cpf.outdir = cpf.p.Results.outdir;
            cpf.outdir_tmp = cpf.p.Results.outdir_tmp;
            if isempty(cpf.outdir_tmp)
                cpf.outdir_tmp = [cpf.outdir, 'tmp/'];
            end
            
            % 0. set preliminaries: meshes, eigenbases and fmaps
            cpf.prep_prelim();
            fprintf('\n');
            
        end
        
        
        function str = get_param( cpf )
            str = sprintf(...
                ['n%d_siz_%s_lsm%.2g_lc%.2g_lsi%.2g_la%.2g_lin%.2g_lcj'...
                '%.2g_lo%.2g_init_%s_delta_%.2g_rhow_%d_alg_%d'],...
                cpf.n, cpf.siz, cpf.lsm, cpf.lc, cpf.lsi, cpf.la, cpf.lin,...
                cpf.lcj, cpf.lo, cpf.optinit, cpf.delta, cpf.rhow, cpf.align_guiding);
        end
        
        % Adds to A the boundary alignment terms
        function add_bdry_align( cpf )
            CM = cpf.M;
            A = cpf.A;
            
            if cpf.align_boundary == 0
                return;
            end
            n = cpf.align_boundary;
            
            bloc = find(CM.ie == 0);
            tb = abs(CM.e2t(bloc,1));
            X = CM.vertices;
            d = X(CM.edges(bloc,1),:)-X(CM.edges(bloc,2),:);
            A.f = [A.f; tb];
            A.d = [A.d; d];
            A.ns = [A.ns; repmat(n, length(tb), 1)];
            
            cpf.A = A;
        end

        % Adds to A the curvature alignment terms
        function add_curv_align( cpf )
            CM = cpf.M;
            A = cpf.A;
            
            if cpf.align_curvature == 0
                return;
            end
            n = cpf.align_curvature;
            [~, ~, ~, damin, ~] = SHAPEOP.guiding_field2( CM );
            
            tb = [1:CM.nf]';
            d = damin;
            A.f = [A.f; tb];
            A.d = [A.d; d];
            A.ns = [A.ns; repmat(n, length(tb), 1)];
            
            cpf.A = A;
        end
        
        % Adds to A the guiding field alignment terms
        function add_guide_align( cpf )
            CM = cpf.M;
            A = cpf.A;
            
            if cpf.align_guiding == 0
                return;
            end
            
            if cpf.align_guiding == 1 
                [cons2, cons4, ~, damin, dmin, RHO, PHI, locse] = ...
                    SHAPEOP.guiding_field( CM );
            elseif cpf.align_guiding == 2
                [cons2, cons4, ~, damin, dmin, RHO, PHI, locse] = ...
                    SHAPEOP.guiding_field2( CM );
            elseif cpf.align_guiding == 3
                % this should not be used for quad meshes
                assert(cpf.n ~= 4);
                [cons2, cons4, cpf.w, damin, dmin, RHO, PHI, locse] = ...
                    cpf.optimized_guiding_field;
            else
                error('?');
            end
            
            cpf.cons2 = cons2; cpf.cons4 = cons4;
            locs2 = find(MESH.normv(cons2) > 1e-5);
            locs4 = find(MESH.normv(cons4) > 1e-5);
            
            % For n=6, alignment is 6 in elliptic, and 2 everywhere else
            if cpf.n == 6
                locs6 = intersect(locs2, locse);
                locs2ne = setdiff(locs2, locse);
                nel = 6;
                A.ns = [A.ns; ...
                    repmat(2, length(locs2ne), 1); 
                    repmat(4, length(locs4), 1);
                    repmat(nel, length(locs6), 1)];
                tb = [locs2ne; locs4; locs6];
                d = [cons2(locs2ne,:); cons4(locs4,:); cons2(locs6,:)];                
            % For n==4, all alignment should be 4-rosy
            elseif cpf.n == 4
                A.ns = [A.ns; ...
                    repmat(4, length(locs2), 1); repmat(4, length(locs4), 1)];
                tb = [locs2; locs4];
                d = [cons2(locs2,:); cons4(locs4,:)];
            end

            A.f = [A.f; tb];
            A.d = [A.d; d];
            
            A.RHO = RHO;
            A.PHI = PHI;
            A.dmin = dmin;
            A.damin = damin;
            
            cpf.A = A;
        end
        
        % Pre-process several computations
        function prep_prelim( cpf )
            % For accumulating energy values
            cpf.E = []; cpf.ES = []; cpf.EPSI = []; 
            cpf.EPSI_C = []; cpf.EPSI_S = []; cpf.EPSI_CJ = []; 
            cpf.EPSI_A = [];
			cpf.EPSI_IN = [];
            cpf.AR = [];
        end
        
        % Compute 2d edge lengths according to the metric g
        function [mel2, el2] = comp_2dmel(cpf)
            CM = cpf.M; Cg = cpf.g;
            
            [E1, E2, E3] = CM.face_edges;
            E1i = CM.EB*E1(:); E2i = CM.EB*E2(:); E3i = CM.EB*E3(:);
            E1l2 = CM.NORMV(E1i, 2, Cg);
            E2l2 = CM.NORMV(E2i, 2, Cg);
            E3l2 = CM.NORMV(E3i, 2, Cg);
            
            mel2 = mean([E1l2; E2l2; E3l2]);            
            el2 = [E1l2, E2l2, E3l2];
        end
        
                
        % Smoothness, Injectivity, Continuity, Alignment, Conjugacy
        % z: 2*nie, u: 2*nf, v: 2*nf
        % ey = [real(z);imag(z);u;v];
        function [ e, g ] = comp_energy_ls( cpf, ey, ~ )
            
            CM = cpf.M; 
            nf = CM.nf;
            nie = CM.nie;
            
            z = ey(1:nie) + 1i*ey(nie+1:2*nie);
            U = ey(2*nie+1:2*nie+2*nf); V = ey(2*nie+2*nf+1:end);
            
            [es, gs_U, gs_V, gs_z] = smoothness_ls(cpf, U, V);
            gs = [gs_z, gs_U, gs_V];
            
            [epsi, gpsi_UV, gpsi_z] = constraints_ls(cpf, U, V, z);
            gpsi = [gpsi_z, gpsi_UV];
            
            % energy
            a = cpf.lsm; 
            % Globally scale energy values to avoid very small values
            C = cpf.optC;
            e = C*[sqrt(a)*es; epsi];
            
            if nargout > 1
                 % smoothness terms
                 g = C*[sqrt(a)*gs;gpsi];
            end
            
            % stats
            cpf.ES = [cpf.ES; es'*es]; cpf.EPSI = [cpf.EPSI; epsi'*epsi]; 
            x = reshape(U,[],2); y = reshape(V,[],2);
            J = [0,-1;1,0]; s = dot(x',-J*y')./vecnorm(x')./vecnorm(y');

            cpf.AR = [cpf.AR,s'];
            % For debugging
            cpf.fvals.es = es; cpf.fvals.epsi = epsi;
            cpf.grads.gs = gs; cpf.grads.gpsi = gpsi;
            
        end

        % objective: 2nie x 1, gradient 2nie x (2nie + 4nf)
        % For each interior edge:
        % real((dr_i^{-1}(Jeij))^N) - real((dr_j^{-1}(Jeij))^N),
        % imag((dr_i^{-1}(Jeij))^N) - imag((dr_j^{-1}(Jeij))^N),
        function [es,gs_U,gs_V,gs_z] = smoothness_ls(cpf, U, V)            
            n = cpf.n;
            
            %% Init
            CM = cpf.M; x = CM.vertices; nf = CM.nf;
            e2v = CM.edges; nie = CM.nie;
            inner = CM.inner_edges;
            UV = [U;V];
            
            %% Objective
            
            % edges in R^3
            eij_3d = x(e2v(inner,1),:) - x(e2v(inner,2),:);
            
            % edges in respective basis
            eij_i = CM.EBt1*eij_3d(:);
            eij_j = CM.EBt2*eij_3d(:);
            
            % Rotate
            Jeij_i = [-eij_i(nie+1:end); eij_i(1:nie)];
            Jeij_j = [-eij_j(nie+1:end); eij_j(1:nie)];
            
            % U,V of neighboring triangles
            zs = sparse(size(CM.e2t1,1),size(CM.e2t1,2));
            Ai = [CM.e2t1, zs; zs, CM.e2t1];
            Aj = [CM.e2t2, zs; zs, CM.e2t2];
            % Ai, Aj : 4nie x 4nf
            UVi = Ai * UV; 
            UVj = Aj * UV;
            
            
            % For each:
            % UV: 4nie x 1, eij: 2nie x 1, n \in [1,..,6]
            % out:  2nie x 1
            % grad: 2nie x 4nie
            
            if nargout == 1
                outi = CPF.ww(UVi, Jeij_i, n);
                outj = CPF.ww(UVj, Jeij_j, n);
            else % if nargout > 1
                [outi,gradi] = CPF.ww(UVi, Jeij_i, n);
                [outj,gradj] = CPF.ww(UVj, Jeij_j, n);
                
                % each: 2nie * 4nf
                gradi = gradi * Ai;
                gradj = gradj * Aj;
            end
            
            es = cpf.smw*(outi - outj);

            if nargout > 1
                %% Gradient wrt z
                gs_z = sparse(2*nie,2*nie);
                gs_UV = cpf.smw*(gradi - gradj);
                gs_U = gs_UV(:,1:2*nf);
                gs_V = gs_UV(:,2*nf+1:end);
            end
        end        

        
        function [epsi,gpsi_UV, gpsi_z] = constraints_ls(cpf, U, V, z)
            nie = cpf.M.nie;
            
            if nargout == 1
                epsi_c = continuity(cpf, [U;V], z, cpf.cont_n); % Continuity         
                epsi_s = sizing(cpf, U, V);                     % Sizing           
                epsi_cj = conjugacy(cpf, U, V);                 % Conjugacy
                epsi_co = orthogonality(cpf, U, V);             % Orthogonality
                epsi_a = alignment(cpf, [U;V]);                 % Alignment     
                epsi_in = injectivity(cpf, U, V);               % Injectivity
            else % if nargout > 1
                % Continuity
                [epsi_c, gpsi_c_UV, gpsi_c_z] = continuity(cpf, [U;V], z, cpf.cont_n);         

                % Sizing
                [epsi_s, gpsi_s_U, gpsi_s_V] = sizing(cpf, U, V);            
                % Sizing independent of z
                gpsi_s_z = sparse(size(gpsi_s_U,1),2*nie);          
                gpsi_s_UV = [gpsi_s_U, gpsi_s_V];

                % Conjugacy
                [epsi_cj, gpsi_cj_U, gpsi_cj_V] = conjugacy(cpf, U, V);            
                % conjugacy independent of z
                gpsi_cj_z = sparse(size(gpsi_cj_U,1),2*nie);          
                gpsi_cj_UV = [gpsi_cj_U, gpsi_cj_V];

                % orthogonality
                [epsi_co, gpsi_co_U, gpsi_co_V] = orthogonality(cpf, U, V);            
                % orthogonality independent of z
                gpsi_co_z = sparse(size(gpsi_co_U,1),2*nie);          
                gpsi_co_UV = [gpsi_co_U, gpsi_co_V];
                
                % Alignment
                [epsi_a, gpsi_a_UV] = alignment(cpf, [U;V]);            
                gpsi_a_z = sparse(size(gpsi_a_UV,1),2*nie);          

                % Injectivity
                [epsi_in, gpsi_in_U, gpsi_in_V] = injectivity(cpf, U, V);
                gpsi_in_z = sparse(size(gpsi_in_U,1),2*nie);          
                gpsi_in_UV = [gpsi_in_U, gpsi_in_V];
            end
            
            a = sqrt(cpf.lc); b = sqrt(cpf.lsi); c = sqrt(cpf.la); ...
                d = sqrt(cpf.lin); e = sqrt(cpf.lcj); f = sqrt(cpf.lo);
            epsi = [a*epsi_c; b*epsi_s; e*epsi_cj; f*epsi_co; ...,
                c*epsi_a; d*epsi_in];
            if nargout > 1
                gpsi_UV = [a*gpsi_c_UV; b*gpsi_s_UV; ...
                    e*gpsi_cj_UV; f*gpsi_co_UV; ...,
                    c*gpsi_a_UV; d*gpsi_in_UV];
                gpsi_z = [a*gpsi_c_z; b*gpsi_s_z; ...
                    e*gpsi_cj_z; f*gpsi_co_z; ...
                    c*gpsi_a_z; d*gpsi_in_z];
            end
            
            % stats
            cpf.EPSI_C = [cpf.EPSI_C; epsi_c'*epsi_c]; 
            cpf.EPSI_S = [cpf.EPSI_S; epsi_s'*epsi_s]; 
            cpf.EPSI_CJ = [cpf.EPSI_CJ; epsi_cj'*epsi_cj];
            cpf.EPSI_CO = [cpf.EPSI_CO; epsi_co'*epsi_co];            
            cpf.EPSI_A = [cpf.EPSI_A; epsi_a'*epsi_a];
            cpf.EPSI_IN = [cpf.EPSI_IN; epsi_in'*epsi_in];
            
            % For debugging
            cpf.fvals.epsi_c = epsi_c; 
            cpf.fvals.epsi_s = epsi_s;
            cpf.fvals.epsi_cj = epsi_cj;
            cpf.fvals.epsi_co = epsi_co;            
            cpf.fvals.epsi_a = epsi_a;
            cpf.fvals.epsi_in = epsi_in;
            if nargout > 1
                cpf.grads.gpsi_c = [gpsi_c_z, gpsi_c_UV]; 
                cpf.grads.gpsi_s = [gpsi_s_z, gpsi_s_UV]; 
                cpf.grads.gpsi_cj = [gpsi_cj_z, gpsi_cj_UV]; 
                cpf.grads.gpsi_co = [gpsi_co_z, gpsi_co_UV]; 
                cpf.grads.gpsi_a = [gpsi_a_z, gpsi_a_UV]; 
                cpf.grads.gpsi_in = [gpsi_in_z, gpsi_in_UV]; 
            end
        end

        % objective: 4nie x 1, gradient 4nie x (2nie + 4nf)
        % For each interior edge:
        % real((dr_i^{-1}(eij))^N) - real(z_ij),
        % imag((dr_i^{-1}(eij))^N) - imag(z_ij),
        % real((dr_j^{-1}(eij))^N) - real(z_ij),
        % imag((dr_j^{-1}(eij))^N) - imag(z_ij),
        function [epsi_c, gpsi_c_UV, gpsi_c_z] = continuity(cpf, UV, z, n)
            if nargin < 4
                n = cpf.n;
            end
            
            %% Init
            CM = cpf.M; x = CM.vertices; 
            e2v = CM.edges; nie = CM.nie;
            inner = CM.inner_edges;
            
            %% Objective
            
            % edges in R^3
            eij_3d = x(e2v(inner,1),:) - x(e2v(inner,2),:);
            
            % edges in respective basis
            eij_i = CM.EBt1*eij_3d(:);
            eij_j = CM.EBt2*eij_3d(:);
            % U,V of neighboring triangles
            zs = sparse(size(CM.e2t1,1),size(CM.e2t1,2));
            Ai = [CM.e2t1, zs; zs, CM.e2t1];
            Aj = [CM.e2t2, zs; zs, CM.e2t2];
            % Ai, Aj : 4nie x 4nf
            UVi = Ai * UV; 
            UVj = Aj * UV;
            
            
            % For each:
            % UV: 4nie x 1, eij: 2nie x 1, n \in [1,..,6]
            % out:  2nie x 1
            % grad: 2nie x 4nie
            
            if nargout == 1
                outi = CPF.ww(UVi, eij_i, n);
                outj = CPF.ww(UVj, eij_j, n);
            else % if nargout > 1
                [outi,gradi] = CPF.ww(UVi, eij_i, n);
                [outj,gradj] = CPF.ww(UVj, eij_j, n);
                
                % each: 2nie * 4nf
                gradi = gradi * Ai;
                gradj = gradj * Aj;
            end
            
            epsi_c = [...
                outi - [real(z);imag(z)];
                outj - [real(z);imag(z)]];

            if nargout > 1
                %% Gradient wrt z
                zs = sparse(nie,nie); 
                os = speye(nie,nie);
                gpsi_c_z = [...
                    -os, zs; ...
                    zs, -os;...
                    -os, zs; ...
                    zs, -os];

                gpsi_c_UV = [gradi ; gradj];            
            end
        end
        
        
        % objective: 2nf x 1, gradient: 2nf x 4nf
        function [epsi_s, gpsi_s_U, gpsi_s_V] = sizing(cpf, U, V) 
            % Sizing metric, 2nf x 2nf
            g = cpf.g;
            nf = cpf.M.nf;
            
            % f_i(u_{i}) = u_{i}^T G_{i} u_{i} - 1 [where G_{i} is 2x2, and u_{i} is 2x1]
            % d(f_i)/d(u_{i}) = 2 G_{i} u_{i} [is 2x1]
            % gradf : nf x 2nf
            % gradf_i = d(f_i)/d(u_{i}) = 2 (G u)_{i}
            epsi_s = [MESH.NORMV(U, 2, g).^2 - 1; MESH.NORMV(V, 2, g).^2 - 1];
            if nargout > 1
                gu = g*U; gum = 2*[spdiags(gu(1:nf),0,nf,nf),spdiags(gu(nf+1:2*nf),0,nf,nf)];
                gv = g*V; gvm = 2*[spdiags(gv(1:nf),0,nf,nf),spdiags(gv(nf+1:2*nf),0,nf,nf)];
                zs = sparse(nf,2*nf);            
                gpsi_s_U = [gum;zs];
                gpsi_s_V = [zs;gvm];
            end
        end

        % objective: nf x 1, gradient: 2*nf x 4nf
        function [epsi_cj, gpsi_cj_U, gpsi_cj_V] = conjugacy(cpf, U, V) 
            % Sizing metric, 2nf x 2nf
            SO = cpf.SO;
            nf = cpf.M.nf;
            
            % USE shape operator, not g!
            Sv = SO * V; Su = SO * U;
            uSv = U .* (Sv);
            % Only for interior faces            
            epsi_cj = cpf.Bf*(uSv(1:nf) + uSv(nf+1:end));

            if nargout > 1
                gpsi_cj_U = cpf.Bf*[spdiags(Sv(1:nf),0,nf,nf), ...
                    spdiags(Sv(nf+1:end),0,nf,nf)];
                gpsi_cj_V = cpf.Bf*[spdiags(Su(1:nf),0,nf,nf), ...
                    spdiags(Su(nf+1:end),0,nf,nf)];
            end
        end
        
        
        function [epsi_co, gpsi_co_U, gpsi_co_V] = orthogonality(cpf, U, V) 
            if cpf.orth_obj == 1
                [epsi_co, gpsi_co_U, gpsi_co_V] = orthogonality1(cpf, U, V);
            elseif cpf.orth_obj == 2
                [epsi_co, gpsi_co_U, gpsi_co_V] = orthogonality2(cpf, U, V);
            else
                error('wrong orth objective');
            end
        end        
        
        % objective: nf x 1, gradient: 2*nf x 4nf
        function [epsi_co, gpsi_co_U, gpsi_co_V] = orthogonality1(cpf, U, V) 
            % Sizing metric, 2nf x 2nf
            nf = cpf.M.nf;
            
            uSv = U .* V;
            epsi_co = uSv(1:nf) + uSv(nf+1:end);
                   
            if nargout > 1
                gpsi_co_U = [spdiags(V(1:nf),0,nf,nf), ...
                    spdiags(V(nf+1:end),0,nf,nf)];
                gpsi_co_V = [spdiags(U(1:nf),0,nf,nf), ...
                    spdiags(U(nf+1:end),0,nf,nf)];
            end
        end

        % objective: nf x 1, gradient: 2*nf x 4nf
        function [epsi_co, gpsi_co_U, gpsi_co_V] = orthogonality2(cpf, U, V) 
            % Sizing metric, 2nf x 2nf
            nf = cpf.M.nf;
            
            %  f = (x'y)/|x|/|y|
            %   \frac{\partial f}{\partial x} = 
            % 1/(|x||y|) y - 1/(|x|^{3} |y\|) y'x  x
            nU = MESH.normv(reshape(U,[],2)); nV = MESH.normv(reshape(V,[],2));
            nUnV = nU.*nV;
            uSv = U .* V; uSv = uSv(1:nf) + uSv(nf+1:end);
            epsi_co = uSv ./ nUnV;
                   
            if nargout > 1
                nU3nVuSv = nU.^3.*nV./uSv; nV3nUuSv = nV.^3.*nU./uSv;
                ddU1 = 1./nUnV.*V(1:nf) - 1./nU3nVuSv.*U(1:nf);
                ddU2 = 1./nUnV.*V(nf+1:end) - 1./nU3nVuSv.*U(nf+1:end);
                ddV1 = 1./nUnV.*U(1:nf) - 1./nV3nUuSv.*V(1:nf);
                ddV2 = 1./nUnV.*U(nf+1:end) - 1./nV3nUuSv.*V(nf+1:end);
                gpsi_co_U = [spdiags(ddU1,0,nf,nf), ...
                    spdiags(ddU2,0,nf,nf)];
                gpsi_co_V = [spdiags(ddV1,0,nf,nf), ...
                    spdiags(ddV2,0,nf,nf)];
            end
        end
        
        % objective: nf x 1, gradient: nf x 4nf
        function [epsi_a, gpsi_a_UV] = alignment(cpf, UV) 
            if isempty(cpf.A.Ai)
                epsi_a = [];
                gpsi_a_UV = [];
                return
            end
            
            % Indexes which faces we need, 4na x 4nf
            Ai = cpf.A.Ai;
            % Indexes which outputs we need, na x 2na
            Bi = cpf.A.Bi;
            
            % Number of alignment constraints x 4nf
            % Ai: 4na x 4nf, ea: 2na x 1, ns: na x 1
            UVa = Ai * UV; 
            ea = cpf.A.d;
            ns = cpf.A.ns;

            % UVa: 4na x 1, ea: 2na x 1, ns: na x 1
            % out:  2na x 1
            % grad: 2na x 4na
            if nargout == 1
                out = CPF.ww(UVa, ea, ns);
            else % if nargout > 1
                [out,grad] = CPF.ww(UVa, ea, ns);
            end
            
            % Bi: na x 2na (taking only the imaginary parts)
            epsi_a = Bi * out;
            
            if nargout > 1
                % na x 4nf
                gpsi_a_UV = Bi * grad * Ai;
            end
        end
        
        
        % objective: nf x 1, gradient: nf x 4nf
        function [epsi_in, gpsi_in_U, gpsi_in_V] = injectivity(cpf, U, V)

            nf = cpf.M.nf;
            J = [0,-1;1,0];  
            sconst = cpf.inj_br_s;
            epsilon = cpf.inj_br_eps; 

            U2 = reshape(U,[],2)'; V2 = reshape(V,[],2)';     
            s = dot(U2,-J*V2);
            U_prev = reshape(cpf.U_prev,[],2)'; V_prev = reshape(cpf.V_prev,[],2)';     
            s_prev = dot(U_prev,-J*V_prev); 

            si_gal = sconst*s_prev;

            x = s-epsilon;
            g = (1./(si_gal.^3)).*x.^3-(3./(si_gal.^2)).*x.^2+(3./(si_gal)).*x;
            
            phi = zeros(nf,1);
            phi(x<=0) = inf;
            phi(x>0) = 1./g(x>0)-1;
            phi(x>=si_gal) = 0;
            epsi_in = phi;

            if nargout > 1
                dg = (3./(si_gal.^3)).*x.^2-(6./(si_gal.^2)).*x+(3./(si_gal));

                dphi = zeros(nf,1);
                dphi(x<=0) = inf;
                dphi(x>0) = -(1./(g(x>0).^2)).*dg(x>0);
                dphi(x>=si_gal) = 0;

                gu1 = dphi.*(V2(2,:)');
                gu2 = -dphi.*(V2(1,:)');
                gv1 = -dphi.*(U2(2,:)');
                gv2 = dphi.*(V2(1,:)');

                gpsi_in_U = [spdiags(gu1,0,nf,nf),spdiags(gu2,0,nf,nf)];
                gpsi_in_V = [spdiags(gv1,0,nf,nf),spdiags(gv2,0,nf,nf)];
                
            end
            cpf.U_prev = U; cpf.V_prev = V;
        end
        
        
        function [ g ] = comp_gradient( cpf, ey, ~ )
            [~,g] = cpf.comp_energy( ey );
        end
        
        % Initial solution
        function [z0, U0, V0] = opt_init(cpf)
            CM = cpf.M;
            
            % Rand init
            if cpf.optinit == 'r'
                U0 = rand(2*CM.nf,1);
                U0 = MESH.normalize_vf(U0,2);
                nU0g = MESH.NORMV(U0, 2, cpf.g);
                U0 = reshape(U0,[],2)./repmat(nU0g,1,2);
                U0 = U0(:);
                % V0 = J U0
                V0 = [-U0(CM.nf+1:end);U0(1:CM.nf)];
                
                % Normalize fields wrt g
                nU0g = MESH.NORMV(U0, 2, cpf.g);
                U0 = reshape(U0,[],2)./repmat(nU0g,1,2);
                U0 = U0(:);
                nV0g = MESH.NORMV(V0, 2, cpf.g);
                V0 = reshape(V0,[],2)./repmat(nV0g,1,2);
                V0 = V0(:);
                
            % U0 given, assumes U0 has already been set (as an extrinsic
            % vector). U0 is not normalized to length 1.
            elseif cpf.optinit == 'u'
                U0 = CM.EB*cpf.U0(:);
                % V0 = J U0
                V0 = [-U0(CM.nf+1:end);U0(1:CM.nf)];                
            % Smooth
            elseif cpf.optinit == 's'
                op = -CM.D*CM.G;
                [ey,~] = eigs(op,[],10,1e-5);
                U0 = CM.G*ey(:,2);
                U0 = MESH.normalize_vf(U0,3);
                U0 = CM.EB*U0(:);
                % V0 = J U0
                V0 = [-U0(CM.nf+1:end);U0(1:CM.nf)];                

                % Normalize fields wrt g
                nU0g = MESH.NORMV(U0, 2, cpf.g);
                U0 = reshape(U0,[],2)./repmat(nU0g,1,2);
                U0 = U0(:);
                nV0g = MESH.NORMV(V0, 2, cpf.g);
                V0 = reshape(V0,[],2)./repmat(nV0g,1,2);
                V0 = V0(:);
                
            % Guiding field (where available as a 2-rosy) smoothed to all
            % shape for n=6, and curvature directions where not umbilics
            % smoothed to all shape for n=4
            elseif cpf.optinit == 'k'    
                % must align to guiding field for this initialization
                assert(cpf.align_guiding > 0);
                cons2 = cpf.cons2; cons4 = cpf.cons4;
                locs2 = find(MESH.normv(cons2) > 1e-5);
                locs4 = find(MESH.normv(cons4) > 1e-5);
    
                % Guiding fields from connected components (1) or damin (2)
                if cpf.align_guiding == 1 || cpf.align_guiding == 2
                    % For hexes, use guiding field
                    if cpf.n == 6 || cpf.n == 2
                        % If there are parabolic regions, use guding field for 
                        % initialization
                        if ~isempty(locs2)
                            w = smooth_n(cpf, cons2, 2);
                        % Otherwise, curvature is uniform (elliptic or hyperbolic)
                        % If there are locations without umbilics, use dmin there
                        % otherwise will return smoothest 2-field
                        else
                            w0 = zeros(cpf.M.nf,3);
                            w0(locs4,:) = cpf.A.dmin(locs4,:);
                            w = smooth_n(cpf, w0, 2);
                        end
                    elseif cpf.n == 4
                        % locations are mutually exclusive
                        w = smooth_n(cpf, cons2 + cons4, 4, cpf.A.PHI);
                    end
                    cpf.w = w;
                % Optimized guiding field
                elseif cpf.align_guiding == 3
                    assert(~isempty(cpf.w));
                    w = cpf.w;
                    if isempty(locs2) && isempty(locs4)
                        w = smooth_n(cpf, cons2, 2);
                        cpf.w = w;
                    end
                end
                    
                cpf.show_w;
                
                U0 = CM.EB*w(:);
                V0 = [-U0(CM.nf+1:end);U0(1:CM.nf)];                                                
                % Hex alignment change
%               % If aligning V with the prescribed direction use:
%                 V0 = CM.EB*w(:);
%                 U0 = -[-V0(CM.nf+1:end);V0(1:CM.nf)];                                                
                
                % Normalize fields wrt g
                nU0g = MESH.NORMV(U0, 2, cpf.g);
                U0 = reshape(U0,[],2)./repmat(nU0g,1,2);
                U0 = U0(:);
                nV0g = MESH.NORMV(V0, 2, cpf.g);
                V0 = reshape(V0,[],2)./repmat(nV0g,1,2);
                V0 = V0(:);
            % Guiding field (where available as a 2-rosy), damin (where not)
            % KVF sizing
            elseif cpf.optinit == 'l'               
                [cons2, cons4, ~, damin, dmin] = SHAPEOP.guiding_field( CM );
                w = damin;
                locs = find(MESH.normv(cons2)>1e-5);
                w(locs,:) = cons2(locs,:);
                locs = find(MESH.normv(cons4)>1e-5);
                w(locs,: ) = dmin(locs,: );
                U0 = CM.EB*w(:);
                V0 = [-U0(CM.nf+1:end);U0(1:CM.nf)];                                
                
                % Normalize fields wrt kvf
                qu = 1; qv = 1;
                nu = 2*pi*max(cpf.M.vertices(:,1)) / qu;
                nv = 2*pi*max(cpf.M.vertices(:,3)) / qv;
                normu = MESH.normv(cpf.M.Iv2f*cpf.M.vertices(:,1:2));
                normu = 2*pi*normu / nu;
                normv = max(cpf.M.vertices(:,3));
                normv = 2*pi*normv / nv;
                
                nU0 = MESH.NORMV(U0, 2);
                U0 = reshape(U0,[],2)./repmat(nU0,1,2).*normu;
                U0 = U0(:);
                nV0 = MESH.NORMV(V0, 2);
                V0 = reshape(V0,[],2)./repmat(nV0,1,2).*normv;
                V0 = V0(:);
            else
                error('Invalid optinit');
            end

            % pullback edges for all interior edges (both triangles)
            [z0_i, z0_j] = cpf.pullback_edges(U0, V0);
            
            % Set z0 as average of pullback edges of neighboring triangles
            z0 = (z0_i + z0_j)/2;
            z0 = [real(z0); imag(z0)];
            
            % Initialize for injectivity
            cpf.U_prev = U0; cpf.V_prev = V0;  
            J = [0,-1;1,0];
            cpf.inj_br_eps = ...
                min(dot(reshape(U0,[],2), (-J*reshape(V0,[],2)')',2))' / 5;
        end
        
        % pullback a single triangle
        function e2d = pullback_all_triangles(cpf, U, V)
            CM = cpf.M;
            [U,V] = cpf.vfsi(U,V);
            
            Ei = CM.EB*[CM.E1(:), CM.E2(:), CM.E3(:)];
            e2d = zeros(CM.nf*2,3);
            for i=1:3
                e2d(:,i) = CPF.gg([U;V],Ei(:,i));
            end            
        end
        
        % compute the angle defect in the pulled back metric. 
        % Ks : real curvature
        % SA : angle sum
        function [Ks,SA] = angle_defect(cpf, U, V)
            CM = cpf.M;
            nv = CM.nv; nf = CM.nf;
            T = cpf.M.triangles;
            
            e2d = cpf.pullback_all_triangles(U,V);
            
            % edge lengths and angles
            L1 = MESH.normv(reshape(e2d(:,1),[],2));
            L2 = MESH.normv(reshape(e2d(:,2),[],2));
            L3 = MESH.normv(reshape(e2d(:,3),[],2));
            A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
            A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
            A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
            A = [A1,A2,A3];
            A = acos(A);

            % Gaussian curvature
            I = reshape(T,nf*3,1);
            J = ones(nf*3,1);
            S = reshape(A,nf*3,1);
            SA = sparse(I,J,S,nv,1);    
            SA = full(SA);
            Ks = 2*pi*ones(nv,1) - SA;
                        
            Ks(CM.edges(cpf.M.ie == 0,:)) = 0;
        end
        
        % pullback edges for all interior edges (both triangles)
        % z_i, z_j: pulled back edges as complex numbers to the power n, nie x 1
        % eij_i_2d, eij_j_2d: pulled back edges as 2d vectors, nie x 2
        function [z_i, z_j, eij_i_2da, eij_j_2da] = pullback_edges(cpf, U, V)
            CM = cpf.M;
            inner = CM.inner_edges;

            [U,V] = cpf.vfsi(U,V);
            
            z_i = zeros(CM.nie, 1); z_j = zeros(CM.nie,1);
            eij_i_2da = zeros(CM.nie, 2); eij_j_2da = zeros(CM.nie, 2);
            for i=1:CM.nie
                ti = abs(CM.e2t(inner(i),1)); tj = abs(CM.e2t(inner(i),2));
                uvi = [U(ti+[0,CM.nf]),V(ti+[0,CM.nf])];
                uvi = inv(uvi);
                uvj = [U(tj+[0,CM.nf]),V(tj+[0,CM.nf])];
                uvj = inv(uvj);
                eij = CM.vertices(CM.edges(inner(i),1),:) - CM.vertices(CM.edges(inner(i),2),:);
                eij = eij';
                eij_i = CM.EB(ti+[0,CM.nf],ti+[0, CM.nf, 2*CM.nf])*eij;
                eij_j = CM.EB(tj+[0,CM.nf],tj+[0, CM.nf, 2*CM.nf])*eij;
                eij_i_2d = uvi*eij_i;
                eij_j_2d = uvj*eij_j;
                z_i(i) = (eij_i_2d(1) + 1i*eij_i_2d(2)).^cpf.n;
                z_j(i) = (eij_j_2d(1) + 1i*eij_j_2d(2)).^cpf.n;
                
                eij_i_2da(i,:) = eij_i_2d';
                eij_j_2da(i,:) = eij_j_2d';
            end
        end
        

        % iterate LM with increased lc
        function opt_power_cfs( cpf )
            if cpf.profile
                profile on;
            end
            
            CM = cpf.M; 
            
            tic;
            [cpf.z0, cpf.U0, cpf.V0] = opt_init(cpf);
                        
            ey = [cpf.z0;cpf.U0;cpf.V0];

            cpf.U0 = reshape(cpf.M.EBI*cpf.U0(:),[],3);
            cpf.V0 = reshape(cpf.M.EBI*cpf.V0(:),[],3);
           
            opts = [];
            opts.display = 'iter'; % NOT USED
            opts.optTol = cpf.abseps;  % NOT USED
            opts.progTol = cpf.abseps; % NOT USED
            opts.funcTol = cpf.function_tol;
            opts.MaxFunEvals = 1000; % NOT USED
            opts.MaxIter = 1000; % NOT USED
            % Run warm up iterations
            for i=1:cpf.iter_opt
                [epf,output] = OPT.iter(@cpf.comp_energy_ls,...
                                    @cpf.comp_gradient,...
                                    cpf.opttool,ey,opts);
                cpf.lc = cpf.lc*cpf.iter_opt_mult;
                ey = epf;
            end
            [epf,output] = OPT.iter(@cpf.comp_energy_ls,...
                                @cpf.comp_gradient,...
                                cpf.opttool,ey,opts);
            
            cpf.mtime = toc;
            fprintf('\tComputation took: %f sec\n',cpf.mtime);
            
            % stats    
            cpf.epf = epf;
            cpf.miter = output.iterations;
            
            % extract U,V vector fields
            cpf.z = epf(1:2*CM.nie);
            cpf.U = reshape(CM.EBI*epf(2*CM.nie+1:2*CM.nie+2*CM.nf),[],3);
            cpf.V = reshape(CM.EBI*epf(2*CM.nie+2*CM.nf+1:end),[],3);
            
            if cpf.profile
                profile off;
            end
        end
                       
        % Store the resulting texture coordinates in OBJ 
        function store_obj( cpf )
            CM_c = cpf.M_c;
            objname = [cpf.outdir cpf.M.name '.obj'];
            MESH_IO.wobj(objname, CM_c.vertices, CM_c.triangles, ...
                CM_c.texture_vertices, CM_c.texture_triangles, 'cross.png');                   
        end
        
        
        % Main algorithm, including all the required sub-steps :)
        function solve( cpf )            
            % optimize for smooth & aligned power coordinate fields
            fprintf('Optimize for smooth & aligned power coordinate fields\n');
            cpf.opt_power_cfs();
            fprintf('\n');
            
            cpf.plot_energy_terms();

            if cpf.figs
                CM = cpf.M;
                figure; hold on; title('initial vector fields');
                MESH_VIS.vf(CM, cpf.U0, 'color','r');
                MESH_VIS.vf(CM, cpf.V0, 'color','w');
                figure; hold on; title('final vector fields');
                MESH_VIS.vf(CM, cpf.U, 'color','r');
                MESH_VIS.vf(CM, cpf.V, 'color','w');
            end
            
            % Cannot continue to parameterization for n < 4
            if cpf.n < 4
                return
            end
            
            % parameterize and mesh using optimized VFs
            [cpf.Mr, cpf.Mrd] = cpf.uv2mesh(cpf.U, cpf.V);
            
            % Planarize if needed
            if cpf.planar
                if cpf.dualize
                    cpf.Mp = cpf.planarize(cpf.Mrd);
                else
                    cpf.Mp = cpf.planarize(cpf.Mr);
                end
                
                % save the triangulated hex
                % note that this code uses a basic triangulation and can result in intersecting triangles
                outmesh_planar_tri = [cpf.Mp.name '_tri.off'];
                MESH_IO.woffp(outmesh_planar_tri,cpf.Mp.vertices,cpf.Mp.faces);
            end
            
            % save results to OBJ files
            fprintf('Store results in OBJ files\n');
            cpf.store_obj();

            % plot the result
            if cpf.figs
                figure; MESH_VIS.mesh( cpf.Mr ); title('primal');
                if cpf.dualize
                    figure; MESH_VIS.mesh( cpf.Mrd ); title('dual');
                end
                if cpf.planar
                    figure; MESH_VIS.mesh( cpf.Mp ); title('planar');
                end
            end
            
            % Print poisson error
            fprintf('Poisson error (before int) %.2g\n', cpf.perr);
            % Print num sing
            fprintf('%d singularities\n', nnz(cpf.S));
        end
                     
        % parameterize and mesh using U V
        % outputs regular mesh and its dual mesh if needed
        % matching is computed here and given as input to cutter.        
        function [Mr, Mrd] = uv2mesh(cpf, U, V)
            CM = cpf.M;

            %% 1. UV to grads
            fprintf('Convert to gradients\n');
            [gradu, gradv] = ...
                CPF.uv2grads(CM.EB*U(:),CM.EB*V(:),CM.nf);
            % In extrinsic basis
            cpf.gradu = reshape(CM.EBI*gradu,[],3);
            cpf.gradv = reshape(CM.EBI*gradv,[],3);            
                        
            %% 2. Compute matching
            [cpf.mtch,cpf.S] = compute_matching_and_sing(cpf, U, V);
            cpf.vis_matching_all(cpf.U,cpf.V,cpf.mtch,'before combing');
            cpf.vis_sing(cpf.S);
            
            %% 3. cut mesh with directional and prepare for parameterization
            fprintf('Cut mesh\n');
            [cpf.graduc, cpf.gradvc, cpf.M_c, cpf.M_c_ni, cpf.S, ...
                cpf.Fc, cpf.mtchc, cpf.pd] = ...
                    cpf.cut_and_comb(cpf.gradu, cpf.gradv, cpf.mtch, cpf.S);
                
            %% combed U,V
            graduc = reshape(cpf.M.EB*cpf.graduc(:),[],2);
            gradvc = reshape(cpf.M.EB*cpf.gradvc(:),[],2);
            [Uc, Vc] = cpf.grads2uv(graduc, gradvc, cpf.M.nf);
            cpf.Uc = reshape(cpf.M.EBI*Uc(:),[],3); 
            cpf.Vc = reshape(cpf.M.EBI*Vc(:),[],3);
            
            %% Vis matching after combing
            cpf.vis_matching_all(cpf.Uc, cpf.Vc, cpf.mtchc,...
                'after combing',1);

            %% 4. parameterize without integer translations
            % 5. parameterize with integer translations
            % 6. Mesh (and dualize if needed)
            fprintf('Parameterize and mesh\n');
            [Mr,Mrd] = cpf.param_and_mesh;
        end
        
        % For every non-boundary edge, compute the smallest rotation which
        % is a multiple of 2*pi/n and matches a triangle to its neighbor
        % Also computes singularities
        function [mtch,S,SA,K,vtheta] = compute_matching_and_sing(cpf, U, V)
            [~, ~, eij_i, eij_j] = cpf.pullback_edges(U, V);
            
            z_i = eij_i(:,1) + 1i*(eij_i(:,2));
            z_j = eij_j(:,1) + 1i*(eij_j(:,2));

            theta = angle(z_j./z_i);                       
                        
            % Compute singularities
            [~, SA] = cpf.angle_defect(U,V);

            CM = cpf.M; nv = CM.nv;
            theta_all = zeros(CM.ne,1); 
            theta_all(CM.inner_edges) = theta;

            % Reminder after removing 2*pi/n jumps
            theta_rem_all = ...
                theta_all - round(theta_all / (2*pi/cpf.n))*2*pi/cpf.n;

            vtheta = zeros(nv,1);
            for i=1:nv
                if CM.bv(i), continue; end;
                to_v = nonzeros(CM.v2e(:,i));
                to_v_agrees_with_ti = sign(CM.e2t(to_v,1));

                from_v = nonzeros(CM.v2e(i,:)); 
                from_v_agrees_with_ti = sign(CM.e2t(from_v,1));

                vtheta(i) = ...
                    sum(to_v_agrees_with_ti.*theta_rem_all(to_v)) ...
                    - sum(from_v_agrees_with_ti.*theta_rem_all(from_v));
            end

            K = 2*pi - (SA+vtheta);
            S = round(K/(2*pi/cpf.n));           
            S(CM.bv==1) = 0; K(CM.bv==1) = 0;
            assert(norm(K - S*2*pi/cpf.n) < 1e-10);

            % Note that matching depends on orientation!!
            mtch = round(theta_all / (2*pi/cpf.n));                                    
            % Shift to [0,..,cpf.n-1];
            locs = find(mtch<0);
            mtch(locs) = mtch(locs)+cpf.n;
        end
        
        % Vis matching will UV, grads and rawfield
        function vis_matching_all(cpf, U, V, mtch, str, show_bdry)
            if nargin < 6
                show_bdry = 0;
            end
            cpf.vis_matching(U, V, mtch, ...
                ['matching, UV, ' str], show_bdry, cpf.S);

            Ui = cpf.M.EB*U(:); Vi = cpf.M.EB*V(:);
            [gradui, gradvi] = cpf.uv2grads(Ui,Vi,cpf.M.nf);
            gradu = reshape(cpf.M.EBI*gradui,[],3);
            gradv = reshape(cpf.M.EBI*gradvi,[],3);
            cpf.vis_matching(gradu, gradv, mtch, ...
                ['matching, grads ' str], show_bdry, cpf.S);            
            
            F = cpf.grads2rawfield(gradu, gradv);
            cpf.vis_matching(...
                F(:,1:3*cpf.n/2), F(:,3*cpf.n/2+1:3*cpf.n), ...
                mtch,...
                ['matching, rawfield ' str], show_bdry, cpf.S);            
        end
        
        % Visualize the matching
        function vis_matching(cpf, U, V, mtch, str, show_bdry, S)
            if ~cpf.figs
                return
            end
            
            mtch = mtch(cpf.M.inner_edges);
            
            if nargin < 7
                S = [];
            end        
            if nargin < 6
                show_bdry = 0;
            end
            
            figure; hold on; title(str);
            MESH_VIS.mvf(cpf.M, [U,V], 'f', zeros(cpf.M.nf,1));

            if cpf.n == 2
                cols = ['w','m'];
                lgnd = {'0','1'};
            elseif cpf.n == 4
                cols = ['w','r','m','k'];
                lgnd = {'0','1','2','-1'};
            else
                cols = ['w','r','b','m','c','k'];
                lgnd = {'0','1','2','3','-2','-1'};
            end
            ms = [0:cpf.n-1]; 
            xx = zeros(cpf.M.nie,1);
            for i=1:length(ms)
                locs = find(mtch == ms(i));
                lw = 4; if i==1, lw = 1; end;
                if ~isempty(locs)
                    MESH_VIS.edges(cpf.M, ...
                        cpf.M.edges(cpf.M.inner_edges(locs),:), ...
                        'color',cols(i), 'linewidth', lw);
                    xx(locs) = 1;
                end
                MESH_VIS.text_line_annotation(i,'string',...
                    [lgnd{i} ' '], 'color',cols(i));                    
            end
            
            if show_bdry
                MESH_VIS.bdry(cpf.M_c,'color','y','linewidth',1)
            end   
            if ~isempty(S)
                MESH_VIS.singular_pts(cpf.M, .01, S);
            end
        end

        function vis_sing(cpf, S)
            if ~cpf.figs
                return
            end
            
            figure;MESH_VIS.mesh(cpf.M);
            MESH_VIS.singular_pts(cpf.M, .01, cpf.S); title('singularities');           
        end
        
        % use libdirectional to cut the mesh and prepare for
        % parameterization
        % Input:
        % + gradu, gradv: gradients
        % + mtch: edge based matching 
        %
        % Does:  
        % matching 
        %            -> effort (sum of rotation anlges per vertex)
        % effort + matching 
        %            -> singVertices + singIndices (effort_to_indices)
        % singVertices 
        %            -> cut edges (cut_mesh_with_singularities)
        % cut edges + field + matching 
        %            -> combed field + combed matching (comb)
        % combed matching + singVertices + cut edges 
        %            -> cut mesh + paramdata (setup_parameterization)
        % 
        % Output:
        % + graduc, gradvc: combed gradients
        % + M_c, M_c_ni: cut mesh (for final param, and for non int param, resp.)
        % + S: singularities indices
        % + Fc: combed rawfield
        % + mtchc: combed matching (on the cut mesh)
        % + outpd: serialized paramdata file name (for sending to param_and_mesh)
        function [graduc, gradvc, M_c, M_c_ni, S, Fc, mtchc, outpd] = ...
                cut_and_comb(cpf, gradu, gradv, mtch, S)
            CM = cpf.M;
            outmesh = [cpf.outdir CM.name '.out'];
            outmesh_tmp = [cpf.outdir_tmp CM.name '.out'];
            infield = [outmesh_tmp '.rawfield'];
            inmatching = [outmesh_tmp '.mtc'];
            insing = [outmesh_tmp '_s'];

            outmesh_cut = [outmesh_tmp '_cut']; 
            outpd = [outmesh_tmp '_pd'];
            outmatching = [outmesh_tmp '.mtc.c'];             
            outmatching_uncombed = [outmesh_tmp '.mtc.uc']; % FOR DEBUGGNIG
            outfield = [outmesh_tmp '.rawfield.c'];
            outsing = [outmesh_tmp '_s.c']; 
            
            F = cpf.grads2rawfield(gradu, gradv);
            MESH_IO.wrawfield(infield, F);
            
            if ~isempty(mtch)
                [EV, EF, FE, ori] = CM.tolibdirectional;
                fmtch = ori.*mtch;
                MESH_IO.wmtch(inmatching, EV, EF, FE, cpf.n, fmtch);                        
            end
            if ~isempty(S)
                locs = find(S);
                MESH_IO.wsing(insing, locs, S(locs), cpf.n);                                        
            end

            % Cutter
            %   -h,--help                             Print this help message and exit 
            %   -i,--input TEXT:FILE REQUIRED         The path to the OFF file containing the input mesh. 
            %   -f,--input_field TEXT:FILE REQUIRED   The path to the rawfield file containing the vector field. 
            %   -M,--input_matching TEXT:FILE         File with predefined matching. 
            %   -s,--input_sing TEXT:FILE             The path to the sing file (singularities). 
            %   -o,--output TEXT                      The path to the output OFF file. 
            %   -U,--output_serialized_pd TEXT        The output serialized pd 
            %   -m,--output_combed_matching TEXT      Output combed maching. 
            %   -a,--output_matching TEXT             Output uncombed maching. 
            %   -F,--output_combed_filed TEXT         Output combed filed. 
            %   -S,--output_sing_filed TEXT           Path to singularities information. (sign file fomrat). 
            cmdline = sprintf(['!"%s%s" -i "%s%s.off" -f %s ' ...
                ' -o %s.off -U %s -F %s -a %s -m %s -S %s.txt'],...
                cpf.cutterdir, cpf.cutterbin,...
                cpf.datadir, CM.name, infield, ...
                outmesh_cut, outpd, outfield, ...
                outmatching_uncombed, outmatching, outsing);
            if ~isempty(mtch)
                cmdline = sprintf('%s -M %s', cmdline, inmatching);
            end
            if ~isempty(S)
                cmdline = sprintf('%s -s %s', cmdline, insing);
            end
            eval( cmdline );                        
            
            % Output field
            Fc = MESH_IO.rrawfield(outfield);
            [graduc, gradvc] = cpf.rawfield2grads(Fc);
                      
            % Cut meshes (integer and non_integer)
            M_c = MESH(outmesh_cut);
            M_c_ni = MESH(outmesh_cut);
            
            % matching and singularities
            [EV, EF, FE, n, Mt] = MESH_IO.rmtch(outmatching);
            assert(n == cpf.n, 'wrong n when reading matching');
            [mtchc, ori] = CM.fromlibdirectional(EV, EF, FE, Mt);
            locsneg = find(ori < 0);
            mtchc(locsneg) = cpf.n - mtchc(locsneg);
            mtchc = mod(mtchc,6);         
            
            [Sind,Sval] = MESH_IO.rsing([outsing '.txt']);            
            S = zeros(CM.nv,1);
            S(Sind) = Sval;
            S(CM.bv==1) = 0;  % Ignore bdry vertices
        end
        
        
        % Convert parameterization gradients gradu, gradv to rawfield, 
        % which is composed from the gradients of the unit root functions of degree n
        function F = grads2rawfield(cpf, gradu, gradv);
            n = cpf.n;
            
            % recall that D_U, D_V are invariant to change in metric, and
            % U, V are the pushforward of \hat u, \hat v.
            % D_U u = <U, gradu> = <\hat u, \hat u> = 1
            % D_U v = <U, gradv> = <\hat u, \hat v> = 0
            % D_V u = <V, gradu> = <\hat v, \hat u> = 0
            % D_V v = <V, gradv> = <\hat v, \hat v> =  1
            % ==>
            % [gradu; gradv] [U, V] = [1, 0; 0, 1] ==> 
            % [gradu; gradv] = [U,V]^{-1}
            % 
            % Repeat the same for gradh_1 gradh_2 gradh_3 (e.g. for n=6)
            % D_U h_1 = <U, gradh_1> = <\hat u, \hat h_1> = 1
            % D_U h_2 = <U, gradh_2> = <\hat u, \hat h_2> = cos(pi/3)
            % D_U h_3 = <U, gradh_3> = <\hat u, \hat h_3> = -cos(pi/3)
            % D_V h_1 = <V, gradh_1> = <\hat v, \hat h_1> = 0
            % D_V h_2 = <V, gradh_2> = <\hat v, \hat h_2> = sin(pi/3)
            % D_V h_3 = <V, gradh_3> = <\hat v, \hat h_3> = sin(pi/3)
            % 
            % rn = [[1; cos(pi/3); -cos(pi/3)], [0; sin(pi/3); sin(pi/3)]];
            % [gradh_1; gradh_2; gradh_3]_3x2 [U, V]_2x2 = rn_3x2 ==> 
            % [gradh_1; gradh_2; gradh_3] =  rn [U,V]^{-1} ==>
            % [gradh_1; gradh_2; gradh_3] =  rn [grad_u; grad_v] ==>
            %
            % In general, the size of rn is nx2            
            a = linspace(0,2*pi,n+1)';
            rn = [-sin(a), cos(a)];
            Fc = cell(1,n);
            for i=1:n
                Fc{1,i} = gradu*rn(i,1) + gradv*rn(i,2);
            end
            F = cell2mat(Fc);
        end
 
        
        % Convert rawfield (which is composed from the gradients of the 
        % unit root functions of degree n) to [gradu, gradv] - gradients of
        % parameterization functions
        % See comments on grads2rawfields
        function [graduc, gradvc] = rawfield2grads(cpf, F, CM)
            if nargin < 3
                CM = cpf.M;
            end
            n = cpf.n;
            assert(size(F,1) == CM.nf);
            assert(size(F,2) == n*3);
            if n < 4
                error('cannot construct gradients from rawfield for n < 4');
            end
            
            a = linspace(0,2*pi,n+1)';
            rn = [-sin(a), cos(a)];
            
            % first 2 rows of rn should be enough to reconstruct gradients
            rni = inv(rn(1:2,:));   
            % [gradh_1; gradh_2; gradh_3] =  rn [grad_u; grad_v] ==>
            % [gradh_1; gradh_2] =  rn(1:2,:) [grad_u; grad_v] ==>
            % rni [gradh_1; gradh_2] = [grad_u; grad_v]
            grad1 = F(:,1:3); grad2 = F(:,4:6);
            graduc = grad1*rni(1,1) + grad2*rni(1,2);
            gradvc = grad1*rni(2,1) + grad2*rni(2,2);
        end
        
        % Convert parameterization functions (which are composed from the 
        % unit root functions of degree n) to [u, v] - the parameterization functions
        % See comments on grads2rawfields
        function [u, v] = nfuncs2uv(cpf, V, CM)
            if nargin < 3
                CM = cpf.M;
            end
            n = cpf.n;
            assert(size(V,1) == CM.nv);
            assert(size(V,2) == n);
            if n < 4
                error('cannot construct gradients from rawfield for n < 4');
            end
            
            a = linspace(0,2*pi,n+1)';
            rn = [-sin(a), cos(a)];
            
            
            % first 2 rows of rn should be enough to reconstruct functions
            rni = inv(rn(1:2,:));   
            
            A = rni * V(:,1:2)';
            u = A(1,:)'; v = A(2,:)';
        end
         
        
        function vis_param(cpf, M_c, U, V)
            % Visualize
            if cpf.figs
                figure; MESH_VIS.vf(cpf.M, U, 'color','k');
                hold on; MESH_VIS.vf(cpf.M, V, 'color','r');            
                MESH_VIS.func(M_c, M_c.texture_vertices(:,1)); colormap(lines);        

                % Show texture mesh
                M_c.texture_triangles = M_c.triangles;
                figure; MESH_VIS.texture_mesh(M_c);            
            end
        end
       
        
        % Parameterize, given combed field, combed matching, cut mesh,
        % Singularities and pd structure. All computed in cut_and_comb.
        function [Mr, Mrd] = param_and_mesh(cpf)
            CM = cpf.M;
            name = CM.name;           
            outmesh = [cpf.outdir name '.out'];
            outmesh_tmp = [cpf.outdir_tmp name '.out'];
            combedfield = [outmesh_tmp '.rawfield.c'];
            combedmatching = [outmesh_tmp '.mtc.c']; 
            combedsing = [outmesh_tmp '_s.c']; 
            mesh_cut = [outmesh_tmp '_cut']; 
            pd = [outmesh_tmp '_pd'];            
            
            param_out = [cpf.outdir_tmp name '-param.bin'];
            poisson_mat_out = [cpf.outdir_tmp name '_poisson.mat'];
            
             % Cancel out normalization
            if cpf.auto_lr
                nn = zeros(cpf.M.nf,cpf.n);
                F = cpf.grads2rawfield(cpf.gradu, cpf.gradv);
                for i=1:cpf.n
                    nn(:,i) = MESH.normv(F(:,3*(i-1)+1:3*i));
                end
                navg = mean(nn(:));
                cpf.lr = 1/(navg*cpf.M.bbd);
            end
            
            %% Parameterize
%           -h,--help                             Print this help message and exit 
%           -i,--input TEXT:FILE REQUIRED         The path to the OFF file containing the input mesh. 
%           -f,--input_field TEXT:FILE REQUIRED   The path to the rawfield file containing the vector field. 
%           -U,--input_serialized_pd TEXT:FILE    Input PD do be deserialized. 
%           -C,--input_cut_mesh TEXT:FILE         Input cut mesh 
%           -s,--input_sing TEXT:FILE             The path to the sing file (singularities). 
%           -R,--length_ratio FLOAT:NONNEGATIVE REQUIRED 
%                                                 Length ratio. 
%           -M,--input_matching TEXT:FILE         File with predefined matching. 
%           -o,--output TEXT                      The path to the serialized PD. 
%           -L,--output_matlab TEXT=poisson.mat   File name for the matlab output. 
%           -m,--output_matching TEXT             Output combed matching. 
%           -F,--output_combed_filed TEXT         Output combed filed. 
%           -S,--output_sing_filed TEXT           Path to singularities information of the combed field. (sign file fomrat). 
%           -c,--output_cut_mesh TEXT             Output cut mesh. 
%           -N,--normalize_field                  Keep the boundary point at the boundary of the original triangulation. 
%           -I,--use_integers                     Rounding on/off. 
%           -A,--is_combed                        Assume the input field and matching is already combed.
            cmdline = sprintf(['!"%s%s" -i "%s%s.off" -f %s '...
                '-M %s -s %s.txt -C %s.off -U %s '...
                '-R %d -N -A '...
                '-o %s -L %s'], ...
                cpf.paramerdir, cpf.paramerbin,...
                cpf.datadir, name, combedfield, ...
                combedmatching, combedsing, mesh_cut, pd, cpf.lr, ...
                param_out, poisson_mat_out);
            eval( cmdline );   

            %% BarrierPoissonSolve
            fprintf('Barrier Poisson Solve\n');
            BarrierPoissonSolve;
            
            % Poisson error before integers
            cpf.perr = norm(A*x0 - b)/norm(A*x0);
            cpf.x2cornerMat = x2CornerMat;
                                    
            % Convert to texture coords (with only 2 values)
            cpf.M_c_ni.texture_vertices = cpf.x2tv(x0);
            cpf.M_c_ni.texture_triangles = F;

            if success == 1
                cpf.M_c.texture_vertices = cpf.x2tv(xcurr);
                cpf.M_c.texture_triangles = F;
            end
            
            if cpf.figs
                figure; title('parameterization non integer = black, integer=red');
                MESH_VIS.texture_mesh(cpf.M_c_ni,'edgecolor','k');
                if success == 1
                    hold on
                    MESH_VIS.texture_mesh(cpf.M_c,'edgecolor','r');
                end
            end
            
            %% Check int param success
            if success ~= 1
                error('BarrierPoissonSolve failed, cannot find integer param.');
            end
            

            outmesh_ns = [outmesh_tmp '_ns']; 
            outcf = [outmesh_tmp '_cf.csv']; 
            outvf = [outmesh_tmp '_vf.csv']; 

            % Output parameterization result for mesher
            fullx = xcurr(:,1);
            writematrix(fullx, outvf);

            vfunc = cpf.x2cornerMat * xcurr;
            vfunc = reshape(vfunc, cpf.n, [])';            
            
            cfunc = zeros(CM.nf, cpf.n*3);
            for i=1:CM.nf
              for j=1:3
                  cfunc(i:i,N * (j - 1) + 1:N * (j - 1) + N) =...
                      vfunc(cpf.M_c.triangles(i,j),:);
              end
            end
            writematrix(cfunc, outcf);
            
            %% Mesh
            cmdline = sprintf('!"%s%s" "%s%s.off" %s.off %s.off %s.off %s %s %s', ...
                cpf.mesherdir, cpf.mesherbin,...
                cpf.datadir, cpf.M.name, mesh_cut, outmesh_ns, outmesh, ...
                param_out, ...
                outcf, outvf);
            fprintf('Mesher\n');
            eval( cmdline );    

            
            Mr = MESHP(outmesh);                        
            Mrd = [];
            if cpf.dualize
                outmesh_dual = [outmesh '_d'];
                cmdline = sprintf('!"%s%s" %s.off %s.off -b CLIPPED_OLD', ...
                    cpf.dualizerdir, cpf.dualizerbin,...
                    outmesh, outmesh_dual);
                fprintf('Dualizer\n');
                eval( cmdline );    
                Mrd = MESHP(outmesh_dual);            
            end            
        end
        
        % Call planarization
        function Mp = planarize(cpf, M)
            outmesh_planar = [M.name '_p'];

            % Usage: Planarize_bin [OPTIONS] 1 2 3 4 [5] [6] [7]
            % Positionals:
            %   1 TEXT:FILE REQUIRED                  The path to the OFF file containing the input mesh.
            %   2 TEXT:FILE REQUIRED                  The path to the OFF file containing a triangulation of the mesh to be planarized.
            %   3 TEXT REQUIRED                       The path to the output OFF file.
            %   4 FLOAT:NONNEGATIVE REQUIRED          Max. planarity error threshold (in percent) for stopping the optimization.
            %   5 TEXT                                The path to the output CSV file containing normals.
            %   6 TEXT                                The path to a CSV containing planarity per face for each iteration.
            %   7 INT=100                             Max. number of iterations.
            % Options:
            %   -h,--help                             Print this help message and exit
            %   -i,--input TEXT:FILE REQUIRED         The path to the OFF file containing the input mesh.
            %   -t,--input_trig TEXT:FILE REQUIRED    The path to the OFF file containing a triangulation of the mesh to be planarized.
            %   -o,--output TEXT REQUIRED             The path to the output OFF file.
            %   -T,--threshold FLOAT:NONNEGATIVE REQUIRED
            %                                         Max. planarity error threshold (in percent) for stopping the optimization.
            %   -n,--output_normals TEXT              The path to the output CSV file containing normals.
            %   -P,--output_planarity TEXT            The path to a CSV containing planarity per face for each iteration.
            %   -m,--max_iter INT=100                 Max. number of iterations.
            %   -S,--snap_boundary=0                  Keep the boundary points at the boundary of the original triangulation.
            %   -Q,--quad_mode=0                      Add extra energy part that refines the quad mesh.
            %   -s,--sing_symm_mode=0                 Add extra energy part that enforces symmetry of odd degree faces.
            cmdline = sprintf(...
                '!"%s%s" %s.off "%s.off" %s.off %g',...
                cpf.planarizerdir, cpf.planarizerbin,...
                M.name, [cpf.datadir cpf.M.name], ...
                outmesh_planar, cpf.mp);
            if cpf.n==4
                cmdline = sprintf('%s -Q',cmdline);
            end
            if cpf.planarize_sing_cor
               cmdline = sprintf('%s -s',cmdline);
            end
            fprintf('Planarizer\n');
            eval(cmdline);       
            Mp = MESHP(outmesh_planar);
        end
           
        % convert parameterization output to texture vertices
        % Uses x2CornerMat, that's computed in parameterize
        % Note that output might be translated from input, which 
        % may harm integer-ness.
        function tv = x2tv(cpf, x)
            vfunc = cpf.x2cornerMat * x;
            vfunc = reshape(vfunc, cpf.n, [])';            
            [uu,vv] = cpf.nfuncs2uv(vfunc, cpf.M_c);
            tv = [uu,vv];
        end
   
        % Smooth a direction field d using GODF, with n
        % d is assumed to be unit length, and non-zero locations are
        % treated as constraints
        % if we is given, it is used for weighing the constraints
        function w = smooth_n(cpf, d, n, we)
            nf = cpf.M.nf;
            if nargin < 4
                we = [];
            end
            
            locs = find(MESH.normv(d) > 1e-5);            
            nl = size(locs,1);
            Aeq = ...
                sparse([1:2*nl], [locs, nf+locs], ones(2*nl,1),2*nl, 2*nf);
            % -> local basis -> power n -> locs
            beq = reshape(cpf.ff(cpf.M.EB*d(:),n),[],2);
            beq = beq(locs,:); 
            beq = beq(:);

            % weigh by we if given
            if ~isempty(we)
                we = [we(locs);we(locs)];
                We = spdiags(we, 0, 2*nl, 2*nl);
                Aeq = We*Aeq;
                beq = We*beq;
            end
            
            C = cpf.M.godf(n);
            d = zeros(nf*2,1);

            % If there are no constraints, use the eigenvector
            if size(Aeq,1) > 0
                x = lsqlin(C, d, [], [], Aeq, beq);
            else
                [x,~] = eigs(C, 1, 'SM');
            end

            % -> sqrt n -> extrinsic
            w = reshape(cpf.M.EBI*cpf.ff(x,1/n),[],3);
            % Normalize
            w = MESH.normalize_vf(w);
        end
        
        % Smooth a direction field d using GODF, with n=2
        % d2 are n=2 constraints
        % d4 are n=4 constraints for w^2
        function w = smooth_n2_n4(cpf, d2, d4)
            nf = cpf.M.nf;
            
            locs2 = find(MESH.normv(d2) > 1e-5);
            locs4 = find(MESH.normv(d4) > 1e-5);
            nl2 = size(locs2,1);
            nl4 = size(locs4,1);
            
            if nl2 == 0
                error('no N=2 constraints');
            end

            % initialize with n=2 constraints only
            w0 = smooth_n(cpf, d2, 2);
            
            if nl4 == 0
                w = w0;
                return
            end
            
            Aeq2 = ...
                sparse([1:2*nl2], [locs2, nf+locs2], ones(2*nl2,1),2*nl2, 2*nf);
            % -> local basis -> power n -> locs
            beq2 = reshape(cpf.ff(cpf.M.EB*d2(:),2),[],2);
            beq2 = beq2(locs2,:); 
            beq2 = beq2(:);
            
            Aeq4 = ...
                sparse([1:2*nl4], [locs4, nf+locs4], ones(2*nl4,1),2*nl4, 2*nf);
            % -> local basis -> power n -> locs
            beq4 = reshape(cpf.ff(cpf.M.EB*d4(:),4),[],2);
            beq4 = beq4(locs4,:); 
            beq4 = beq4(:);
            
            [~,L2] = cpf.M.godf(2);
            
            % normalize operators
            [~,dd] = eigs(L2'*L2,2,'SM');
            L2_units = dd(1,1);
            Aeq4_units = 2*nl4;
            Aeq2_units = 2*nl2;
            
            L2 = 1./sqrt(L2_units)*L2;
            Aeq4 = Aeq4 / (sqrt(Aeq4_units)); beq4 = beq4 / (sqrt(Aeq4_units));
            if nl2 ~= 0
                Aeq2 = Aeq2 / (sqrt(Aeq2_units)); beq2 = beq2 / (sqrt(Aeq2_units));
            end
            
            gfo.L2 = L2; gfo.Aeq2 = Aeq2; gfo.Aeq4 = Aeq4;
            gfo.beq2 = beq2; gfo.beq4 = beq4;
            
            w0 = cpf.ff(cpf.M.EB*w0(:),2);
            w = OPT.lsqnonlin(...
                @(x) guiding_field_objective(cpf, x, gfo), w0);
            
            % -> sqrt n -> extrinsic
            w = reshape(cpf.M.EBI*cpf.ff(w,1/2),[],3);
            % Normalize
            w = MESH.normalize_vf(w);            
        end
        
        % x*'L2*x + ||Aeq2*x - beq2||^2 + ||Aeq4*x^2 - beq4||^2
        function [e,g] = guiding_field_objective(cpf, ey, gfo)
            L2 = gfo.L2; Aeq2 = gfo.Aeq2; beq2 = gfo.beq2;
            Aeq4 = gfo.Aeq4; beq4 = gfo.beq4;
            
            if nargin > 1
                [ey2, g4] = cpf.ff(ey,2);
            else
                ey2 = cpf.ff(ey,2);
            end
            
            w2 = 1e5; w4 = 1e5;
            
            e = [L2*ey; ...
                sqrt(w2)*(Aeq2*ey - beq2); ...
                sqrt(w4)*(Aeq4*ey2 - beq4)];
            
            if nargin > 1
               g = [L2; ...
                   sqrt(w2)*Aeq2; ...
                   sqrt(w4)*Aeq4*g4];                   
            end            
        end
        
        function plot_energy_terms( cpf )
            if cpf.figs
                figure; hold on; title('energy'); set(gcf, 'WindowStyle', 'docked');

                a = cpf.lsm; b = cpf.lc; c = cpf.lsi; d = cpf.la; e = cpf.lin; f = cpf.lcj; h = cpf.lo;
                l = round(length(cpf.ES)/2);
                nn = min(length(cpf.ES), length(cpf.EPSI));
                has_alignment = ...
                    cpf.align_boundary + cpf.align_curvature + cpf.align_guiding;
                if has_alignment
                    has_alignment = ~isempty(cpf.A.f);
                end
                if has_alignment == 0 
                    cpf.EPSI_A = zeros(nn,1);
                end
                plot([a*cpf.ES(1:nn), b*cpf.EPSI_C(1:nn), c*cpf.EPSI_S(1:nn),...
                    f*cpf.EPSI_CJ(1:nn), d*cpf.EPSI_A(1:nn), e*cpf.EPSI_IN(1:nn), h*cpf.EPSI_CO(1:nn)], ...
                    'linewidth',3);
                ylim([0,a*cpf.ES(l)]); 
                legend('smoothness','continuity','size','conjugacy','alignment','injectivity','orthogonality');

                hold off
            end
        end
        
        function [cons2, cons4, w, damin, dmin, RHO, PHI, locse] = ...
                optimized_guiding_field(cpf)
            % Calls guiding_field2 and returns a 2-rosy that interpolates the
            % cons2 and cons4 constraints. Also returns the full field.
            [cons2, cons4, ~, damin, dmin, RHO, PHI, locse] ...
                = SHAPEOP.guiding_field2(cpf.M);
            
            locs2 = find(MESH.normv(cons2) > 1e-5);
            locs4 = find(MESH.normv(cons4) > 1e-5);

            if cpf.ell_cons_guiding == 0
                locs4 = setdiff(locs4, locse);
            end
            
            if cpf.figs
                ds = {cons2(locs2,:), cons4(locs4,:)};
                fs = {locs2, locs4};
                ns = [2,4];
                figure; MESH_VIS.mesh(cpf.M); hold on;
                MESH_VIS.vf(cpf.M, ds,'locs', fs, 'nRosy', ns,'color', 'r');
                title('constraints, before optimization=red, after=blue');
            end
            
            % If there are parabolic constraints
            if ~isempty(locs2)
                w = cpf.smooth_n2_n4(cons2, cons4);
            % Curvature is uniform (elliptic or hyperbolic)
            % use dmin
            else
                w0 = zeros(cpf.M.nf,3);
                w0(locs4,:) = dmin(locs4,:);
                w = smooth_n(cpf, w0, 2);                
            end
                       
            
            % Switch 2-constraints to interpolated field
            % Also add 4-constraints, and remove all 4-constraints
            cons2(locs2,:) = w(locs2,:);
            cons2(locs4,:) = w(locs4,:);            
            cons4(locs4,:) = 0;            

            if cpf.figs
                locs2 = find(MESH.normv(cons2) > 1e-5);
                locs4 = find(MESH.normv(cons4) > 1e-5);
                ds = {cons2(locs2,:), cons4(locs4,:)};
                fs = {locs2, locs4};
                ns = [2,4];
                MESH_VIS.vf(cpf.M, ds,'locs', fs, 'nRosy', ns,'color', 'b');
            end
        end
        
        function show_w(cpf, str)
            if nargin < 2
                str = '';
            end
            if cpf.figs
                if cpf.n == 6
                    figure; 
                    MESH_VIS.mvf(cpf.M,[cpf.cons2, cpf.w]);
                elseif cpf.n == 4
                    figure; 
                    MESH_VIS.mvf(cpf.M,[cpf.cons2+cpf.cons4, cpf.w],...
                        'nrosy',4);
                end
                title([str ' Initial U from smoothly completed '...
                    'guiding field'], 'interpreter','none');
            end            
        end
   
        function show_planarity(cpf, M)
            if nargin < 2
                M = cpf.Mp;
            end
            pl = M.planarity_general;
            figure; MESH_VIS.func(M, pl);
            title(['planarity ' M.name]);
        end
        
        function show_alignment(cpf, A, M)
            if nargin < 3
                M = cpf.M;
            end
            if cpf.figs
                locs2 = find(A.no == 2);
                locs4 = find(A.no == 4);
                locs6 = find(A.no == 6);

                ds = {A.do(locs2,:), A.do(locs4,:), A.do(locs6,:)};
                fs = {A.f(locs2), A.f(locs4), A.f(locs6)};
                ns = [2,4,2];
                if ~isempty(M)
                    figure; MESH_VIS.mesh(M); hold on;
                end
                MESH_VIS.vf(cpf.M, ds,'locs', fs, 'nRosy', ns,...
                    'color', 'r', 'mesh', 0);
                title('Alignment constraints');
            end
        end
        
        function show_final_fields(cpf, M)
            if nargin < 2
                M = cpf.M;
            end
            if cpf.figs
                figure; MESH_VIS.mesh(M); hold on;
                MESH_VIS.vf(cpf.M, cpf.U, 'color', 'r', 'mesh', 0);
                MESH_VIS.vf(cpf.M, cpf.V, 'color', 'b', 'mesh', 0);
                title('Final fields');
            end
        end

        function show_initial_fields(cpf, M)
            if nargin < 2
                M = cpf.M;
            end
            if cpf.figs
                figure; MESH_VIS.mesh(M); hold on;
                MESH_VIS.vf(cpf.M, cpf.U0, 'color', 'r', 'mesh', 0);
                MESH_VIS.vf(cpf.M, cpf.V0, 'color', 'b', 'mesh', 0);
                title('Initial fields');
            end
        end
        
        function show_final_energy(cpf, name)
            nf = cpf.M.nf; 
            f = eval(['cpf.fvals.' name]);
            if length(f) == 2*nf
                f = reshape(f,[],2);
                for j=1:2
                    figure; MESH_VIS.func(cpf.M, f(:,j).^2);
                    title(['final ' name ' energy ' num2str(j)]);
                end
            else
                error('not supported yet');
            end
        end
        
        function show_curvature(cpf)
            if cpf.figs
                e = cpf.k_eps; r = cpf.delta;
                [~,S] = SHAPEOP.curvature_sizing(cpf.M, e, r);
                figure; MESH_VIS.func(cpf.M,S.kminf); 
                title([cpf.name ' kmin'], 'interpreter','none');
                figure; MESH_VIS.func(cpf.M,S.kmaxf); 
                title([cpf.name ' kmax'], 'interpreter','none');                
                figure; MESH_VIS.func(cpf.M,abs(S.kmaxf./S.kminf)); 
                title([cpf.name ' kmax/kmin'], 'interpreter','none');
                caxis(gca,[-2,2]);
            end            
        end
        
        % Input: U,V either extrinsic, intrinsic or flattened
        % Output: U,V intrinsic
        function [U,V] = vfsi(cpf,U,V)
            if size(U,2) == 3
                U = cpf.M.EB*U(:);
                V = cpf.M.EB*V(:);
            elseif size(U,2) == 2
                assert(size(U,1) == cpf.M.nf);
                assert(size(V,1) == cpf.M.nf);
                U = U(:); V = V(:);
            elseif size(U,2) == 1
                assert(size(U,1) == 2*cpf.M.nf);
                assert(size(V,1) == 2*cpf.M.nf);
            else
                error('?');
            end                      
        end                    
        
    end
    
    methods (Static)
   
        % [gradu(i,:); gradv(i,:)] = [U(i,:)', V(i,:)']^{-1}
        function [gradu, gradv] = uv2grads(U,V,nf)
            Is = ones(nf,1); zs = zeros(nf,1);
            UVi_1 = CPF.gg([U(:); V(:)], [Is; zs]);
            UVi_2 = CPF.gg([U(:); V(:)], [zs; Is]);
            gradu = [UVi_1(1:nf); UVi_2(1:nf)];
            gradv = [UVi_1(nf+1:end); UVi_2(nf+1:end)];
        end
        
        % [U(i,:)', V(i,:)'] = [gradu(i,:); gradv(i,:)]^{-1}
        function [U, V] = grads2uv(gradu,gradv,nf)
            Is = ones(nf,1); zs = zeros(nf,1);
            ggi_1 = CPF.gg([gradu(:,1); gradv(:,1); gradu(:,2); gradv(:,2)], ...
                [Is; zs]);
            ggi_2 = CPF.gg([gradu(:,1); gradv(:,1); gradu(:,2); gradv(:,2)], ...
                [zs; Is]);
            U = reshape(ggi_1,[],2);
            V = reshape(ggi_2,[],2);
        end
      
        % in: 2n x 1
        % out: 2n x 1, f(x,n) = x.^n
        % grad: 2n x 2n
        function [e,g] = ff(in,n)
            s = size(in,1)/2;
            a = in(1:s); b = in(s+1:end);
            c = a+1i*b;
            cn = c.^n;
            e = [real(cn); imag(cn)];
            gc = n.*c.^(n-1);
            if nargout > 1
                g = CPF.c2m(gc);
            end
        end

        % in: 4n x 1, e: 2n x 1
        % out:  2n x 1, g(x,e) = x^{-1}*e (inv of each 2x2 mat)
        % grad: 2n x 4n
        function [out,grad] = gg(in,e)
            in = reshape(in,[],4); e = reshape(e, [], 2);
            n = size(in,1);
            x = in(:,1:2)'; y = in(:,3:4)'; e = e';
            J = [0,-1;1,0];

            s = dot(x,-J*y);
            g = 1./s .* [dot(-y, J' * e); dot(x, J' * e)];
                
            if nargout > 1
                % inv([x y]) = 1/s [-Jy; Jx]
                % Each is 2xn
                % deriv of g(1) wrt x(1) and x(2)
                grg1x = repmat(-g(1,:)./s,2,1) .* (-J*y);
                % deriv of g(2) wrt x(1) and x(2)
                grg2x = repmat(-g(1,:)./s,2,1) .* (J*x);

                % deriv of g(1) wrt y(1) and y(2)
                grg1y = repmat(-g(2,:)./s,2,1) .* (-J*y);
                % deriv of g(2) wrt y(1) and y(2)
                grg2y = repmat(-g(2,:)./s,2,1) .* (J*x);


                % Size 2 x 4n
                gr = [grg1x(1,:),grg1x(2,:),grg1y(1,:),grg1y(2,:);...
                    grg2x(1,:),grg2x(2,:),grg2y(1,:),grg2y(2,:);];

                % Make sparse matrix
                II = [1:n,1:n,1:n,1:n,...
                      n+1:2*n,n+1:2*n,n+1:2*n, n+1:2*n]';
                JJ = [1:4*n, ...
                      1:4*n]';
                SS = [gr(1,:), gr(2,:)]';
                %%% USE FAST MATRIX CONSTRUCTION
                grad = fast_sparse(II,JJ,SS, 2*n, 4*n);
    %            g = 1./s * [-y, x]' * (J' * e);

    %             dgdx - e'*J*y/s^2*[-y x]'*J'
    %             dgdy - (-e'*J*x/s^2*[-y x]'*J')
    %             out = inv([x y])*e;
    %             grx = -[1,0]*g*inv([x y]);
    %             gry = -[0,1]*g*inv([x y]);
            end
            
            g = g';
            out = g(:);
        end
        
        % in: 4n x 1, e: 2n x 1, n \in [1,..,6]
        % out:  2n x 1, w(x,e,n) = f(g(in,e),n)
        % grad: 2n x 4n
        % w = f \circ g
        function [out, grad] = ww(in, e, n)
            if nargout > 1
                [outg, gradg] = CPF.gg(in, e);
                [outf, gradf] = CPF.ff(outg, n);
            else
                outg = CPF.gg(in, e);
                outf = CPF.ff(outg, n);                
            end
            
            out = outf;
            if nargout > 1
                grad = gradf*gradg;
            end
        end
        

        
        % converts a vector of complex numbers to a block matrix
        % with 2x2 blocks of real numbers. 
        % Each complex number c_i becomes:
        % [real(c_i), -im(ci); im(ci), real(ci)]
        % Each block is placed at {(i,i), (i,n+i), (n+i,i), (n+i,n+i) }
        function m = c2m(c)
            n = length(c);
            re = real(c); im = imag(c);
            II = [1:n, 1:n, n+1:2*n, n+1:2*n]';
            JJ = [1:n, n+1:2*n, 1:n, n+1:2*n]';
            SS = [re; -im; im; re];
            m = sparse(II, JJ, SS, 2*n, 2*n);    
        end
        
        % converts a 2x2 block matrix to a vector of real numbers
        % Each block is placed at {(i,i), (i,n+i), (n+i,i), (n+i,n+i) }
        % output is of size 2n x 2, where input matrix has n blocks
        function v = m2v(m)
            n = size(m,1)/2;
            II = [1:n, 1:n, n+1:2*n, n+1:2*n]';
            JJ = [1:n, n+1:2*n, 1:n, n+1:2*n]';
            SS = m(sub2ind(size(m), II, JJ));
            v = [SS(1:2*n), SS(2*n+1:end)];
        end
                   
    end
end
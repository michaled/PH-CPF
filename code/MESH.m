classdef MESH < handle
    
    properties (Access='public')
        
        name
        
        vertices
        texture_vertices
        triangles
        texture_triangles
        edges
        
        nv              % #vertices
        nf              % #faces
        ne              % #edges
        nie             % #interior edges
        
        va              % vertex areas
        ta              % face areas
        ea              % edge areas
        
        VA, IVA         % (lumped) vertices mass matrix and inverse
        TA, ITA         % (lumped) triangles mass matrix and inverse
        EA              % (lumped) edges mass matrix
        
        Iv2f            % interpolate vertices to faces
        If2v            % interpolate faces to vertices
        
        Nf              % normal-per-face
        Nv              % normal-per-vertex
        Nv_IO = 1       % whether to read normals from file
        R               % in-plane rotation operator wrt normal
                
        G               % gradient on vertex functions
        D               % divergence on face fields
        C               % curl of face fields
        RG
        
        mel             % mean edge length
        
        E1, E2, E3      % triangle edges
        F1, F2          % frame per face
        EB, EBI         % an edge based basis for vfs.
        EBt1, EBt2
        
        GE              % gradient on edge functions
        DE              % divergence on face fields
        
        Ie2f            % interpolate edges to faces
        If2e            % interpolate faces to edges
        
        e2t, t2e, v2e   % map edges to triangles, triangles to edges  
                        % vertices to edges
        e2t1, e2t2      % sparse matrices mapping edges to left and 
                        % right triangles on inner edges.

        VFI,VFJ
        VEI,VEJ

        ie              % 1 if interior edge
        inner_edges     % indices of interior edges
        bv              % boundary vertices
        bf              % faces which have boundary vertices
        
        VL, VLn         % smoothness operators, vector field and rosy field
        VLh             % "half" smoothness operator, so ||VLh x||^2 = x' VL x
    end
    
    methods
        
        function [ mesh ] = MESH( meshname, nospectral )
            
            if nargin < 2; nospectral = 0; end;
            if nargin < 1; meshname = 'sphere_s3'; end
                        
            mesh.name = meshname;
            
            [mesh.vertices, mesh.triangles] = MESH_IO.roff([meshname '.off']);
                        
            mesh.nv = size(mesh.vertices,1);
            mesh.nf = size(mesh.triangles,1);
            
            mesh.Nf = face_normals( mesh );
            mesh.ta = face_areas( mesh );
            mesh.Nv = vertex_normals( mesh );
            mesh.va = vertex_areas( mesh );
            
            [mesh.E1,mesh.E2,mesh.E3] = face_edges( mesh );
            mesh.R = rot( mesh );
            
            mesh.G = grad( mesh );
            mesh.D = div( mesh );
            mesh.RG = mesh.R*mesh.G;
            
            mesh.If2v = f2v( mesh );
            mesh.Iv2f = v2f( mesh );
            
            [mesh.edges,mesh.e2t,mesh.t2e, ...
                mesh.e2t1, mesh.e2t2, ...
                mesh.v2e,mesh.ie, ...
                mesh.ne, mesh.nie, mesh.inner_edges, mesh.bv, mesh.bf] = ...
                nc_data( mesh );

            mesh.ea = edge_areas( mesh );
            [mesh.VA,mesh.IVA,mesh.TA,mesh.ITA,mesh.EA] = mass_matrices( mesh );
            
            mesh.GE = grade( mesh );
            mesh.DE = dive( mesh );
            
            mesh.C = curl( mesh );
            
            mesh.If2e = f2e( mesh );
            mesh.Ie2f = e2f( mesh );
            
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            mesh.VFI = I; mesh.VFJ = J;
            
            I = repmat(1:mesh.nie,1,3); I = I(:); 
            J = (1:mesh.nie)'; J = [J; J+mesh.nie; J+2*mesh.nie];
            mesh.VEI = I; mesh.VEJ = J;
            
            [mesh.F1, mesh.F2, mesh.EB, mesh.EBI] = edge_basis( mesh );
            
            % For non manifold meshes
            if ~nospectral
                % smooth frame based on GODF
                [mesh.VL, mesh.VLh] = mesh.godf( 1 );
                [mesh.F1, mesh.F2] = smooth_frame( mesh, 1 ); 

                % GODF for n-RoSy fields
                mesh.VLn = mesh.godf( 4 );
            end
        end
        
        function [ edges, e2t, t2e, ...
                e2t1, e2t2, ...
                v2e, ie, ...
                ne, nie, inner_edges, bv, bf] = nc_data( mesh )
            T = double( mesh.triangles );
            
            
            I = [T(:,2);T(:,3);T(:,1)];
            J = [T(:,3);T(:,1);T(:,2)];
            S = [1:mesh.nf,1:mesh.nf,1:mesh.nf];
            E = sparse(I,J,S,mesh.nv,mesh.nv);
            
            Elisto = [I,J];
            sElist = sort(Elisto,2);
            s = (MESH.normv(Elisto - sElist) > 1e-12);
            t = S'.*(-1).^s;
            [edges,une] = unique(sElist, 'rows');
            
            ne = size(edges,1);
            e2t = zeros(ne,4);
            t2e = zeros(mesh.nf,3);
            ie = zeros(ne,1);
            for m=1:length(edges)
                i = edges(m,1); j = edges(m,2);
                t1 = t(une(m));
                t2 = -(E(i,j) + E(j,i) - abs(t1))*sign(t1);
                e2t(m,1:2) = [t1, t2];
                f = T(abs(t1),:); loc = find(f == (sum(f) - i - j));
                t2e(abs(t1),loc) = m*sign(t1);
                e2t(m,3) = loc;
                if t2 ~= 0
                    f = T(abs(t2),:); loc = find(f == (sum(f) - i - j));
                    t2e(abs(t2),loc) = m*sign(t2);
                    e2t(m,4) = loc;
                    ie(m) = 1;
                end
            end

            v2e = sparse(edges(:,1),edges(:,2),1:length(edges),mesh.nv,mesh.nv);            
            
            ne = size(edges,1);
            nie = sum(ie);
            inner_edges = find(ie);
            bv = zeros(mesh.nv,1);
            bv(edges(ie == 0,:)) = 1;
            bf = zeros(mesh.nf,1);
            bf(sum(ismember(mesh.triangles, find(bv==1)),2)>0) = 1;
            
            t1 = abs(e2t(inner_edges,1)); t2 = abs(e2t(inner_edges,2));
            
            I = 1:2*nie;
            S = ones(2*nie,1);
            e2t1 = sparse(I, [t1; t1+mesh.nf], S, 2*nie, 2*mesh.nf);
            e2t2 = sparse(I, [t2; t2+mesh.nf], S, 2*nie, 2*mesh.nf);
        end
        
        function [ ta ] = face_areas( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            ta = mesh.normv( cross( P1, P2 ) ) / 2;
        end
        
        function [ va ] = vertex_areas( mesh )
            va = full( sum( mass_matrix(mesh), 2 ));
        end
        
        function [ M ] = mass_matrix( mesh )
            T = double( mesh.triangles ); 
            
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            Mij = 1/12*[mesh.ta; mesh.ta; mesh.ta];
            Mji = Mij;
            Mii = 1/6*[mesh.ta; mesh.ta; mesh.ta];
            In = [I;J;I];
            Jn = [J;I;I];
            Mn = [Mij;Mji;Mii];
            M = sparse(In,Jn,Mn,mesh.nv,mesh.nv);
        end
        
        function [ ea ] = edge_areas( mesh )
%             % Hodge star
%             L1 = MESH.normv( mesh.E1 );
%             L2 = MESH.normv( mesh.E2 );
%             L3 = MESH.normv( mesh.E3 );
%             A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
%             A2 = (L3.^2 + L1.^2 - L2.^2) ./ (2.*L3.*L1);
%             A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
%             A = acos( [A1; A2; A3] );
                        
            % *1 - cot weights
            T = double( mesh.triangles );
            I = [T(:,2);T(:,3);T(:,1)];
            J = [T(:,3);T(:,1);T(:,2)];
%             S = 0.5*cot(A);
            S = repmat(mesh.ta/3,3,1);
            In = [I;J];
            Jn = [J;I];
            Sn = [S;S];
            W = sparse(In,Jn,Sn,mesh.nv,mesh.nv);
            ea = zeros(length(mesh.edges),1);
            for i=1:length(ea)
                ea(i) = W(mesh.edges(i,1),mesh.edges(i,2));
                ea(i) = ea(i) + W(mesh.edges(i,2),mesh.edges(i,1));
            end
        end
        
        function [ VA, IVA, TA, ITA, EA ] = mass_matrices( mesh )
            sv = mesh.nv; sf = mesh.nf; se = mesh.ne; 
            VA = spdiags(mesh.va,0,sv,sv);
            IVA = spdiags(1./mesh.va,0,sv,sv);
            TA = spdiags([mesh.ta; mesh.ta; mesh.ta],0,3*sf,3*sf);
            ITA = spdiags(1./[mesh.ta; mesh.ta; mesh.ta],0,3*sf,3*sf);
            EA = spdiags(mesh.ea,0,se,se);
        end
        
        function [ Nf ] = face_normals( mesh )
            
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            Nf = cross( P1, P2 );
            Nf = MESH.normalize_vf( Nf );
        end
        
        function [ Nv ] = vertex_normals( mesh )
            if exist([mesh.name '_normals.off']) == 2
                [Nv, ~] = MESH_IO.roff([mesh.name '_normals.off']);
                fprintf('Read mesh normals from file\n');
                return;
            elseif mesh.Nv_IO & exist([mesh.name '_vertex_normals.csv']) == 2
                Nv = MESH_IO.rNvcsv([mesh.name '_vertex_normals.csv'], mesh.nv);
                fprintf('Read mesh normals from file\n');
                return;
            end
            
            TA = spdiags(mesh.ta,0,mesh.nf,mesh.nf);
            
            I = double( repmat(mesh.triangles(:),3,1) );
            J = repmat(1:3,3*mesh.nf,1); J = J(:);
            
            S = repmat(TA*mesh.Nf,3,1); S = S(:);
            
            Nv = full( sparse(I,J,S,mesh.nv,3) );
            Nv = MESH.normalize_vf( Nv );
        end
        
        function [ E1, E2, E3 ] = face_edges( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            E1 = X(T(:,3),:) - X(T(:,2),:);
            E2 = X(T(:,1),:) - X(T(:,3),:);
            E3 = X(T(:,2),:) - X(T(:,1),:);
            E = [E1; E2; E3];
            
            mesh.mel = mean( MESH.normv( E ) );
        end
        
        function [ R ] = rot( mesh )
            sf = mesh.nf; n = mesh.Nf;
            
            II = repmat((1:sf)',1,2); II = [II II+sf II+2*sf]'; II = II(:);
            JJ1 = (1:sf)'; JJ2 = JJ1+sf; JJ3 = JJ2+sf;
            JJ = [JJ2 JJ3 JJ1 JJ3 JJ1 JJ2]'; JJ = JJ(:);
            SS = [-n(:,3) n(:,2) n(:,3) -n(:,1) -n(:,2) n(:,1)]'; SS = SS(:);
            
            R = sparse(II,JJ,SS,3*sf,3*sf);
        end
        
        function [ G ] = grad( mesh )
            % G corresponds to eq. (3.9) in Polygon mesh processing book
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];

            J = double( mesh.triangles' );
            JJ = [J(:); J(:); J(:)];

            RE1 = mesh.R * mesh.E1(:);
            RE2 = mesh.R * mesh.E2(:);
            RE3 = mesh.R * mesh.E3(:);

            S = [RE1(:) RE2(:) RE3(:)]'; SS = S(:);
            
            G = sparse(II,JJ,SS,3*mesh.nf,mesh.nv);
            ITA = spdiags(.5*repmat(1./mesh.ta,3,1),0,3*mesh.nf,3*mesh.nf);

            G = ITA*G;
        end
        
        function [ D ] = div( mesh )
            % D corresponds to eq. (3.12) in Polygon mesh processing book
            IVA = spdiags(1./mesh.va,0,mesh.nv,mesh.nv);
            TAC = spdiags(repmat(mesh.ta,3,1),0,3*mesh.nf,3*mesh.nf);

            D = - IVA * mesh.G' * TAC;
        end
        
        function [ C ] = curl( mesh )            
            C = - mesh.DE * mesh.R;
        end
        
        function [ GE ] = grade( mesh )
            X = mesh.vertices;
            E = mesh.edges;
            
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];

            J = double( abs(mesh.t2e)' );
            JJ = [J(:); J(:); J(:)];
            
            f1 = repmat( sign(mesh.t2e(:,1)), 1, 3 );
            f2 = repmat( sign(mesh.t2e(:,2)), 1, 3 );
            f3 = repmat( sign(mesh.t2e(:,3)), 1, 3 );
            
            at2e = abs(mesh.t2e);
            EE1 = f1 .* ( X(E(at2e(:,1),1),:) - X(E(at2e(:,1),2),:) );
            EE2 = f2 .* ( X(E(at2e(:,2),1),:) - X(E(at2e(:,2),2),:) );
            EE3 = f3 .* ( X(E(at2e(:,3),1),:) - X(E(at2e(:,3),2),:) );
            
            RE1 = - mesh.R * EE1(:);
            RE2 = - mesh.R * EE2(:);
            RE3 = - mesh.R * EE3(:);
            
            S = [RE1(:) RE2(:) RE3(:)]'; SS = S(:);
            
            GE = sparse(II,JJ,SS,3*mesh.nf,mesh.ne);
            ITA = spdiags(repmat(1./mesh.ta,3,1),0,3*mesh.nf,3*mesh.nf);

            GE = ITA*GE;
        end
        
        function [ DE ] = dive( mesh )
            IEA = spdiags(1./mesh.ea,0,mesh.ne,mesh.ne);
            TAC = spdiags(repmat(mesh.ta,3,1),0,3*mesh.nf,3*mesh.nf);

            DE = - IEA * mesh.GE' * TAC;
        end
        
        function [ op ] = fvf( mesh, vf )
            % op = [v]_{\cdot}^T * grad \in nf x nv
            vf = reshape(vf,[],3);
            
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            V = sparse(I,J,vf,mesh.nf,3*mesh.nf);

            op = V * mesh.G;
        end
        
        function [ op ] = ifvf( mesh, vf )
            % op = I^F_V * [v]_{\cdot}^T * grad \in nv x nv
            
            op = mesh.If2v * mesh.fvf( vf );
        end
        
        function [ op ] = fvfbar( mesh, f )
            % op = [grad f]_{\cdot}^T \in nf x 3nf
            
            gf = mesh.G * f;            
            GF = sparse(mesh.VFI,mesh.VFJ,gf,mesh.nf,3*mesh.nf);
            
            op = GF;
        end
        
        function [ op ] = ifvfbar( mesh, f )
                        
            op = mesh.If2v * mesh.fvfbar( f );
        end
        
        function [ op ] = v2f( mesh )
%             ff = mean( fv(mesh.triangles), 2 );
            
            VA = spdiags(mesh.va,0,mesh.nv,mesh.nv);
            ITA = spdiags(1./mesh.ta,0,mesh.nf,mesh.nf);
            
            op = ITA * mesh.If2v' * VA;
        end
    
        function [ op ] = f2v( mesh )

            IVA = spdiags(1./mesh.va,0,mesh.nv,mesh.nv);
            
            I = double( mesh.triangles ); I = I(:);
            J = repmat(1:mesh.nf,3,1)'; J = J(:);
            S = repmat(mesh.ta/3,3,1);
                        
            op = IVA*sparse(I,J,S,mesh.nv,mesh.nf);
        end
        
        function [ op ] = e2f( mesh )
            EA = spdiags(mesh.ea,0,mesh.ne,mesh.ne);
            ITA = spdiags(1./mesh.ta,0,mesh.nf,mesh.nf);
            
            op = ITA * mesh.If2e' * EA;
        end
        
        function [ op ] = f2e( mesh )
            IEA = spdiags(1./mesh.ea,0,mesh.ne,mesh.ne);
            
            I = double( abs(mesh.t2e) ); I = I(:);
            J = repmat(1:mesh.nf,3,1)'; J = J(:);
            S = repmat(mesh.ta/3,3,1);
                        
            op = IEA*sparse(I,J,S,mesh.ne,mesh.nf);
        end
        
        function [ pvf ] = project_vf( mesh, vf )
            pvf = vf - mesh.Nf .* repmat(dot(mesh.Nf,vf,2),1,3);
        end
        
        function [ NE1, NE2, EB, EBI ] = edge_basis( mesh )
            
            NE1 = MESH.normalize_vf( mesh.E1 );
            NE2 = reshape( mesh.R * NE1(:), [], 3 );
            
            B1 = sparse(mesh.VFI,mesh.VFJ,NE1,mesh.nf,3*mesh.nf);
            B2 = sparse(mesh.VFI,mesh.VFJ,NE2,mesh.nf,3*mesh.nf);
            
            EB = [B1; B2];
            EBI = EB';
            
            t1 = abs(mesh.e2t(mesh.inner_edges,1)); 
            t2 = abs(mesh.e2t(mesh.inner_edges,2));            
            
            B1t1 = sparse(mesh.VEI,mesh.VEJ,NE1(t1,:),mesh.nie,3*mesh.nie);
            B2t1 = sparse(mesh.VEI,mesh.VEJ,NE2(t1,:),mesh.nie,3*mesh.nie);            
            mesh.EBt1 = [B1t1; B2t1];

            B1t2 = sparse(mesh.VEI,mesh.VEJ,NE1(t2,:),mesh.nie,3*mesh.nie);
            B2t2 = sparse(mesh.VEI,mesh.VEJ,NE2(t2,:),mesh.nie,3*mesh.nie);            
            mesh.EBt2 = [B1t2; B2t2];            
        end
        
        function [ F1, F2 ] = smooth_frame( mesh, opt )
            
            if opt == 1
                tic; %fprintf('\tEIGS of GODF Laplacian: ');
                op = (mesh.VL+mesh.VL')/2;
                [f,~] = eigs(op,[],6,'SM');
                %fprintf('%f sec\n',toc);
                
                % Make sure mesh is single copmonent if using 
                % spectral smooth frame
                [~,dd]= eigs(-mesh.D*mesh.G,[],6,'SM');
                assert(dd(2,2) > 1e-10,...
                    'Mesh must be single component for a smooth spectral frame')
                
                F1 = MESH.normalize_vf( reshape( mesh.EBI*f(:,1), [], 3 ) );
            else
                F1 = mesh.F1;
            end
            
            F2 = reshape(mesh.R * F1(:), [], 3 );
            
            B1 = sparse(mesh.VFI,mesh.VFJ,F1,mesh.nf,3*mesh.nf);
            B2 = sparse(mesh.VFI,mesh.VFJ,F2,mesh.nf,3*mesh.nf);
            
            mesh.VLh = mesh.VLh * mesh.EB;
            mesh.EB = [B1; B2]; mesh.EBI = mesh.EB';
            mesh.VLh = mesh.VLh * mesh.EBI;

            mesh.VL = mesh.VLh' * mesh.VLh;
            
            t1 = abs(mesh.e2t(mesh.inner_edges,1)); 
            t2 = abs(mesh.e2t(mesh.inner_edges,2));            
            
            B1t1 = sparse(mesh.VEI,mesh.VEJ,F1(t1,:),mesh.nie,3*mesh.nie);
            B2t1 = sparse(mesh.VEI,mesh.VEJ,F2(t1,:),mesh.nie,3*mesh.nie);            
            mesh.EBt1 = [B1t1; B2t1];

            B1t2 = sparse(mesh.VEI,mesh.VEJ,F1(t2,:),mesh.nie,3*mesh.nie);
            B2t2 = sparse(mesh.VEI,mesh.VEJ,F2(t2,:),mesh.nie,3*mesh.nie);            
            mesh.EBt2 = [B1t2; B2t2];
        
        end
        
        function [ op ] = hodge( mesh )
            op = mesh.D' * mesh.VA * mesh.D + mesh.C' * mesh.EA * mesh.C;
            op = mesh.ITA * op;
            op = mesh.EB * op * mesh.EBI;
        end
        
        % op = oph' * oph so can be used with lsqnonlin
        function [ op, oph ] = godf( mesh, n )
            
            M = mesh;
            
            inner = M.inner_edges;
            
            t1 = abs(M.e2t(inner,1)); t2 = abs(M.e2t(inner,2));
            oe = -ones(M.nie,1); ze = zeros(M.nie,1);
            
            % Connection-specific computations
            EV = M.vertices(M.edges(inner,2),:) - M.vertices(M.edges(inner,1),:);
            EV = MESH.normalize_vf( EV );
            
            IN1 = atan2(dot(EV,M.F2(t1,:),2),dot(EV,M.F1(t1,:),2));
            IN2 = atan2(dot(EV,M.F2(t2,:),2),dot(EV,M.F1(t2,:),2));
            PT = n*(IN2-IN1);
            
            II = repmat((1:M.nie)',2,1); II = repmat([II; II+M.nie],2,1);
            JJ = [t1; t1+M.nf; t1; t1+M.nf; ...
                  t2; t2+M.nf; t2; t2+M.nf];
            SS = [cos(PT); -sin(PT); sin(PT); cos(PT); oe; ze; ze; oe]; 
            CovD = sparse(II,JJ,SS,2*M.nie,2*M.nf);
            
            Ws = spdiags(repmat(sqrt(M.ea(inner)),2,1), 0, 2*M.nie, 2*M.nie );
            oph = Ws * CovD;
            op = oph' * oph;

        end
        
        function f = smooth_delta( mesh, vid, r )
            d = geodesics_heat( mesh, vid );
            f = exp( -d/r ); 
        end
        function f = smooth_deltas( mesh, vids, r )
            d = geodesics_heats( mesh, vids );
            f = exp( -d/r ); 
        end
        
        function [ f, g ] = hodge_decomposition( mesh, u )
            u = reshape(u,[],1);
            
            L = mesh.D * mesh.G;
            divu = mesh.D * u;
            
            f = L \ divu;
            Jgradg = u - mesh.G*f;
            g = - L \ (mesh.D*(mesh.R*Jgradg));
        end
        
        function s = bbd( mesh )
            s = norm(max(mesh.vertices) - min(mesh.vertices));
        end
        
        function [ rvf ] = rotate_vf( mesh, vf )
            rvf = cross( mesh.Nf, vf );
        end
        
        % Converts MESH structures to libdirectional structures
        % ori indicates if the libdirectional dual edge is oriented in 
        % the same direction as the MESH edge (+1 if same, -1 if different)        
        % EF(i,1) is always on the left of [EV(i,1), EV(i,2)].
        % e2t(i,1) and e2t(i,2) are arbitrary, where the sign indicates
        % if their orientation agrees with edges(i,:).
        function [EV, EF, FE, ori] = tolibdirectional(mesh)
            ne = mesh.ne; nf = mesh.nf;
            e2t = mesh.e2t;
            
            EV = mesh.edges;
            EF = zeros(size(e2t,1),2);
            FE = zeros(nf,3);
            
            ori = ones(ne,1);
            for i=1:ne
                for j=1:2
                    if (e2t(i,j) < 0)
                        EF(i,2) = -e2t(i,j);
                        FE(EF(i,2),e2t(i,2+j)) = i;
                        if j == 1, ori(i) = -1; end;
                    else
                        if (e2t(i,j) > 0)
                            EF(i,1) = e2t(i,j);
                            FE(EF(i,1),e2t(i,2+j)) = i;
                            if j == 2, ori(i) = -1; end;
                        end
                    end
                end
            end
            FE = [FE(:,3),FE(:,1:2)];                        
        end
        
        % Converts edge based info from libdirectional structures to 
        % edge based info with MESH structures
        % ori indicates if the libdirectional dual edge is oriented in 
        % the same direction as the MESH edge (+1 if same, -1 if different)        
        % EF(i,1) is always on the left of [EV(i,1), EV(i,2)].
        % e2t(i,1) and e2t(i,2) are arbitrary, where the sign indicates
        % if their orientation agrees with edges(i,:).
        function [Xo, ori] = fromlibdirectional(mesh, EV, EF, FE, Xi)
            ne = mesh.ne; nf = mesh.nf;
            assert(size(EV,1) == ne);
            assert(size(EF,1) == ne);
            assert(size(FE,1) == nf);
            
            Xo = zeros(ne,1);
            ori = ones(ne,1);
            for i = 1:ne
                em = mesh.v2e(EV(i,1), EV(i,2));
                f1 = EF(i,1); f2 = EF(i,2);
                % not boundary
                if f1 > 0 && f2 > 0
                    f1m = mesh.e2t(em,1);
                    f2m = mesh.e2t(em,2);
                    assert(f1 == abs(f1m) || f1 == abs(f2m));
                    assert(f2 == abs(f1m) || f2 == abs(f2m));
                    if f1 == abs(f1m) 
                        fm = f1m;
                    elseif f1 == abs(f2m) 
                        fm = f2m;
                        ori(em) = -1;
                    else
                        error('?');
                    end       
                % Boundary
                else
                    f = max(f1,f2);
                    fm = mesh.e2t(em,1);
                    assert(f == abs(fm));
                    assert(mesh.e2t(em,2) == 0);
                end
                Xo(em) = Xi(i);
            end                        
        end        
    end
    
    methods (Static)
        
        function [ NV ] = NORMV( V, dim, G )
            if nargin < 2
                dim = 3;
            end
            if nargin < 3
                G = speye(size(V,1));
            end
            
            sf = size(V,1)/dim;
            
            NV = V.*(G*V);
            if dim == 3
                NV = sqrt( NV(1:sf,:) + NV(sf+1:2*sf,:) + NV(2*sf+1:end,:) ); 
            elseif dim == 2
                NV = sqrt( NV(1:sf,:) + NV(sf+1:2*sf,:) ); 
            else
                error('?');
            end
            assert(all(NV>1e-14));
        end

        function [ nv ] = normv( v, dim )
            if nargin > 1
                v = reshape(v,[],dim);
            end
            nv = sqrt(sum(v.^2,2));
        end 
        
        function [ nnv ] = normalize_vf( v, dim )
            if nargin > 1
                v = reshape(v, [], dim);
            else
                dim = 3;
            end
            nv = MESH.normv(v); nnv = v; fid = nv > 1e-10;
            nnv(fid,:) = nnv(fid,:) ./ repmat(nv(fid),1,dim);
            if nargin > 1
                nnv = nnv(:);
            end
        end

        % Volume of multi-vector - cross product of every two adjacent
        % fields. Input is 2D (in basis of tangent plane per face)
        function vol = mvf_vol_2d(mvf, nf)
            if size(mvf,2) == 1
                mvf = reshape(mvf,[],nf)';
            end
            
            if size(mvf,1) == nf
                assert(mod(size(mvf,2), 2) == 0);
                n = round(size(mvf,2) / 2);
                mvf = reshape(mvf, [], n);
            else
                n = size(mvf,2);
            end
            
            vol = zeros(nf,n);
            Nf = repmat([0,0,1],nf,1);
            for i=1:n
                ip1 = i+1; if ip1 > n, ip1 = 1; end;
                v1 = reshape(mvf(:,i),[],2);
                v2 = reshape(mvf(:,ip1),[],2);
                v1 = [v1,zeros(nf,1)];
                v2 = [v2,zeros(nf,1)];
                vol(:,i) = dot(cross(v1,v2,2),Nf,2);
            end
        end        
        
        % Volume of multi-vector - cross product of every two adjacent
        % fields. Input is 3D.
        function vol = mvf_vol(mvf, Nf)
            nf = size(Nf,1);
            if size(mvf,2) == 1
                mvf = reshape(mvf,[],nf)';
            end
            
            if size(mvf,1) == nf
                assert(mod(size(mvf,2),3) == 0);
                n = round(size(mvf,2) / 3);
                mvf = reshape(mvf, [], n);
            else
                n = size(mvf,2);
            end
            
            vol = zeros(nf,n);
            for i=1:n
                ip1 = i+1; if ip1 > n, ip1 = 1; end;
                v1 = reshape(mvf(:,i),[],3);
                v2 = reshape(mvf(:,ip1),[],3);
                vol(:,i) = dot(cross(v1,v2,2),Nf,2);
            end
        end
        
        function F = normalize_func(min_new,max_new,f)
            n = size(f,2);
            
            F = zeros(size(f));
            for i = 1:n
                fnew = f(:,i) - min(f(:,i));
                F(:,i) = (max_new-min_new)*fnew/max(fnew) + min_new;
            end
        end
        
        function [ r ] = rmse( f1, f2 )
            mf2 = mean(f2);
            r = sum((f1-f2).^2) / sum((f2-mf2).^2);
        end

        function normalize_mesh(inname, outname, bbd_out)
            if nargin < 3
                bbd_out = 1;
            end
            m = MESH(inname);
            me = mean(m.vertices);
            xx = m.vertices - me;
            MM = max(m.vertices); mm = min(m.vertices);
            bbd = norm(MM - mm);
            xx = xx / bbd * bbd_out;
            MESH_IO.woff([outname '.off'], xx, m.triangles);
        end

        
    end
end
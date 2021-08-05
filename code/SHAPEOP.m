classdef SHAPEOP
        
    properties (Constant)
        shapeop_func = @SHAPEOP.shape_operator_FTF;
    end
    
    methods (Static)
        
        function [ Nv ] = vertex_normals( mesh )
            
            TA = spdiags(mesh.ta,0,mesh.nf,mesh.nf);
            
            I = double( repmat(mesh.triangles(:),3,1) );
            J = repmat(1:3,3*mesh.nf,1); J = J(:);
            S = repmat(TA*mesh.Nf,3,1); S = S(:);
            
            Nv = full( sparse(I,J,S,mesh.nv,3) );
            Nv = MESH.normalize_vf( Nv );
        end
        
        function [ V, D ] = shape_operator( mesh )
            [V,D] = SHAPEOP.shape_operator_R( mesh );
        end
        
        % similar to Rusinkiewicz, 2004
        % "Estimating curvatures and their derivatives on triangle meshes"
        % V will be in the abs(kmin) direction
        function [ V, D  ] = shape_operator_R( mesh )
            
            HITAR = .5 * sqrt( mesh.ITA ) * mesh.R;
            
            Nv = SHAPEOP.vertex_normals( mesh );

            gE = cell(3,1);
            gE{1} = reshape(HITAR * mesh.E1(:),[],3);
            gE{2} = reshape(HITAR * mesh.E2(:),[],3);
            gE{3} = reshape(HITAR * mesh.E3(:),[],3);

            V = zeros(mesh.nf,3);
            D = zeros(mesh.nf,1);
            for j = 1:mesh.nf
                s = 0;
                for i = 1:3
                    s = s + Nv(mesh.triangles(j,i),:)' * gE{i}(j,:);
                end
                
                % project into the tangent plane
                p = eye(3) - mesh.Nf(j,:)'*mesh.Nf(j,:);
                s = p * s;
                
                % symmetrize
                s = (s+s')/2;
                
                % compute eigenvectors
                [v,d] = eig( s ); dd = diag(d);
                
                % eliminate normal eigenvector
                e = [mesh.E1(j,:)' mesh.E2(j,:)'];
                [~,i] = min( MESH.normv( v'*e ) );
                v(:,i) = []; dd(i) = [];
                
                [~,i] = min( dd );
                V(j,:) = v(:,i)';
                D(j) = ( dd(1) - dd(2) ).^2;
            end
        end
        
        function [dminf, dmaxf, kminf, kmaxf] = shape_operator_R_all( mesh )
            
            HITAR = .5 * sqrt( mesh.ITA ) * mesh.R;
            
            Nv = SHAPEOP.vertex_normals( mesh );

            gE = cell(3,1);
            gE{1} = reshape(HITAR * mesh.E1(:),[],3);
            gE{2} = reshape(HITAR * mesh.E2(:),[],3);
            gE{3} = reshape(HITAR * mesh.E3(:),[],3);

            dminf = zeros(mesh.nf,3); dmaxf = dminf;
            kminf = zeros(mesh.nf,1); kmaxf = kminf;
            for j = 1:mesh.nf
                s = 0;
                for i = 1:3
                    s = s + Nv(mesh.triangles(j,i),:)' * gE{i}(j,:);
                end
                
                % project into the tangent plane
                p = eye(3) - mesh.Nf(j,:)'*mesh.Nf(j,:);
                s = p * s;
                
                % symmetrize
                s = (s+s')/2;
                
                % compute eigenvectors
                [v,d] = eig( s ); dd = diag(d);
                
                % eliminate normal eigenvector
                e = [mesh.E1(j,:)' mesh.E2(j,:)'];
                [~,i] = min( MESH.normv( v'*e ) );
                v(:,i) = []; dd(i) = [];
                
                [~,i] = min( dd );
                dminf(j,:) = v(:,i)';
                dmaxf(j,:) = v(:,3-i)';
                kminf(j) = dd(i);
                kmaxf(j) = dd(3-i);
            end
        end
        
        % based on Cohen-Steiner and Morvan, 2003
        % Seems to be the only one that gives correct results on a sphere
        function [dminf, dmaxf, kminf, kmaxf] = shape_operator_CSM( mesh )
            X = mesh.vertices; T = double( mesh.triangles );
            
            % adjacency matrix
            I = [T(:,2);T(:,3);T(:,1)];
            J = [T(:,3);T(:,1);T(:,2)];
            S = [1:mesh.nf,1:mesh.nf,1:mesh.nf];
            E = sparse(I,J,S,mesh.nv,mesh.nv);

            % vertex-based shape operator
            SX = zeros(9,mesh.nv);
            for i = 1:mesh.nv
                [~,n,~] = find( E(i,:) );
                sx = 0;
                for j = 1:numel(n)
                    e = X(i,:) - X(n(j),:);
                    
                    ei = mesh.v2e(min(i,n(j)),max(i,n(j)));
                    ni = mesh.Nf(abs( mesh.e2t( ei, 1 ) ),:);
                    nj = mesh.Nf(abs( mesh.e2t( ei, 2 ) ),:);
                    beta = acos( dot(ni,nj) );
                    
                    sx = sx + beta*(e'*e)/norm(e);
                end
                sx = .5*sx/mesh.va(i);
                SX(:,i) = sx(:);
            end
            
            dminf = zeros(mesh.nf,3); dmaxf = dminf;
            kminf = zeros(mesh.nf,1); kmaxf = kminf;
            
            % face-based shape operator
            for j = 1:mesh.nf
                st = reshape( sum( SX(:,T(j,:)), 2 ) / 3, [], 3 );
                
                % project into the tangent plane
                p = eye(3) - mesh.Nf(j,:)'*mesh.Nf(j,:);
                st = p * st;
                
                % symmetrize
                st = (st+st')/2;
                
                % compute eigenvectors
                [v,d] = eig( st ); v = real(v); d = real(d);
                dd = diag(d);
                
                % eliminate normal eigenvector
                e = [mesh.E1(j,:)' mesh.E2(j,:)'];
                [~,i] = min( MESH.normv( v'*e ) );
                v(:,i) = []; dd(i) = [];
                
                [~,i] = min( dd );
                dminf(j,:) = v(:,i)';
                dmaxf(j,:) = v(:,3-i)';
                kminf(j) = dd(i);
                kmaxf(j) = dd(3-i);
            end
        end
        
        % Operator used for FTF, gradient of vertex normals
        function [dminf, dmaxf, kminf, kmaxf] = shape_operator_FTF( mesh )
            T = mesh.triangles;
            ITA = repmat(1./mesh.ta,1,3);

            gE = cell(3,1);
            gE{1} = .5 * ITA .* mesh.rotate_vf( mesh.E1 );
            gE{2} = .5 * ITA .* mesh.rotate_vf( mesh.E2 );
            gE{3} = .5 * ITA .* mesh.rotate_vf( mesh.E3 );

            dminf = zeros(mesh.nf,3); dmaxf = dminf;
            kminf = zeros(mesh.nf,1); kmaxf = kminf;

            for j = 1:mesh.nf
                s = 0;
                for i = 1:3
                    s = s + mesh.Nv(T(j,i),:)'*gE{i}(j,:);
                end
                % project into the tangent plane
                p = eye(3) - mesh.Nf(j,:)'*mesh.Nf(j,:);
                s = p * s;
                % symmetrize
                s = (s+s')/2;

                % compute eigenvectors
                [v,d] = eig( s ); v = real(v); d = real(d);
                dd = diag(d);
                [abs_dd, i] = sort(abs(dd));
                
                if all(abs_dd>1e-15) || (abs_dd(1)<1e-15 && all(abs_dd(2:3)>1e-15)) % all ev ~= 0 or 1 ev = 0
                    % eliminate normal eigenvector
                    e = [mesh.E1(j,:)' mesh.E2(j,:)'];
                    [~,i] = min( MESH.normv( v'*e ) );
                    v(:,i) = []; dd(i) = [];
                    
                    [~,i] = sort( dd );
                    v = v(:,i); dd = dd(i);
                    
                elseif all(abs_dd<1e-15) % all ev = 0
                    v1 = mesh.normalize_vf(mesh.E1(j,:));
                    v2 = cross(v1, mesh.Nf(j,:));
                    v2 = v2./sqrt(sum(v2.^2,2));
                    v = [v1', v2']; dd = [0, 0];
                    
                elseif all(abs_dd(1:2)<1e-15) && (abs_dd(3)>1e-15) % 2 ev = 0
                    v = v(:,i); v1 = v(:,3); dd1 = dd(3);
                    assert(MESH.normv(v1'*mesh.Nf(j,:)')<1e-15);
                    v2 = cross(v1', mesh.Nf(j,:));
                    v2 = v2./sqrt(sum(v2.^2,2));
                    v = [v1, v2']; dd = [dd1, 0];
                    [~,i] = sort( dd );
                    v = v(:,i); dd = dd(i);
                end
                
                dminf(j,:) = v(:,1)';
                dmaxf(j,:) = v(:,2)';
                kminf(j) = dd(1);
                kmaxf(j) = dd(2);

            end
            
        end

        % return a 2fx2f matrix, where each 2x2 matrix is the shape operator 
        function SO = shapeop(mesh)
            nf = mesh.nf;
            [dminf, dmaxf, kminf, kmaxf] = SHAPEOP.shapeop_func(mesh);
            dminfb = reshape(mesh.EB*dminf(:),[],2);
            dmaxfb = reshape(mesh.EB*dmaxf(:),[],2);
            v11 = dminfb(:,1); v12 = dminfb(:,2);
            v21 = dmaxfb(:,1); v22 = dmaxfb(:,2);
            v = [spdiags(v11,0,nf,nf), spdiags(v12,0,nf,nf);...
                 spdiags(v21,0,nf,nf), spdiags(v22,0,nf,nf)];
            d11 = kminf; d12 = kmaxf;
            d = [spdiags(d11,0,nf,nf), sparse(nf,nf); ...
                 sparse(nf,nf), spdiags(d12,0,nf,nf)];
            SO = v'*d*v;            
        end
        
        % return a 2fx2f matrix, where each 2x2 matrix is such that
        % <v,v>_{g_i} = v'*g_i*v for a single face i and a vector v \in R^2
        % in the basis of that face. In addition we have S'*g*S = Id, where
        % S is a 2x2 matrix whose columns are the curvature directions
        % dmin, dmax, scaled by sqrt(r)/sqrt(kmins), sqrt(r)/sqrt(kmaxs), 
        % where kmins = sqrt(kmin^2 + e^2), and kmaxs = sqrt(kmax^2 + e^2)
        % respectively, and r 
        % 1/e will be the max element size
        % S - struct with dmin, dmax, kmin, kmax
        function [g, S] = curvature_sizing(mesh, e, r)
            if nargin < 2
                e = 1e-6;
            end
            nf = mesh.nf;
            [dminf, dmaxf, kminf, kmaxf] = SHAPEOP.shapeop_func(mesh);
            dminfb = reshape(mesh.EB*dminf(:),[],2);
            dmaxfb = reshape(mesh.EB*dmaxf(:),[],2);
            kminfs = sqrt(kminf.^2 + e^2);
            kmaxfs = sqrt(kmaxf.^2 + e^2);
            v11 = dminfb(:,1); v12 = dminfb(:,2);
            v21 = dmaxfb(:,1); v22 = dmaxfb(:,2);
            v = [spdiags(v11,0,nf,nf), spdiags(v12,0,nf,nf);...
                 spdiags(v21,0,nf,nf), spdiags(v22,0,nf,nf)];
            d11 = kminfs/r; d12 = kmaxfs/r;

            d = [spdiags(d11,0,nf,nf), sparse(nf,nf); ...
                 sparse(nf,nf), spdiags(d12,0,nf,nf)];
            g = v'*d*v;
            S.dminf = dminf; S.dmaxf = dmaxf; S.kminf = kminf; S.kmaxf = kmaxf;
            S.d11 = d11; S.d12 = d12;
            S.kminfs = kminfs; S.kmaxfs = kmaxfs;
        end
        
        % Unit testing for curvature_sizing
        function check_curvature_sizing(mesh, e, r)
            [g,S] = SHAPEOP.curvature_sizing(mesh, e, r);            
            dminf = S.dminf; dmaxf = S.dmaxf;
            kminf = S.kminf; kmaxf = S.kmaxf;

            dminfb = reshape(mesh.EB*dminf(:),[],2);
            dmaxfb = reshape(mesh.EB*dmaxf(:),[],2);
            kminfs = sqrt(kminf.^2 + e^2);
            kmaxfs = sqrt(kmaxf.^2 + e^2);
            
            % Scale 
            dminfb = spdiags(sqrt(r./kminfs), 0, mesh.nf, mesh.nf)*dminfb;
            dmaxfb = spdiags(sqrt(r./kmaxfs), 0, mesh.nf, mesh.nf)*dmaxfb;
            
            % Check identity
            dminfb = dminfb(:); dmaxfb = dmaxfb(:);
            sf = mesh.nf;
            NV = dminfb.*(g*dminfb); 
            NV = NV(1:sf,:) + NV(sf+1:2*sf,:); 
            assert(max(abs(NV - 1)) < 1e-3);

            NV = dminfb.*(g*dmaxfb); 
            NV = NV(1:sf,:) + NV(sf+1:2*sf,:); 
            assert(max(abs(NV)) < 1e-3);

            NV = dmaxfb.*(g*dminfb); 
            NV = NV(1:sf,:) + NV(sf+1:2*sf,:); 
            assert(max(abs(NV)) < 1e-3);

            NV = dmaxfb.*(g*dmaxfb); 
            NV = NV(1:sf,:) + NV(sf+1:2*sf,:); 
            assert(max(abs(NV - 1)) < 1e-3);
        end
        
        % Generates a guiding field that is (mostly) aligned with the
        % minimum absolute curvature direction
        % Where cons2 is not zero, the U direction should be aligned with
        % either cons2 or -cons2
        % Where cons4 is not zero, the U or V directions should be aligned with
        % cons4 or -cons4
        function [cons2, cons4, lost_ind, damin, dmin, RHO, PHI,...
            locse, locsps, locsh] = guiding_field(mesh)
            rho_min = 0.001;         % min confidence for non-planar
            rho_minp = 0.05;         % min super confidence for parabolic non planar
            delta_phip = 0.01;       % distance from pi/4 for parabolic
            delta_phie = 0.2;        % distance from pi/2 for umbilic

            
            m = mesh;
            [dmin, dmax, kmin, kmax] = SHAPEOP.shapeop_func(m);
            
            PHI = abs(atan((kmin + kmax)./(kmin - kmax)));
            RHO = sqrt(kmin.*kmin + kmax.*kmax);
            RHO = RHO./max(RHO);

            dmin = MESH.normalize_vf(dmin);
            dmax = MESH.normalize_vf(dmax);

            %%
            nf = m.nf;
            [~,iik] = sort(abs([kmin, kmax]),2);
            kamin = zeros(nf,1); kamax = kamin;
            damin = zeros(nf,3); damax = damin;
            for i=1:nf
                if iik(i,1) == 1
                    kamin(i) = kmin(i);
                    kamax(i) = kmax(i);

                    damin(i,:) = dmin(i,:);
                    damax(i,:) = dmax(i,:);
                else
                    kamin(i) = kmax(i);
                    kamax(i) = kmin(i);

                    damin(i,:) = dmax(i,:);
                    damax(i,:) = dmin(i,:);
                end
            end

            % compute connected components of constrained
            op = m.godf(1);
            G = graph((abs(op(1:m.nf,1:m.nf))>0));
            locsps = find((RHO > rho_minp) .* (abs(PHI - pi/4) < delta_phip));       % green
            locsh = find((RHO > rho_min) .* (PHI < pi/4 - delta_phip));             % blue
            locspq = find((RHO > rho_min & RHO < rho_minp).*(abs(PHI - pi/4) < delta_phip));   % cyan
            locse = find((RHO > rho_min) .* (PHI < pi/2 - delta_phie) .* (PHI > pi/4 + delta_phip));    % red
            locs = union(locsps, union(locsh,union(locspq, locse)));
            H = subgraph(G, locs);
            bins = conncomp(H);
            bbins = zeros(m.nf,1);
            bbins(locs) = bins;
            
            % For each CC decide on either dmin or dmax by majority vote of parabolic
            % regions
            constraints = zeros(nf,3);
            locs_voted = [];
            for i=1:max(bbins)
                locsi = find(bbins == i); 
                locsip = intersect(locsi, locsps);
                if length(locsip) > 1
                    ak = m.ta.*abs(kamax);
                    choicei = ...
                        (sum(ak(locsip).*(iik(locsip,1)==1)) > sum(ak(locsip).*(iik(locsip,1)==2)));
                    if choicei == 1
                        constraints(locsi,:) = dmin(locsi,:);
                    else
                        constraints(locsi,:) = dmax(locsi,:);
                    end
                    locs_voted = union(locs_voted, locsip);
                end
            end

            lost_ind = zeros(m.nf,1);
            lost_ind(locs_voted,:) = (MESH.normv(constraints(locs_voted,:)-damin(locs_voted,:)) > 1e-10);

            cons_2_locs = find(MESH.normv(constraints) > 0);
            cons_4_locs = setdiff(locs, cons_2_locs);

            cons2 = constraints; 
            cons4 = zeros(nf,3); cons4(cons_4_locs,:) = damin(cons_4_locs,:);            
        end

        % Generates a guiding field that is (mostly) aligned with the
        % minimum absolute curvature direction
        % Where cons2 is not zero, the U direction should be aligned with
        % either cons2 or -cons2
        % Where cons4 is not zero, the U or V directions should be aligned with
        % cons4 or -cons4
        % No connected components 
        function [cons2, cons4, lost_ind, damin, dmin, RHO, PHI, ...
                locse, locsp, locsh] = ...
                guiding_field2(mesh)
            % Note that both RHO and PHI are scale invariant
            rho_min = 0.01;         % min confidence for non-planar
            delta_phip = 0.05;       % distance from pi/4 for parabolic
            delta_phie = 0.2;        % distance from pi/2 for umbilic

            
            m = mesh;
            [dmin, dmax, kmin, kmax] = SHAPEOP.shapeop_func(m);
            
            PHI = abs(atan((kmin + kmax)./(kmin - kmax)));
            RHO = sqrt(kmin.*kmin + kmax.*kmax);
            RHO = RHO./max(RHO);

            dmin = MESH.normalize_vf(dmin);
            dmax = MESH.normalize_vf(dmax);

            %%
            nf = m.nf;
            [~,iik] = sort(abs([kmin, kmax]),2);
            kamin = zeros(nf,1); kamax = kamin;
            damin = zeros(nf,3); damax = damin;
            for i=1:nf
                if iik(i,1) == 1
                    kamin(i) = kmin(i);
                    kamax(i) = kmax(i);

                    damin(i,:) = dmin(i,:);
                    damax(i,:) = dmax(i,:);
                else
                    kamin(i) = kmax(i);
                    kamax(i) = kmin(i);

                    damin(i,:) = dmax(i,:);
                    damax(i,:) = dmin(i,:);
                end
            end
            
            % Elliptic without umbilics
            locse = find((RHO > rho_min) ...            % Not planar
                .* (PHI < pi/2 - delta_phie) ...        % Not umbilic
                .* (PHI > pi/4 + delta_phip));          % Eliptic

            % Hyperbolic 
            locsh = find((RHO > rho_min) ...            % Not planar
                .* (PHI < pi/4 - delta_phip));          % Hyperbolic
            
            % parabolic
            locsp = find((RHO > rho_min) ...            % Not planar
                .*(abs(PHI - pi/4) < delta_phip));      % Parabolic
            
            
            cons2 = zeros(m.nf,3); cons4 = cons2;

            % N=2 constraints: parabolic, damin
            cons2(locsp,:) = damin(locsp,:);
            
            % N=4, constraints: hyperbolic non umbilic and elliptic non
            % umbilic, damin
            cons4(locsh,:) = damin(locsh,:);
            cons4(locse,:) = damin(locse,:);
                        
            % remove constraints on faces which have boundary vertices,
            % curvature is not reliable there
            cons2(mesh.bf==1,:) = 0;
            cons4(mesh.bf==1,:) = 0;
            
            % To be compatible with prev version
            lost_ind = [];
        end
                
        % Porting of the libigl implementation
        function [dmin, dmax, kmin, kmax] = shape_operator_IGL(mesh)
            radius = 5;
            useKring = true;
            type = 'FACE';
            
            assert(radius > 2, 'shapeop_igl: The radius has to be bigger than 2!')
            
            if strcmp(type, "VERTEX")
                dmin = zeros(mesh.nv, 3);
                dmax = zeros(mesh.nv, 3);
                kmin = zeros(mesh.nv, 1);
                kmax = zeros(mesh.nv, 1);
                
                cc = CURVATURE_CALCULATOR;
                cc.init(mesh);
                
                cc.sphereRadius = radius;
                if useKring
                    cc.st = 'K_RING_SEARCH';
                    cc.kRing = radius;
                end
                cc.computeCurvature();
                
                dmin = cc.curvDir(1: end, 1:3);
                dmin ./ sqrt(sum(dmin.^2,2));
                dmax = cc.curvDir(1: end, 4:6);
                dmax ./ sqrt(sum(dmax.^2,2));
                kmin = cc.curv(:, 1);
                kmax = cc.curv(:, 2);
                status = cc.curvatureComputed;
            elseif strcmp(type, "FACE")
                [VCBary, FCBary] = false_barycentric_subdivision(mesh.triangles, mesh.vertices);
                
                fakemesh = MESH( [mesh.name, '-fake']);%, VCBary, FCBary);
                cc = CURVATURE_CALCULATOR;
                cc.init(fakemesh);
                
                cc.sphereRadius = radius;
                if useKring
                    cc.st = 'K_RING_SEARCH';
                    cc.kRing = radius;
                end
                cc.computeCurvature();
                
                dmin = cc.curvDir(mesh.nv + 1: end, 1:3);
                dmin ./ sqrt(sum(dmin.^2,2));
                dmax = cc.curvDir(mesh.nv + 1: end, 4:6);
                dmax ./ sqrt(sum(dmax.^2,2));
                kmin = cc.curv(mesh.nv + 1: end, 1);
                kmax = cc.curv(mesh.nv + 1: end, 2);
                
                % some postprocessing
                for i  = 1 : mesh.nf
                    dmax(i,:) = dmax(i,:) -  (dot(dmax(i,:), mesh.Nf(i,:))) * mesh.Nf(i,:);
                    dmax(i,:) = dmax(i,:) ./ norm(dmax(i,:));
                    
                    dmin(i,:) = dmin(i,:) -  (dot(dmin(i,:), mesh.Nf(i,:))) * mesh.Nf(i,:);
                    dmin(i,:) = dmin(i,:) -  (dot(dmin(i,:), dmax(i,:))) * dmax(i,:);
                    dmin(i,:) = dmin(i,:) ./ norm(dmin(i,:));
                    
                    if dot(cross(dmin(i,:), dmax(i,:)), mesh.Nf(i,:)) < 0
                        dmin(i,:) = -dmin(i,:);
                    end
                end
                status = cc.curvatureComputed;
            else
                throw(MException('shapeop_igl:BadType', 'The type argument can be either VERTEX or FACE!'));
            end
        end
        
        % based on Cohen-Steiner and Morvan, 2003, operator is 2-ring smoothed
        % Amir Vaxman's implementation
        function [dminf, dmaxf, kminf, kmaxf]=shape_operator_AV(mesh)
            X = mesh.vertices;
            T = mesh.triangles;
            edges = mesh.edges;
            e2t = mesh.e2t;
            ieid = find(mesh.ie);
            ne = mesh.ne; nv = mesh.nv; nf = mesh.nf;
            Nf = mesh.Nf;
            
            ev = X(edges(:,2),:)-X(edges(:,1),:);
            
            %computing angles 
            cosa = zeros(ne,1);
            signa = zeros(ne,1);

            tposi = (e2t(ieid,1)>0) + 2*(e2t(ieid,2)>0);
            tnegi = 3-tposi;
            post_ind = sub2ind(size(e2t), ieid, tposi);
            negt_ind = sub2ind(size(e2t), ieid, tnegi);

            Nft1 = Nf(abs(e2t(ieid,1)),:);
            Nft2 = Nf(abs(e2t(ieid,2)),:);

            Nft_pos = Nf(abs(e2t(post_ind)),:);
            Nft_neg = Nf(abs(e2t(negt_ind)),:);

            cosa(ieid) = dot(Nft1,Nft2,2);
            signa(ieid) = sign(dot(cross(Nft_pos, Nft_neg), ev(ieid,:),2));

            ang = real(acos(cosa)).*signa; %pi-

            % Compute shape operator per edge
            Se = zeros(ne,9);           
            ev = MESH.normalize_vf(ev);
            for i = 1:ne
                Se(i,:) = reshape(0.5*ang(i)*ev(i,:)'*ev(i,:),1,9);
            end

            % Average from edges to vertices            
            EVAdj = sparse(double(repmat((1:ne)', 1,2)), double(edges), ...
                        ones(ne,2), ne, nv);
            Sv = EVAdj'*Se;
            
            % Smooth on 2-ring vertices
            OneRingSum = ...
                speye(nv) + ...
                sparse(double(T),double([T(:,2);T(:,3);T(:,1)]), ones(3*nf,1), nv,nv) + ...
                sparse(double([T(:,2);T(:,3);T(:,1)]), double(T),ones(3*nf,1), nv,nv);            
            TwoRingSum = OneRingSum * OneRingSum;
            [I,J] = find(TwoRingSum);
            TwoRingSum = sparse(I,J,ones(length(I),1), nv, nv);

            t2v = sparse(double(repmat((1:nf)',1,3)), ...
                double(reshape(T,nf*3,1)), ones(nf*3,1)/3, nf, nv);
            Arv = t2v' * mesh.ta;
                        
            TRVorArea = TwoRingSum * Arv;
            TRSv = TwoRingSum * Sv;
 
            Sv = TRSv ./ repmat(TRVorArea,1,9);

            dminf=zeros(nf,3);
            dmaxf=zeros(nf,3);
            kminf=zeros(nf,1);
            kmaxf=zeros(nf,1);
            Sp=(Sv(T(:,1),:)+Sv(T(:,2),:)+Sv(T(:,3),:))/3;
            for i = 1:nf
                CurrSv = reshape(Sp(i,:),3,3);
                [V,D] = eig(CurrSv);

                %weeding out the normal
                Eigs = real([diag(abs(D)) diag(D) V']);
                Eigs = sortrows(Eigs,1);
                Eigs = Eigs(2:3,2:end);
                Eigs = sortrows(Eigs,1);

                dminf(i,:) = Eigs(2,2:end);
                %removing normal component
                dminf(i,:) = dminf(i,:)-dot(dminf(i,:),Nf(i,:),2)*Nf(i,:);
                dminf(i,:) = dminf(i,:)./normv(dminf(i,:));
                dmaxf(i,:) = cross(Nf(i,:), dminf(i,:));                  
                kminf(i,:) = Eigs(1,1);
                kmaxf(i,:) = Eigs(2,1);
            end
        end
    end    
end


classdef MESH_VIS

    
    properties
    end
    
    methods (Static)
        
        function mesh(M,varargin)
            p = inputParser;
            p.KeepUnmatched=true;
            addOptional(p,'FaceColor','w');
            addOptional(p,'EdgeColor','k');
            addOptional(p,'CW',[]);
            addOptional(p,'UVt',[]);
            addOptional(p,'FUVt',[]);
            addOptional(p,'Texture','cross.png');
            addOptional(p,'SingularPts',0);
            addOptional(p,'SingularPtsR',20);
            addOptional(p,'FaceAlpha',.8);
            addOptional(p,'locs',{[1:M.nf]'});
            addOptional(p,'LineWidth',.5);
            addOptional(p,'docked',1);

            parse(p,varargin{:});
            
            locs = p.Results.locs{1};

            if isprop(M,'triangles')
                faces = M.triangles;
            elseif isprop(M,'quads')
                faces = M.quads;
            elseif isprop(M,'faces')
                faces = M.faces;
            else
                error('?');
            end
            
            if isempty( p.Results.UVt )
                if isempty(p.Results.CW) == 1
                    patch('faces',faces(locs,:),'vertices',M.vertices, ...
                        'FaceColor',p.Results.FaceColor, ...
                        'EdgeColor',p.Results.EdgeColor, ...
                        'FaceAlpha',p.Results.FaceAlpha, ...
                        'linewidth',p.Results.LineWidth);
                else
                    CW = 'CData';
                    if size(p.Results.CW,2) == 3; CW = 'FaceVertexCData'; end;

                    patch('faces',faces(locs,:),'vertices',M.vertices, ...
                        CW,p.Results.CW,...
                        'FaceColor','interp', ...
                        'EdgeColor',p.Results.EdgeColor, ...
                        'FaceAlpha',p.Results.FaceAlpha, ...
                        'linewidth',p.Results.LineWidth);
                end

                if p.Results.SingularPts == 1
                    MESH_VIS.singular_pts( M, p.Results.SingularPtsR );
                    %                         camlight; lighting phong; material dull;
                end
                    
            else
                I = imread(p.Results.Texture);
                
                patcht2(M.triangles,M.vertices,...
                    p.Results.FUVt,p.Results.UVt,I);
            end
            
            
            cameratoolbar; cameratoolbar('SetCoordSys','none');
            axis equal; axis off;
            if p.Results.docked
                set(gcf,'windowstyle','docked');
            end
        end
        
        function [ M2 ] = nc_mesh(M)
            M2.vertices = ( M.vertices(M.edges(:,1),:) + ...
                M.vertices(M.edges(:,2),:) )/2;
            M2.triangles = abs(M.t2e);
            
            M2.nf = M.nf;
            M2.nv = M.ne;
            
            M2.om = M;
            M2.om.vertices = M.vertices-M.Nv*M.mel/10;
        end
        
        function texture_mesh(M,varargin)
            p = inputParser;
            addOptional(p,'FaceColor','w');
            addOptional(p,'FaceAlpha',.5);
            addOptional(p,'EdgeColor','k');
            parse(p,varargin{:});
            
            patch('faces',M.texture_triangles,'vertices',M.texture_vertices, ...
                'FaceColor',p.Results.FaceColor, ...
                'EdgeColor',p.Results.EdgeColor, ...
                'FaceAlpha',p.Results.FaceAlpha);    
            
            axis equal; axis on
            set(gcf,'windowstyle','docked');            
        end
        
        function func(M,f,varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addOptional(p,'EdgeColor','none');
            addOptional(p,'Caxis','auto');
            addOptional(p,'View',[0 1 0]);
            addOptional(p,'Dock',1);
            addOptional(p,'Colormap','jet');
            addOptional(p,'OpenGL',0);
            addOptional(p,'LightPos',[0 0 1]);
            addOptional(p,'func_type',[]);
            parse(p,varargin{:});
            
            szf = size(f);
            if szf(1) == szf(2) && norm(f - diag(diag(f)),'fro') < 1e-8
                f = diag(f);
            end
            
            FC = 'w';
            if (isempty(p.Results.func_type) && szf(1) == M.nv) || ...
                    (strcmp(p.Results.func_type,'v') && szf(1) == M.nv)
                FC = 'interp';
            elseif (isempty(p.Results.func_type) && szf(1) == M.nf) || ...
                    (strcmp(p.Results.func_type,'f') && szf(1) == M.nf)
                FC = 'flat';
            elseif (isempty(p.Results.func_type) && szf(1) == M.ne) || ...
                    (strcmp(p.Results.func_type,'e') && szf(1) == M.ne)
                M = MESH_VIS.nc_mesh( M );
                patch('faces',M.om.triangles,'vertices',M.vertices, ...
                    'FaceColor',FC, ...
                    'EdgeColor',p.Results.EdgeColor);
                FC = 'interp';
            else
                error('bad function size');
            end
            
            CW = 'CData';
            if szf(2) == 3; CW = 'FaceVertexCData'; end;
            
            if any(strcmp(properties(M), 'triangles')) == 1
                patch('faces',M.triangles,'vertices',M.vertices, ...
                    CW,f, ...
                    'FaceColor',FC, ...
                    'EdgeColor',p.Results.EdgeColor, ...
                        'FaceLighting','phong');
            elseif any(strcmp(properties(M), 'quads')) == 1
                patch('faces',M.quads,'vertices',M.vertices, ...
                    CW,f, ...
                    'FaceColor',FC, ...
                    'EdgeColor',p.Results.EdgeColor);
            else
                patch('faces',M.faces,'vertices',M.vertices, ...
                    CW,f, ...
                    'FaceColor',FC, ...
                    'EdgeColor',p.Results.EdgeColor);
            end
            
            cameratoolbar; cameratoolbar('SetCoordSys','none');
            axis equal; axis off;
            colorbar; caxis(p.Results.Caxis);
            colormap(p.Results.Colormap);
            
            if p.Results.OpenGL == 1
                camlight('headlight');
                lighting phong;
                material dull;
                %                 light('Position',p.Results.LightPos,'Style','infinite');
            end
            
            if p.Results.Dock == 1; set(gcf, 'WindowStyle', 'docked'); end;
        end
        
        function vf(M,vf,varargin)
            if ~iscell(vf)
                vfs = {vf};
            else
                vfs = vf;
            end
            nvf = length(vfs);
            for i=1:nvf
                vfs{i} = reshape(vfs{i},[],3);
            end
            
            p = inputParser;
            addOptional(p,'FaceColor','w');
            addOptional(p,'EdgeColor','k');
            addOptional(p,'f',MESH.normv(vfs{1}));
            addOptional(p,'mesh',1);
            addOptional(p,'Color','k');
            addOptional(p,'Dock',0);
            addOptional(p,'nRosy',1);
            addOptional(p,'locs',{[1:M.nf]'});
            addOptional(p,'linewidth',1);
            addOptional(p,'threeD',0);
            addOptional(p,'Scale',0.5);

            parse(p,varargin{:});
            
            locs = p.Results.locs;
            for i=1:size(locs)
                assert(size(locs{i},1) == size(vfs{i},1),...
                    'locs should be the same length as vf');
            end
            
            func = p.Results.f;
            if length(locs) > 1 || size(locs{1},1) ~= M.nf
                func = [];
            end
            
            hold on;
            if p.Results.mesh == 1 
                if isempty(func) 
                    MESH_VIS.mesh(M, varargin{:});
                else
                    MESH_VIS.func(M, func);
                end
            end
            
                       
            X = M.vertices;
            if any(strcmp(properties(M), 'triangles')) == 1
                T = M.triangles;
                Xm = (X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;
            else
                Q = M.quads;
                Xm = (X(Q(:,1),:)+X(Q(:,2),:)+X(Q(:,3),:)+X(Q(:,4),:))/4;
            end
            
            
            scale = p.Results.Scale;
            allvfs = [];
            allxs = [];
            for j=1:nvf
%                 R = speye(3*M.nf,3*M.nf);
%                 if (p.Results.nRosy(j) > 1)
%                     for i = 1:(4/p.Results.nRosy(j)); R = M.R*R; end;
%                 end
                rlocs = [locs{j},M.nf+locs{j},2*M.nf+locs{j}];
% 
%                 % works only for n=2/4
%                 for i = 1:p.Results.nRosy(j)
%                     allxs = [allxs; [Xm(locs{j},1), Xm(locs{j},2), Xm(locs{j},3)]];
%                     allvfs = [allvfs; [vfs{j}(:,1), vfs{j}(:,2), vfs{j}(:,3)]];
% 
%                     vfs{j} = reshape(R(rlocs,rlocs)*vfs{j}(:),[],3);
%                 end
                n = p.Results.nRosy(j);
                a = linspace(0,2*pi,n+1)';
                vf2d = reshape(M.EB([locs{j},M.nf+locs{j}],rlocs)*vfs{j}(:),[],2);

                for i=1:n
                    vf2drc = exp(a(i)*1i)*(vf2d(:,1)+vf2d(:,2)*1i);
                    vf2dr = [real(vf2drc), imag(vf2drc)];
                    vfr = reshape(M.EBI(rlocs,[locs{j},M.nf+locs{j}])*vf2dr(:),[],3);
                    allxs = [allxs; ...
                        [Xm(locs{j},1), Xm(locs{j},2), Xm(locs{j},3)]];
                    allvfs = [allvfs; ...
                        [vfr(:,1), vfr(:,2), vfr(:,3)]];
                end
            end
            % requires quiver3D
            if p.Results.threeD
                locs = find(MESH.normv(allvfs) > 1e-5);
                locs = locs(1:20);
                quiver3D( allxs(locs,:), allvfs(locs,:), ...
                    p.Results.Color);                
            else
                quiver3( allxs(:,1), allxs(:,2), allxs(:,3), ...
                    allvfs(:,1), allvfs(:,2), allvfs(:,3), ...
                    scale, 'Color', p.Results.Color, 'linewidth', p.Results.linewidth);
            end            
            cameratoolbar; cameratoolbar('SetCoordSys','none');
            axis equal; axis off;
            colormap(jet);
            if p.Results.Dock == 1; set(gcf, 'WindowStyle', 'docked'); end;
            
            %             view([0 1 0]);
            
            hold off;
        end

        function mvf(M,F,varargin)
            n = size(F,2)/3;
            idx = reshape(1:n*3,3,n)';
            cols = {'r','w','k','m','b','y','c','g'};
            for i=1:n
                MESH_VIS.vf(M, {F(:, idx(i,:))}, varargin{:}, 'Color', cols{i}); 
                hold on
            end
        end
        

        function bdry(M, varargin)
            be = M.edges(M.ie == 0,:);
            MESH_VIS.edges(M, be, varargin{:});
        end
        
        
        function edges(M, edges, varargin)            
            hold on;
            X = [M.vertices(edges(:,1),1)'; M.vertices(edges(:,2),1)'];
            Y = [M.vertices(edges(:,1),2)'; M.vertices(edges(:,2),2)'];
            Z = [M.vertices(edges(:,1),3)'; M.vertices(edges(:,2),3)'];
            line(X,Y,Z,varargin{:});
            hold off
        end
        
        function text_line_annotation(i, varargin)
            annotation('textarrow','headstyle','none',...
                'textcolor','k',...
                'position', [0.05 1-i*0.05 .05 0],...
                'verticalalignment', 'middle',...
                'linewidth',3, varargin{:});            
        end

        function xyplane(bb,zc)
            if nargin < 2, zc = 0; end;
            s = (bb(2)-bb(1))/100;
            [x y] = meshgrid(bb(1):s:bb(2)); 
            z = ones(size(x, 1))*zc; 
            surf(x, y, z,'facealpha',.2,'FaceColor','b','edgecolor','none');            
        end
        
        
        function set_camera(ca,cam)
            set(ca, 'PlotBoxAspectRatio',cam.pba);
            set(ca, 'DataAspectRatio',cam.dar);
            set(ca, 'CameraViewAngle',cam.cva);
            set(ca, 'CameraUpVector',cam.cuv);
            set(ca, 'CameraTarget',cam.ct);
            set(ca, 'CameraPosition',cam.cp);
        end
        
        function cam = get_camera(ca)
            cam.pba = get(ca, 'PlotBoxAspectRatio');
            cam.dar = get(ca, 'DataAspectRatio');
            cam.cva = get(ca, 'CameraViewAngle');
            cam.cuv = get(ca, 'CameraUpVector');
            cam.ct = get(ca, 'CameraTarget');
            cam.cp = get(ca, 'CameraPosition');
        end
        
        function singular_pts( M, r, s )
            X = M.vertices;
            
            SI = find(s);
            
            [SX,SY,SZ] = sphere;
            a = max(max(X)-min(X))*r;
            SX = SX*a; SY = SY*a; SZ = SZ*a;
            
            [cm,cM] = caxis(gca);
            hold on;
            for i = 1:numel(SI)
                xs = X(SI(i),:);
                surf(SX+xs(1),SY+xs(2),SZ+xs(3),'FaceColor',...
                    MESH_VIS.choose_col(s(SI(i))),...
                    'EdgeColor','none','FaceAlpha',1);
            end
            caxis(gca,[cm,cM]);
            hold off;
        end
        
        function col = choose_col(ind)
            if ind == -3
                col = 'm';
            elseif ind == -2
                col = 'c';
            elseif ind == -1
                col = 'r';
            elseif ind == 1
                col = 'k';
            elseif ind == 2
                col = 'b';
            elseif ind == 3
                col = 'g';
            else
                col = 'y';
            end
        end
        
        function cm2 = interp_cm( cm, sz )
            n = numel(cm(:,1));
            xq = linspace(0,1,sz);
            cm2 = zeros(sz,3);
            for i = 1:3
                cm2(:,i) = interp1(cm(:,1),cm(:,i+1),xq);
            end
        end
    end
    
end


classdef MESH_IO
    
    properties
    end
    
    methods (Static)
        
        %%%%%%%%
        % Read %
        %%%%%%%%

        function [X,T] = robj(filename)
            OBJ = read_wobj( filename );
            X = OBJ.vertices;
            T = OBJ.objects(1).data.vertices;
        end
        
        function [X,T] = rply(filename)
            [T,X] = plyread( filename, 'tri' );
        end        
        
        % Can read off file of: only triangles or only quads (depending on
        % d)
        function [X,T] = roff(filename, d)
            if nargin < 2
                d = 3;
            end
            
            fid = fopen(filename,'r');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid);   % -1 if eof
            
            if strcmp(str(1:4), 'COFF')
                [X,T,~] = readCoff(filename,4); % assume 4 color channels
                return;
            end
            
            if ~strcmp(str(1:3), 'OFF')
                error('The file is not a valid OFF one.');
            end
           
            str = fgets(fid);
            sizes = sscanf(str, '%d %d', 2);
            while length(sizes) ~= 2
                str = fgets(fid);
                sizes = sscanf(str, '%d %d', 2);
            end
            nv = sizes(1);
            nf = sizes(2);
            
            % Read vertices
            [X,cnt] = fscanf(fid,'%lf %lf %lf\n', [3,nv]);
            if cnt~=3*nv
                warning('Problem in reading vertices.');
            end
            X = X';
            
            if d == 3
                [T,~] = fscanf(fid,'3 %ld %ld %ld\n', [3,inf]);
            elseif d == 4
                [T,~] = fscanf(fid,'4 %ld %ld %ld %ld\n', [4,inf]);
            else
                error('unsupported faces');
            end
                    
            T = T'+1;
            T = double(T);
            
            fclose(fid);
        end
        
        function [ V ] = rnrosy(filename)
            fid = fopen(filename,'r');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid);
            sizes = sscanf(str, '%d %d', 2);
            nv = sizes(1);
            
            [V,~] = fscanf(fid,'%lf %lf %lf', [3,nv]);
            
            fclose(fid);
        end
        
        function [ XF ] = rffield(filename)
            fid = fopen( [filename '.ffield'], 'r' );
            if( fid == -1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid); % Thresholds
            str = fgets(fid); % ??
            str = fgets(fid); % number of faces
            nf = sscanf(str,'%d');
            str = fgets(fid); % k1 k2 ...
            
            % read cross field
            [XF,cnt] = fscanf(fid,'%lf %lf %lf %lf %lf %lf %lf %lf\n', [8,nf]);
            XF = XF'; XF = XF(:,3:end);
            
            fclose( fid );
        end

        % Read libdirectional rawfield
        function F = rrawfield(filename)
            % First line is metadata
            d = importdata(filename ,'\t', 1);
            if ~isfield(d,'data')
                d = importdata(filename ,' ', 1);
            end
                
            F = d.data;
        end
           
        % Read singularities, format: %vid, %sing index,
        % vertex ids are zero based
        function [Sind, Sval] = rsing(filename)
            fid = fopen(filename,'r');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid);
            sizes = sscanf(str, '%d %d', 2);
            n = sizes(1); ns = sizes(2); 
            
            [X,~] = fscanf(fid,'%d %d\n', [2,ns]);

            Sind = X(1,:)'+1;
            Sval = X(2,:)';
            
            fclose(fid);
        end
        
        % Read libdirectional matching
        % EV, EF, FE - ligdirectional structures, see MESH.fromlibdirectional
        % n - n-symmetry
        % Mt, an array of length EV with the matching
        function [EV, EF, FE, n, Mt] = rmtch(filename)
            fid = fopen(filename,'r');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            str = fgets(fid);
            sizes = sscanf(str, '%d %d %d', 3);
            n = sizes(1); ne = sizes(2); nf = sizes(3);
            
            [X,~] = fscanf(fid,'%d %d %d %d %d\n', [5,ne]);
            EF = X(1:2,:)' + 1;
            EV = X(3:4,:)' + 1;
            Mt = X(5,:)';

            [FE,~] = fscanf(fid,'%d %d %d\n', [3,nf]);
            FE = FE' + 1;
            
            fclose(fid);            
        end
        
        
        %%%%%%%%%
        % Write %
        %%%%%%%%%

        % image files
        function wfigs(filename,mesh,F,varargin)
        
            p = inputParser;
            addOptional(p,'Titles',[]);
            addOptional(p,'Colormap','jet');
            addOptional(p,'CAxis','auto');
            addOptional(p,'CAxisV',[]);
            addOptional(p,'Plot','func');
            addOptional(p,'Montage',1);
            addOptional(p,'View',[]);
            addOptional(p,'Zoom',1);
            addOptional(p,'Resolution',1024);
            addOptional(p,'Camera',[]);
            addOptional(p,'OpenGL',0);
            addOptional(p,'LightPos',[0 0 1]);
            addOptional(p,'UVt',[]);
            addOptional(p,'FUVt',[]);
            addOptional(p,'SingularPts',0);
            addOptional(p,'SingularPtsR',20);
            addOptional(p,'AspectRatio',1);
            addOptional(p,'ARmode',[]);
            addOptional(p,'EdgeColor','k');
            parse(p,varargin{:});
            
            if strcmpi(p.Results.Plot,'mesh') == 1, sz = 1; 
            elseif strcmpi(p.Results.Plot,'func') == 1, sz = size(F,2); end;
            
            for i = 1:sz
                figure;
                if strcmp(p.Results.Plot,'mesh') == 1
                    
                    MESH_VIS.mesh( mesh, ...
                                   'CW',F,...
                                   'UVt',p.Results.UVt,...
                                   'FUVt',p.Results.FUVt,...
                                   'SingularPts',p.Results.SingularPts,...
                                   'SingularPtsR',p.Results.SingularPtsR,...
                                   'EdgeColor',p.Results.EdgeColor);
                    
                elseif strcmp(p.Results.Plot,'func') == 1; 
                    
                    MESH_VIS.func( mesh, F(:,i),...
                                   'Colormap', p.Results.Colormap, 'EdgeColor', p.Results.EdgeColor );
                    colorbar off;
                    if strcmp(p.Results.CAxis,'scaled') == 1
                        ca = [min(min(F)),max(max(F))];
                        caxis( ca );
                    elseif strcmp(p.Results.CAxis,'manual') == 1
                        caxis(p.Results.CAxisV);
                    end
                    set(gca,'xlim',[min(mesh.vertices(:,1)),max(mesh.vertices(:,1))]);
                    set(gca,'ylim',[min(mesh.vertices(:,2)),max(mesh.vertices(:,2))]);
                    set(gca,'zlim',[min(mesh.vertices(:,3)),max(mesh.vertices(:,3))]);
                    set(gca, 'LooseInset', get(gca, 'TightInset'));
                end
                
                if isempty(p.Results.Camera) == 0
                    MESH_VIS.set_camera(gca,p.Results.Camera);
                end
                
                if p.Results.OpenGL == 1
%                     lighting gouraud; % lighting phong;
                    material dull;
%                     light('Position',p.Results.LightPos,...
%                           'Style','local','Color',[.8,.8,.8]);
                    camlight('headlight');
                end
                
                if numel(p.Results.Titles) > 0
                    title(p.Results.Titles{i},'FontSize',25,...
                          'interpreter','latex',...
                          'Units', 'normalized', 'Position', [0.5, 1, 0]); 
                      ax = gca;
                      ax.CameraViewAngle = 10;
                end
                
%                 zoom(p.Results.Zoom);
                
                if sz > 1
                    pngname = sprintf('%s_%04d.png',filename,i);
                else
                    pngname = sprintf('%s.png',filename);
                end
                
                if strcmp(p.Results.ARmode, 'M')
                    AspectRatio = (min(mesh.vertices(:,1))-max(mesh.vertices(:,1)))/(min(mesh.vertices(:,2))-max(mesh.vertices(:,2)));
                else
                    AspectRatio = p.Results.AspectRatio;
                end
                
                dpi = get(0, 'ScreenPixelsPerInch');
                in = p.Results.Resolution/dpi;

                fig = gcf;
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 AspectRatio*in in];
                set(gca,'DataAspectRatio',[1 1 1])
                fig.PaperPositionMode = 'manual';
                print(pngname,'-dpng','-r0')
                close all;
            end
            
            % montage
            if p.Results.Montage == 1
                rows = 1;
                cols = size(F,2) / rows;

                A = [];
                for i = 1:rows
                    b = [];
                    for j = 1:cols
                        pngname = sprintf('%s_%04d.png',filename,cols*(i-1)+j);
                        a = imread(pngname);
                        b = cat(2,b,a);
                    end
                    A = cat(1,A,b);
                end
                imwrite(A,[filename '.png']);
                delete( sprintf('%s_*.png',filename) );
                close all;
            end
%             set(gca,'FontSize',15)
        end
        
        % geometry files
        function wobj(filename,X,T,UV,T_uv,imfname)
            usetexturecoords = 0;
            usematerial = 0;
            objid = 1;
            
            if nargin > 3
                usetexturecoords = 1;
            end
            if nargin > 5
                usematerial = 1;
            end
            if usematerial
                 % Set material
                material(1).type='newmtl';
                material(1).data='default';
                material(2).type='Ka';
                material(2).data=[0.8 0.4 0.4];
                material(3).type='Kd';
                material(3).data=[0.8 0.4 0.4];
                material(4).type='Ks';
                material(4).data=[1 1 1];
                material(5).type='illum';
                material(5).data=2;
                material(6).type='Ns';
                material(6).data=27;
                material(6).type='map_Kd';
                material(6).data=imfname;
                OBJ.material = material;
                OBJ.objects(objid).type='usemtl';
                OBJ.objects(objid).data='default';
                objid = objid + 1;
            end
                        
            OBJ.vertices = X;
            OBJ.objects(objid).data.vertices = T;
            OBJ.objects(objid).type = 'f';
            if usetexturecoords 
                OBJ.vertices_texture = UV;
                OBJ.objects(objid).data.texture = T_uv;            
            end
            
            write_wobj( OBJ, filename );
        end
        
        function woff(filename,X,T,varargin)
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end

            nv = size(X,1);
            nf = size(T,1);

            fprintf(fid, 'OFF\r\n');
            fprintf(fid, '%d %d 0\r\n', nv, nf);

            fprintf(fid, '%.9f %.9f %.9f\r\n', X');
            fprintf(fid, '3 %d %d %d\r\n', (T-1)');

            fclose(fid);
        end

        function woffq(filename,X,Q,varargin)
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end

            nv = size(X,1);
            nf = size(Q,1);

            fprintf(fid, 'OFF\r\n');
            fprintf(fid, '%d %d 0\r\n', nv, nf);

            fprintf(fid, '%.9f %.9f %.9f\r\n', X');
            fprintf(fid, '4 %d %d %d %d\r\n', (Q-1)');

            fclose(fid);
        end        
        
        function woffp2(filename, X, F, varargin)
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            nv = size(X,1);
            nf = size(F,1);

            ff = {}; 
            for i=1:nf
                ff{end+1} = [];
                for j=1:size(F,2)
                    if ~isnan(F(i,j))
                        ff{end}(end+1) = F(i,j);            
                    end
                end
            end

            fprintf(fid, 'OFF\r\n');
            fprintf(fid, '%d %d 0\r\n', nv, nf);

            fprintf(fid, '%.9f %.9f %.9f\r\n', X');
            
            for i=1:nf
                vv = [length(ff{i}); ff{i}'-1];
                fprintf(fid, '%d ', vv);
                fprintf(fid, '\r\n');
            end
            fclose(fid);
        end
        
        
        
        function woffp(filename, X, F, varargin)
            ff = {}; nf = size(F,1);
            for i=1:nf
                ff{end+1} = [];
                for j=1:size(F,2)
                    if ~isnan(F(i,j))
                        ff{end}(end+1) = F(i,j);            
                    end
                end
            end

            newv = [];
            nv = size(X,1);
            tri = zeros(0,3);
            for i=1:nf
                lf = length(ff{i});
                nver = mean(X(ff{i},:));
                newv(end+1,:) = nver;
                for j=1:lf
                    jp1 = j+1; if jp1 > lf; jp1 = 1; end;
                    tt = [ff{i}(j), ff{i}(jp1), nv+i];
                    tri(end+1,:) = tt;
                end
            end
            MESH_IO.woff(filename,[X;newv],tri);            
        end
        
        function wnrosy(filename,M,V,varargin)
            p = inputParser;
            addOptional(p,'n',1);
            parse(p,varargin{:});
            
%             fid = fopen(filename,'w');
%             if( fid==-1 )
%                 error('Cannot open the file.');
%                 return;
%             end
%             
%             nv = size(V,1);
%             fprintf(fid,'%d %d\n', nv, p.Results.n);
%             fprintf(fid,'%.9f %.9f %.9f\r\n', V');
%             
%             fclose(fid);
            % Based on Danielle's code
            export_to_lic(filename,M,V);
        end
        
        function wvf(filename, V)
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            fprintf(fid,'%.9f %.9f %.9f\r\n', V');
            
            fclose(fid);
            
        end

        function wframe(filename, V1, V2)
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end
            
            fprintf(fid,'%.9f %.9f %.9f %.9f %.9f %.9f\r\n', [V1,V2]');
            
            fclose(fid);
        end

        % For libdirectional. F can be any set of vectors per face.
        function wrawfield(filename, F)
            fid = fopen(filename,'wt');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end

            fprintf(fid, '%d %d\n', size(F, 2) / 3, size(F,1));
            fprintf(fid, [repmat('%f\t', 1, size(F, 2)) '\n'], F');
            fclose(fid);
        end
        
        % For libdirectional, 
        % EV, EF, FE - ligdirectional structures, see MESH.tolibdirectional
        % n - n-symmetry
        % Mt, an array of length EV with the matching
        function wmtch(filename, EV, EF, FE, n, Mt )
            ne = size(EV,1);
            nf = size(FE,1);
            
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end

            fprintf(fid, '%d %d %d\n', n, ne, nf);
            fprintf(fid, '%d %d %d %d %d\n', [EF-1, EV-1, Mt]');
            fprintf(fid, '%d %d %d\n', FE'-1);

            fclose(fid);        
        end
        
        % write singularities, format: %vid, %sing index,
        % vertex ids are zero based
        function wsing(filename, Sind, Sval, n);
            ns = size(Sind,1);
            fid = fopen(filename,'w');
            if( fid==-1 )
                error('Cannot open the file.');
                return;
            end

            fprintf(fid, '%d %d\n', n, ns);
            fprintf(fid, '%d %d\n', [Sind-1, Sval]');

            fclose(fid);        
        end
        
        function wvtk(filename,X,T,varargin)
            p = inputParser;

            addOptional(p,'Xn',[]);
            addOptional(p,'func',[]);
            addOptional(p,'func_str','vc');
            addOptional(p,'timestep',0);
            addOptional(p,'CM',[]);

            parse(p,varargin{:});
            
            nv = length(X);
            nf = length(T);

            fid = fopen(filename,'w');

            % standard header
            fprintf(fid,'# vtk DataFile Version 3.0\n');
            fprintf(fid,'vtk output\n');
            fprintf(fid,'ASCII\n');
            fprintf(fid,'DATASET POLYDATA\n');

            % if available, add time annotation
            if isempty( p.Results.timestep ) == 0
                fprintf(fid,'FIELD FieldData 1\n');
                fprintf(fid,'TIME 1 1 double\n');
                fprintf(fid,'%g\n',p.Results.timestep);
            end

            % geometry data
            fprintf(fid,'POINTS %d float\n', nv);
            fprintf(fid,'%g %g %g\n', X');
            fprintf(fid,'POLYGONS %d %d\n', nf, 4*nf);
            fprintf(fid,'3 %d %d %d\n', T'-1);
            fprintf(fid,'\n');

            % add point data
            if isempty( p.Results.func ) == 0 || isempty( p.Results.Xn ) == 0
                fprintf(fid,'POINT_DATA %d\n', nv);

                if isempty( p.Results.func ) == 0
                    fprintf(fid,'SCALARS %s float 1\n', p.Results.func_str);
                    lt = 'default'; 
                    if isempty( p.Results.CM ) == 0; lt = 'lt'; end;
                    fprintf(fid,'LOOKUP_TABLE %s\n', lt);
                    fprintf(fid,'%g\n', p.Results.func);
                end

                if isempty( p.Results.Xn ) == 0
                    fprintf(fid,'NORMALS point_normals float\n');
                    fprintf(fid,'%g %g %g\n', p.Results.Xn');
                end

                fprintf(fid,'\n');
            end

            if isempty( p.Results.CM ) == 0
                fprintf(fid,'LOOKUP_TABLE lt %d\n', length(p.Results.CM));
                fprintf(fid,'%g %g %g 1\n', p.Results.CM');
                fprintf(fid,'\n');
            end

            fclose(fid);
        end
        
        % video files
        function wavi(dirname,filename,varargin)
            
            p = inputParser;
            addOptional(p,'FrameRate',24);
            addOptional(p,'Quality',100);
            addOptional(p,'ext','png');
            parse(p,varargin{:});
            
            outputVideo = VideoWriter( [dirname filename(1:end-1) '.avi'] );
            
            outputVideo.FrameRate = p.Results.FrameRate;
            outputVideo.Quality = p.Results.Quality;

            open(outputVideo);
            
            D = dir([dirname filename '*.' p.Results.ext]);
            N = length(D);
            for i = 1:N
                img = imread( [dirname D(i).name], ...
                              'BackgroundColor', [1,1,1,] );
                writeVideo(outputVideo,img);
            end
            
            close(outputVideo);
        end
        
        function wmpeg4(dirname,filename,varargin)
            p = inputParser;
            addOptional(p,'FrameRate',24);
            addOptional(p,'Quality',100);
            addOptional(p,'ext','png');
            parse(p,varargin{:});
            
            outputVideo = VideoWriter( [dirname filename(1:end-1) '.mp4'], ...
                                       'MPEG-4' );
            
            outputVideo.FrameRate = p.Results.FrameRate;
            outputVideo.Quality = p.Results.Quality;

            open(outputVideo);
            
            D = dir([dirname filename '*.' p.Results.ext]);
            N = length(D);
            for i = 1:N
                img = imread( [dirname D(i).name],'BackgroundColor','none');
                
                % manual padding to a multiple of 2
                sz = size(img);
                if mod(sz(1),2) == 1; img = [img; 255*ones(1,sz(2),3)]; end;
                
                sz = size(img);
                if mod(sz(2),2) == 1; img = [img 255*ones(sz(1),1,3)]; end;
                
                writeVideo(outputVideo,img);
            end
            
            close(outputVideo);
        end
       
        function montage( filename, imgnames, rows, cols, gap, suff )
            if nargin < 6
                suff = '.png';
            end                
            if nargin < 5
                rows = 1; cols = 4; gap = [0,0]; 
            end

            A = []; tt = 1;
            for i = 1:rows
                b = [];
                for j = 1:cols
                    a = imread( [imgnames{tt} suff] ); tt = tt+1;
                    if j > 1
                        x = size(b,1)-size(a,1);
                        if x < 0
                            % b less columns than a
                            x = abs(x);
                            xt = floor(x/2); xb = x - xt;
                            gapt = 255 * ones(xt, size(b,2),  3, 'uint8');
                            gapb = 255 * ones(xb, size(b,2),  3, 'uint8');
                            b = cat(1,gapt,b,gapb);
                        elseif x > 0
                            % a less columns than b
                            xt = floor(x/2); xb = x - xt;
                            gapt = 255 * ones(xt, size(a,2),  3, 'uint8');
                            gapb = 255 * ones(xb, size(a,2),  3, 'uint8');
                            a = cat(1,gapt,a,gapb);
                        end
                    end
                    b = cat(2,b,a);
                    if j < cols && gap(2) ~= 0
                        gapc = 255 * ones(size(b,1), gap(2), 3, 'uint8');
                        b = cat(2,b,gapc);
                    end
                end
                if i > 1
                    x = size(A,2)-size(b,2);
                    if x < 0
                        % A less columns than b
                        x = abs(x);
                        xl = floor(x/2); xr = x - xl;
                        gapl = 255 * ones(size(A,1), xl, 3, 'uint8');
                        gapr = 255 * ones(size(A,1), xr, 3, 'uint8');
                        A = cat(2,gapl,A,gapr);
                    elseif x > 0
                        % b less columns than A
                        xl = floor(x/2); xr = x - xl;
                        gapl = 255 * ones(size(b,1), xl, 3, 'uint8');
                        gapr = 255 * ones(size(b,1), xr, 3, 'uint8');
                        b = cat(2,gapl,b,gapr);
                    end
                end                    
                A = cat(1,A,b);
                if i < rows && gap(1) ~= 0
                    gapr = 255 * ones(gap(1), size(A,2), 3, 'uint8');
                    A = cat(1,A,gapr);
                end
            end
            imwrite(A,[filename suff]);
        end
        
        function [Nv] = rNvcsv(filename, nv)
            
            Nv = csvread(filename);
            
            if size(Nv,1)~=nv || size(Nv,2)~=3
                warning('Problem in reading vertices.');
            end
        end
    end
    
end


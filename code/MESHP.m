% Polygonal mesh
classdef MESHP < handle
    
    properties (Access='public')
        
        name
        
        vertices
        faces
        
        nv              % #vertices
        nf              % #faces
        
    end
    
    methods
        
        function [ mesh ] = MESHP( name )
                        
            mesh.name = name;

            [mesh.vertices, mesh.faces] = readMesh_off([name '.off']);                                    
            if iscell(mesh.faces)
                C = mesh.faces;
                v = cellfun(@length, C);
                mv = max(v);
                C = cellfun(@(x)[x(1:end), NaN(1,mv-length(x))]',C,'UniformOutput',false);
                mesh.faces = cell2mat(C)';
            end
            
            
            mesh.nv = size(mesh.vertices,1);
            mesh.nf = size(mesh.faces,1);
            
        end
        
        
        function [ fa ] = face_areas( mesh )
            X = mesh.vertices;
            F = mesh.faces;
            
            fa = zeros(mesh.nf,1);
            for i=1:mesh.nf
                f = F(i,:);
                f = f(~isnan(f));
                fa(i) = polygonArea3d(X(f,:));
            end            
        end

        function [planarityFaces] = planarity_general( mesh ) % D
            V = mesh.vertices;
            F = mesh.faces;

            planarityFaces = zeros(size(F, 1), 1);
            for fid = 1:size(F, 1)
                curF = F(fid,:); curF(isnan(curF)) = []; dim = length(curF);
                quadPlan = zeros(dim, 1);
                for d = 1:dim
                    v1 = V(F(fid, d),:);
                    v2 = V(F(fid, mod(d, dim) + 1),:);
                    v3 = V(F(fid, mod(d + 1, dim) + 1),:);
                    v4 = V(F(fid, mod(d + 2, dim) + 1),:);

                    diagCross = cross(v3 - v1, v4 - v2);
                    denom = norm(diagCross) * (norm(v3 - v1) + norm(v4 - v2)) / 2.;

                    if(abs(denom) > 10E-10)
                        quadPlan(d) = dot(diagCross, v2 - v1) / denom;
                    else
                        quadPlan(d) = 0;
                    end
                end
                planarityFaces(fid) = 100 * sqrt(norm(quadPlan)^2 / dim);
            end
        end


    end
    
    methods (Static)
    end
end
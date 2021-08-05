classdef MESHQ < handle
    
    properties (Access='public')
        
        name
        
        vertices
        quads
        
        nv              % #vertices
        nf              % #faces
        
        AM,valence
    end
    
    methods
        
        function [ mesh ] = MESHQ( name )
                        
            mesh.name = name;            
                        
            mesh.nv = size(mesh.vertices,1);
            mesh.nf = size(mesh.quads,1);
            
            mesh.AM = adj_mat( mesh );
            mesh.valence = sum(mesh.AM,2);
        end
        
        function [ AM ] = adj_mat( mesh )
            T = double( mesh.quads );
            
            I = [T(:,2);T(:,3);T(:,4);T(:,1)];
            J = [T(:,3);T(:,4);T(:,1);T(:,2)];
            S = ones(size(I(:),1),1);
            AM = sparse(I,J,S,mesh.nv,mesh.nv);
        end
        
    end
    
    methods (Static)
    end
end
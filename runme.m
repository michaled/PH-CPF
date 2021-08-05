
clearvars; close all; clc;

addpath(genpath(pwd));  % code, data, external

n = 6;                  % 4 - for quads, 6 - hexes

% Set params:
delta = 0.005;          % scaling parameter: max distance from surface, as a fraction of model's bounding box
eta = 0.002;            % scaling parameter: area of max element, as a fraction of total surface area
mp = 1;                 % maximal planarity, in precentage

datadir = 'data/';
outdir = sprintf('results_%s/', datestr(now,'mm-dd-yyyy'));
if ~exist(outdir, 'dir')
    mkdir(outdir);
    mkdir([outdir,'tmp/']);
end

params = paramsConst(n, delta, eta, mp, outdir, datadir);

    
D = dir([datadir '*.off']);
Dnames = {D.name}; [~, I] = sort(lower(Dnames)); DnamesS = Dnames(I)';
dvec = 1:length(Dnames);

pl_cm = load('planarity_colormap.mat'); pl_cm = pl_cm.cm1;

for i = dvec
    name = DnamesS{i};
    name = name(1:end-4);
    disp(name);
        
    cpf = CPF(name);

    cpf.set_param(params);
        
    out_data = [outdir name '_' cpf.get_param '.mat'];

    success = 0; ME = [];
    while ~success
        try  
            cpf.solve();
            success = 1;
        catch ME
            if contains(ME.message,'Barrier') & cpf.iter_opt<3
                warning(['Barrier failed, perr = %.2g, '...
                    'rerunning with iter_opt = %d'], ...
                    cpf.perr, cpf.iter_opt+1);
                cpf.iter_opt = cpf.iter_opt+1;
            else
                success = 2; % didn't really succeed
                fprintf('Mesh %s failed.',name);
                save(out_data,'cpf','ME');
            end
        end
    end
    if success == 2
        continue;
    end
    save(out_data,'cpf');
    
    % Vis results:
    ff = cpf.Mp.planarity_general;
    figure;
    MESH_VIS.func(cpf.Mp, ff, 'edgecolor', 'k'); 
    colormap(pl_cm); caxis([0,mp]);
end




function params = paramsConst(n, delta, eta, mp, outdir, datadir)

params = {'n', n};                          % 6 - for hexes, 4 - for quads
params = [params, {'delta'}, {delta}];      % scaling parameter: max distance from surface, as a fraction of model's bounding box
params = [params, {'eta'}, {eta}];          % scaling parameter: area of max element, as a fraction of total surface area    
params = [params, {'mp'}, {mp}];            % maximal planarity, in precentage

% Energy weights:
params = [params, {'lsm'}, {1}];            % weight for smoothness energy
params = [params, {'lc'}, {1}];             % weight for continuity energy 
params = [params, {'lsi'}, {1}];            % weight for sizing energy 
params = [params, {'la'}, {1}];             % weight for alignment energy 
params = [params, {'lin'}, {1}];            % weight for injectivity energy 
params = [params, {'lcj'}, {1}];            % weight for conjugacy energy 
params = [params, {'lo'}, {1}];             % weight for orthogonality energy 

params = [params, {'siz'}, {'k'}];          % siz - which sizing to use for the vector fields, 'i'-identity, 'k'-curvature, 'h'-hodge

params = [params, {'optinit'}, {'k'}];      % U,V initialization method: 'r' - random, 'k' - guiding field     

params = [params, {'datadir'}, {datadir}];  % data folder
params = [params, {'outdir'}, {outdir}];    % ouutput folder

% Additional optional inputs:
% lr                                          % edge length for parameterization and meshing,
                                            % if not set, it is computed automaticly
% rhow                                        % alignment's weighting
                                            % 0 - no weighing, 1 - weighting with RHO squared, 2 - polynomial weighting
% dualize                                     % whether to dualize the remeshed mesh, relevant only for quads, default 1
% planar                                      % whether to planarize the remeshed mesh, default 1
% figs                                        % whether to open figures while running, default 0



end
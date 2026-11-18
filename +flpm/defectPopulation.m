function pop = defectPopulation(type, varargin)
%FLPM.DEFECTPOPULATION  Statistical descriptor of a defect population.
%
%   The model is driven by population descriptors, not defect lists. Defect
%   sets are drawn from them by flpm.samplePopulation / flpm.monteCarlo.
%   Locations are measured from the free surface; no specimen dimension is
%   needed.
%
%   VOLUMETRIC (pores, lack of fusion; XCT)
%     pop = flpm.defectPopulation('volumetric', ...
%             'Size',  makedist('Lognormal','mu',log(0.06),'sigma',0.2), ... % sqrt(area), mm
%             'Depth', makedist('Gamma','a',2,'b',0.12), ...    % centre depth, mm
%             'CountPerLength', 0.5);                            % per mm of surface
%
%   SURFACE (roughness valleys; profilometry)
%     pop = flpm.defectPopulation('surface', ...
%             'Size', {'GeneralizedPareto','k',-0.29,'sigma',0.0147,'theta',0}, ...
%             'SizeOffset', 0.050, ...                            % EVT threshold u
%             'CountPerLength', 0.37, ...                         % V, per mm
%             'ShieldingFactor', 1.02);                           % Y_R
%
%   A descriptor is a distribution object (makedist/fitdist, including
%   'Kernel'), a data vector (kernel fit), a scalar (constant) or a makedist
%   argument cell.
%
%   Name-value pairs
%     'Size'            volumetric: sqrt(area) of the defect, mm
%                       surface: valley depth Rv, mm                 (required)
%     'SizeOffset'      added to every sampled size (EVT threshold u)     [0]
%     'Depth'           volumetric: depth of the defect centre below the
%                       free surface, mm                             (required)
%     'Count'           defects per simulated section, or
%     'CountPerLength'  defects per mm of free-surface length; the count is
%                       CountPerLength * SurfaceLength, or
%     'CountPerArea'    defects per mm^2 of free surface (3-D sampling); the
%                       count is CountPerArea * SurfaceLength * AxialLength
%     'CountModel'      'fixed' | 'poisson'                          ['fixed']
%     'Threshold'       discard sampled sizes <= u (peaks over threshold)
%     'SizeBounds', 'DepthBounds'  truncation [1e-6 Inf], [0 Inf]
%     'ShieldingFactor' surface: Y_R of the roughness valleys, defined with
%                       the valley depth, K = Y_R s sqrt(pi Rv)
%                       (flpm.roughnessShielding); 1.12 = unshielded notch [1.12]
%     'Clustering'      flpm.clustering(...) or a fit from
%                       flpm.fitClustering: Thomas cluster process for the
%                       locations; the count stays the mean total number of
%                       defects. [] = independent (Poisson) locations     [[]]
%     'Copula'          volumetric: dependence of size and depth,
%                       {'Gaussian', rho}, {'t', rho, nu} or the output of
%                       flpm.testDependence; [] = independent            [[]]
%     'Name'            label                                           [type]
%
%   See also flpm.samplePopulation, flpm.monteCarlo, flpm.roughnessShielding,
%   flpm.clustering, flpm.fitClustering.

type = validatestring(type, {'volumetric', 'surface'});
ip = inputParser;
ip.FunctionName = 'flpm.defectPopulation';
ip.addParameter('Size', []);
ip.addParameter('SizeOffset', 0, @(x) isscalar(x) && x >= 0);
ip.addParameter('Depth', []);
ip.addParameter('Count', []);
ip.addParameter('CountPerLength', []);
ip.addParameter('CountPerArea', []);
ip.addParameter('CountModel', 'fixed');
ip.addParameter('Threshold', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
ip.addParameter('SizeBounds', [1e-6 Inf]);
ip.addParameter('DepthBounds', [0 Inf]);
ip.addParameter('ShieldingFactor', 1.12, @(x) isscalar(x) && x > 0);
ip.addParameter('Clustering', []);
ip.addParameter('Copula', []);
ip.addParameter('Name', type);
ip.parse(varargin{:});
r = ip.Results;

if isempty(r.Size)
    error('flpm:population:size', 'A ''Size'' descriptor is required.');
end
if ~isempty(r.Count) + ~isempty(r.CountPerLength) + ~isempty(r.CountPerArea) ~= 1
    error('flpm:population:count', ...
        'Give exactly one of ''Count'', ''CountPerLength'' or ''CountPerArea''.');
end
cnt = [r.Count r.CountPerLength r.CountPerArea];
if ~(isscalar(cnt) && cnt >= 0)
    error('flpm:population:count', 'The defect count must be a non-negative scalar.');
end
pop.type = type;
pop.name = char(r.Name);
pop.size = flpm.internal.toDistribution(r.Size, 'Size');
pop.sizeOffset = r.SizeOffset;
pop.sizeBounds = r.SizeBounds;
pop.count = r.Count;
pop.countPerLength = r.CountPerLength;
pop.countPerArea = r.CountPerArea;
pop.countModel = validatestring(r.CountModel, {'fixed', 'poisson'});
pop.threshold = r.Threshold;
pop.clustering = [];
if ~isempty(r.Clustering)
    c = r.Clustering;
    if isstruct(c) && isfield(c, 'clustering'), c = c.clustering; end
    if ~(isstruct(c) && isfield(c, 'type') && strcmp(c.type, 'thomas'))
        error('flpm:population:clustering', ...
            '''Clustering'' must come from flpm.clustering or flpm.fitClustering.');
    end
    pop.clustering = c;
end
pop.copula = flpm.internal.toCopula(r.Copula);
if ~isempty(pop.copula) && (~strcmp(type, 'volumetric') || ~isempty(pop.clustering))
    error('flpm:population:copula', ...
        'A size-depth copula needs a volumetric population without clustering.');
end
switch type
    case 'volumetric'
        if isempty(r.Depth)
            error('flpm:population:depth', 'Give a ''Depth'' descriptor for a volumetric population.');
        end
        pop.depth = flpm.internal.toDistribution(r.Depth, 'Depth');
        pop.depthBounds = r.DepthBounds;
        pop.shielding = [];
    case 'surface'
        pop.depth = [];
        pop.depthBounds = [];
        pop.shielding = r.ShieldingFactor;
end
end

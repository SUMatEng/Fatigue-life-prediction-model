function c = clustering(varargin)
%FLPM.CLUSTERING  Thomas (Neyman-Scott) cluster process for defect locations.
%
%   c = flpm.clustering('MeanClusterSize', mu, 'Spread', sigma)
%   c = flpm.clustering('MeanClusterSize', mu, 'Spread', [sigmaPlane sigmaAxial], ...
%                       'ParentDepth', makedist('Exponential', 'mu', 0.4))
%   pop = flpm.defectPopulation('volumetric', ..., 'Clustering', c)
%
%   Defect locations follow a Thomas process: cluster centres (parents)
%   form a Poisson process of intensity kappa; each parent has a Poisson
%   number of defects with mean mu, displaced from it by independent
%   Gaussian offsets of standard deviation sigma. Sizes and all other marks
%   are drawn per defect, independently of the location.
%
%   The population's count (Count, CountPerLength or CountPerArea) is the
%   mean total number of defects, lambda. The parent intensity is
%   kappa = lambda/mu, so clustering changes where defects are, not how
%   many there are. Poisson (independent) locations are the limit mu -> 1
%   and the default when no 'Clustering' is given.
%
%   Name-value pairs
%     'MeanClusterSize'  mu, mean number of defects per cluster, >= 1 (required)
%     'Spread'           sigma, mm (required): a scalar (isotropic), or
%                        [sigmaPlane sigmaAxial] with sigmaPlane in the plane
%                        of the cross-section (depth and position along the
%                        surface) and sigmaAxial along the load axis. For a
%                        bar built vertically, the load axis is the build
%                        direction.
%     'ParentDepth'      descriptor of the depth of the cluster centres (mm);
%                        [] uses the population's 'Depth'              [[]]
%
%   A fit from flpm.fitClustering can be passed instead: its field
%   'clustering' holds this struct.
%
%   See also flpm.fitClustering, flpm.defectPopulation, flpm.samplePopulation.

ip = inputParser;
ip.FunctionName = 'flpm.clustering';
ip.addParameter('MeanClusterSize', [], @(x) isscalar(x) && x >= 1);
ip.addParameter('Spread', [], @(x) isnumeric(x) && any(numel(x) == [1 2]) && all(x > 0));
ip.addParameter('ParentDepth', []);
ip.parse(varargin{:});
r = ip.Results;
if isempty(r.MeanClusterSize) || isempty(r.Spread)
    error('flpm:clustering:input', '''MeanClusterSize'' and ''Spread'' are required.');
end
c.type = 'thomas';
c.meanClusterSize = r.MeanClusterSize;
c.spread = [r.Spread(1), r.Spread(end)];        % [plane axial], mm
c.parentDepth = [];
if ~isempty(r.ParentDepth)
    c.parentDepth = flpm.internal.toDistribution(r.ParentDepth, 'ParentDepth');
end
end

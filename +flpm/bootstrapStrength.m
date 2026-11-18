function B = bootstrapStrength(data, build, mat, varargin)
%FLPM.BOOTSTRAPSTRENGTH  Confidence bands of the strength from the defect data.
%
%   B = flpm.bootstrapStrength(data, build, mat, 'SurfaceLength', L)
%   B = flpm.bootstrapStrength(..., 'NumBootstrap', 200, 'Seed', 1, ...
%                              'AxialLength', Lz, 'Quantiles', [0.1 0.5 0.9])
%
%   data   table (or matrix) of measured defects, one row per defect
%   build  function handle: pops = build(dataRows) fits the population
%          descriptors to a set of rows (e.g. kernel or truncated fits)
%
%   The defect rows are resampled with replacement NumBootstrap times. For
%   each replicate the descriptors are refitted with build and the
%   strength quantiles are evaluated in closed form (flpm.weakestLink,
%   with one arrest map for all replicates). The spread of the quantiles
%   over the replicates is the uncertainty that comes from the finite
%   defect sample alone.
%
%   Output struct
%     quantiles    p levels
%     estimate     quantiles from the full data
%     replicates   NumBootstrap x numel(p)
%     ci95         2 x numel(p), percentile intervals
%     map          the arrest map used

ip = inputParser;
ip.addParameter('SurfaceLength', []);
ip.addParameter('AxialLength', 0);
ip.addParameter('NumBootstrap', 200);
ip.addParameter('Seed', 0);
ip.addParameter('Quantiles', [0.1 0.5 0.9]);
ip.addParameter('Map', []);
ip.addParameter('Options', []);
ip.parse(varargin{:});
o = ip.Results;
q = o.Quantiles(:).';
args = {'SurfaceLength', o.SurfaceLength, 'AxialLength', o.AxialLength, 'Options', o.Options};

W0 = flpm.weakestLink(build(data), mat, args{:}, 'Map', o.Map);
M = W0.map;
if ~iscell(M), Mc = {M}; else, Mc = M; end
lims = cellfun(@(m) m.size([1 end]).', Mc, 'UniformOutput', false);
lims = [max(cellfun(@(l) l(1), lims)), min(cellfun(@(l) l(2), lims))];
B.quantiles = q;
B.estimate = W0.quantile(q);
rs = RandStream('mrg32k3a', 'Seed', o.Seed);
n = size(data, 1);
R = nan(o.NumBootstrap, numel(q));
for b = 1:o.NumBootstrap
    rows = data(randi(rs, n, n, 1), :);
    pops = build(rows);
    try
        W = flpm.weakestLink(pops, mat, args{:}, 'Map', M);
        R(b, :) = W.quantile(q);
    catch e
        if ~strcmp(e.identifier, 'flpm:arrestMap:range'), rethrow(e), end
        % refit produced sizes outside the map: extend it once
        M = rebuild(Mc, mat, o.Options, lims);
        W = flpm.weakestLink(pops, mat, args{:}, 'Map', M);
        R(b, :) = W.quantile(q);
    end
end
B.replicates = R;
B.ci95 = prctile(R, [2.5 97.5], 1);
B.map = M;
end

function M = rebuild(Mc, mat, opts, lims)
[mats, ~] = flpm.internal.thresholdNodes(mat);
M = cellfun(@(m) flpm.arrestMap(m, opts, 'SizeRange', [0.5 2] .* lims), mats, 'UniformOutput', false);
if numel(M) == 1, M = M{1}; end
end

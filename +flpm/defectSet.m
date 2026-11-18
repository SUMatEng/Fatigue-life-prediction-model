function set = defectSet(sizeMm, radiusMm, angle, aspect, varargin)
%FLPM.DEFECTSET  Build a defect set from explicit vectors.
%
%   set = flpm.defectSet(size, radius, angle, aspect)
%   set = flpm.defectSet(size, radius, angle, aspect, 'YFloor', yR)
%
%   Convenience for deterministic studies and for re-running the sampled
%   populations of the 2023-2025 scripts (rows of sim_XCT_*.xlsx). NaN
%   entries (padding of variable-length rows) are removed.
%     size    sqrt(area), mm         radius  radial position of centre, mm
%     angle   rad                    aspect  a/c
ip = inputParser;
ip.addParameter('YFloor', []);
ip.addParameter('Population', []);
ip.parse(varargin{:});
v = {sizeMm(:), radiusMm(:), angle(:), aspect(:)};
n = numel(v{1});
if any(cellfun(@numel, v) ~= n)
    error('flpm:defectSet:length', 'All vectors must have the same length.');
end
keep = ~any(isnan([v{:}]), 2);
set = struct('size', v{1}(keep), 'radius', v{2}(keep), 'angle', v{3}(keep), ...
    'aspect', v{4}(keep));
if ~isempty(ip.Results.YFloor)
    y = ip.Results.YFloor(:) .* ones(n, 1);
    set.yFloor = y(keep);
end
if ~isempty(ip.Results.Population)
    p = ip.Results.Population(:) .* ones(n, 1);
    set.population = p(keep);
end
end

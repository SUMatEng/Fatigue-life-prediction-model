function set = defectSet(sizeMm, depthMm, positionMm, varargin)
%FLPM.DEFECTSET  Build a defect set from explicit vectors.
%
%   set = flpm.defectSet(size, depth, position)
%   set = flpm.defectSet(size, depth, position, 'Axial', y, 'Population', p)
%
%     size      sqrt(area) of each defect, mm
%     depth     depth of the defect centre below the free surface, mm
%     position  position of the centre along the free surface, mm
%               (only relative positions matter; scalar = all at one place)
%     'Axial'   position along the load axis, mm (default: all in one plane)
%   NaN entries (padding of variable-length rows) are removed.
ip = inputParser;
ip.addParameter('Population', []);
ip.addParameter('Axial', []);
ip.parse(varargin{:});
n = numel(sizeMm);
v = {sizeMm(:), depthMm(:) .* ones(n, 1), positionMm(:) .* ones(n, 1)};
if any(cellfun(@numel, v) ~= n)
    error('flpm:defectSet:length', 'All vectors must have the same length.');
end
keep = ~any(isnan([v{:}]), 2);
if ~isempty(ip.Results.Axial)
    keep = keep & ~isnan(ip.Results.Axial(:) .* ones(n, 1));
end
set = struct('size', v{1}(keep), 'depth', v{2}(keep), 'position', v{3}(keep));
if ~isempty(ip.Results.Axial)
    y = ip.Results.Axial(:) .* ones(n, 1);
    set.axial = y(keep);
end
if ~isempty(ip.Results.Population)
    p = ip.Results.Population(:) .* ones(n, 1);
    set.population = p(keep);
end
end

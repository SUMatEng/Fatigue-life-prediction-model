function S = prepareState(defects, geom)
%PREPARESTATE  Convert a defect set (mm, public API) to the SI solver state.
%   defects: struct with column (or row) vectors
%       size    sqrt(area) of each defect, mm
%       radius  radial position of the defect centre from the bar axis, mm
%       angle   angular position, rad
%       aspect  projected aspect ratio a/c (<= 1)
%     optional
%       yFloor      lower bound on Y once the defect breaks the surface
%                   (roughness shielding factor). Default geom.surfaceY.
%       population  population index (bookkeeping only)
n = numel(defects.size);
req = {'size', 'radius', 'angle', 'aspect'};
for k = 1:numel(req)
    if ~isfield(defects, req{k}) || numel(defects.(req{k})) ~= n
        error('flpm:defects:field', ...
            'Defect field "%s" missing or of the wrong length.', req{k});
    end
end
if n == 0
    error('flpm:defects:empty', 'The defect set is empty.');
end
sz = defects.size(:);
if any(~isfinite(sz)) || any(sz <= 0)
    error('flpm:defects:size', 'Defect sizes must be positive and finite.');
end
rad = defects.radius(:);
if any(rad < 0) || any(rad > geom.D/2 + 1e-12)
    error('flpm:defects:radius', ...
        'Radial positions must lie in [0, D/2] = [0, %g] mm.', geom.D/2);
end

S.a = sz / 1000;                              % m
[S.x, S.y] = pol2cart(defects.angle(:), rad / 1000);
S.asp = defects.aspect(:);
S.ext = zeros(n, 1);                          % crack extension (Delta a)
if isfield(defects, 'yFloor') && ~isempty(defects.yFloor)
    S.yfloor = defects.yFloor(:) .* ones(n, 1);
else
    S.yfloor = geom.surfaceY * ones(n, 1);
end
S.id = (1:n)';                                % original index
S.dK = zeros(n, 1);                           % last computed nominal dK
end

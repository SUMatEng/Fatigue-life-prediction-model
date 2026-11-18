function S = prepareState(defects)
%PREPARESTATE  Convert a defect set (mm, public API) to the SI solver state.
%   defects: struct with column (or row) vectors
%       size      sqrt(area) of the whole defect, mm
%       depth     depth of the defect centre below the free surface, mm (>= 0)
%       position  position of the centre along the free surface, mm
%     optional
%       axial     position along the load axis, mm (default 0: one plane)
n = numel(defects.size);
req = {'size', 'depth', 'position'};
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
dep = defects.depth(:);
if any(~isfinite(dep)) || any(dep < 0)
    error('flpm:defects:depth', 'Defect depths must be finite and >= 0.');
end
S.a   = sz / 1000;                 % sqrt(area) of the disc, m
S.z   = dep / 1000;                % centre depth, m
S.s   = defects.position(:) / 1000;
if isfield(defects, 'axial') && ~isempty(defects.axial)
    S.y = defects.axial(:) .* ones(n, 1) / 1000;
else
    S.y = zeros(n, 1);
end
S.ext = zeros(n, 1);               % crack extension Delta a
S.id  = (1:n)';                    % original index
S.dK  = zeros(n, 1);               % last nominal dK
end

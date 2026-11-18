function [S, merged] = coalesce(S, rScale)
%COALESCE  Merge defects whose interaction circles overlap (contact rule).
%   [S, merged] = flpm.internal.coalesce(S, rScale)
%   Two defects interact when the distance between centres is smaller than
%   the sum of their interaction radii r = rScale*a (rScale = 1 reproduces the
%   published code, where r = sqrt(area)). Interacting defects are grouped
%   transitively and each group is replaced by one defect (Acta Mater. 2022,
%   Sec. 2.2; code by T. Becker 2022):
%     size      sqrt(sum a_i^2)            (total equivalent area)
%     position  sum(a_i x_i)/sum(a_i)      (sqrt(area)-weighted centre)
%     aspect    sum(a_i asp_i)/sum(a_i)
%     extension max(Delta a_i)            (closure build-up state)
%     yFloor    max(yFloor_i)
%     id        id of the member with the largest dK in the last increment
%   The merged defect takes the place of that member; the others are
%   removed. Merging is repeated until no overlaps remain.
merged = false;
while numel(S.a) > 1
    r = rScale * S.a;
    dx = S.x - S.x.';
    dy = S.y - S.y.';
    ov = sqrt(dx.^2 + dy.^2) < (r + r.');
    ov(1:size(ov,1)+1:end) = false;
    if ~any(ov(:))
        return
    end
    merged = true;
    comp = components(ov);
    keep = true(numel(S.a), 1);
    for g = unique(comp(:)).'
        m = find(comp == g);
        if numel(m) < 2, continue, end
        [~, k] = max(S.dK(m));          % first maximum, as in the code
        rep = m(k);
        w = S.a(m);
        S.x(rep) = sum(w .* S.x(m)) / sum(w);
        S.y(rep) = sum(w .* S.y(m)) / sum(w);
        S.asp(rep) = sum(w .* S.asp(m)) / sum(w);
        S.a(rep) = sqrt(sum(w.^2));
        S.ext(rep) = max(S.ext(m));
        S.yfloor(rep) = max(S.yfloor(m));
        S.dK(rep) = max(S.dK(m));
        keep(m(m ~= rep)) = false;
    end
    S = select(S, keep);
end
end

function comp = components(A)
% Connected components of a symmetric logical adjacency matrix.
n = size(A, 1);
comp = zeros(n, 1);
c = 0;
for s = 1:n
    if comp(s), continue, end
    c = c + 1;
    stack = s;
    comp(s) = c;
    while ~isempty(stack)
        v = stack(end); stack(end) = [];
        nb = find(A(:, v) & comp == 0);
        comp(nb) = c;
        stack = [stack; nb]; %#ok<AGROW>
    end
end
end

function S = select(S, keep)
f = {'a', 'x', 'y', 'asp', 'ext', 'yfloor', 'id', 'dK'};
for k = 1:numel(f)
    S.(f{k}) = S.(f{k})(keep);
end
end

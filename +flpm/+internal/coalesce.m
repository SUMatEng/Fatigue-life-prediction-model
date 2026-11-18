function [S, merged] = coalesce(S)
%COALESCE  Merge defects whose discs touch or overlap.
%   [S, merged] = flpm.internal.coalesce(S)
%   Each defect is a disc of radius r = sqrt(area)/sqrt(pi) centred at
%   (position s, depth z). Discs that touch are grouped transitively and
%   each group is replaced by one disc with
%     area      sum of the member areas (sqrt(area) = sqrt(sum a_i^2))
%     centre    area-weighted centroid
%     extension largest member extension (closure build-up state)
%     id        id of the member with the largest dK
%   This is repeated until no discs touch.
merged = false;
while numel(S.a) > 1
    r = S.a / sqrt(pi);
    ov = hypot(S.s - S.s.', S.z - S.z.') <= (r + r.');
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
        [~, k] = max(S.dK(m));
        rep = m(k);
        w = S.a(m).^2;
        S.s(rep) = sum(w .* S.s(m)) / sum(w);
        S.z(rep) = max(sum(w .* S.z(m)) / sum(w), 0);
        S.a(rep) = sqrt(sum(w));
        S.ext(rep) = max(S.ext(m));
        S.dK(rep) = max(S.dK(m));
        keep(m(m ~= rep)) = false;
    end
    f = fieldnames(S);
    for q = 1:numel(f)
        S.(f{q}) = S.(f{q})(keep);
    end
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

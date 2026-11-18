function [S, merged] = coalesce(S)
%COALESCE  Flaw interaction after BS 7910 and re-characterisation.
%   [S, merged] = flpm.internal.coalesce(S)
%
%   Every defect is idealised as a circular planar flaw of radius
%   r = sqrt(area)/sqrt(pi) in the plane normal to the load through its
%   centre (axial coordinate S.y), centred at (position S.s, depth S.z).
%   For each pair (see docs/MODEL.md for sources and caveats):
%
%   1. Alignment (non-coplanar flaws). Flaws in parallel planes a distance
%      H = |y1 - y2| apart are treated as coplanar if the direct distance
%      between their nearest tips, D = sqrt(H^2 + g^2) with g the in-plane
%      gap, satisfies D <= 0.5 (l1 + l2), l = 2r the flaw length.
%   2. Combination (coplanar flaws), with the gaps between their containing
%      rectangles along the surface (sx) and in depth (sz):
%        overlapping rectangles            interact
%        side by side (sz <= 0)            interact if sx <= min(l1, l2)
%        stacked in depth (sx <= 0)        interact if sz <= a1 + a2 (a = r)
%        otherwise                         no interaction
%   3. Re-characterisation. Interacting flaws (grouped transitively) are
%      replaced by the rectangle containing all of them; the new flaw is the
%      ellipse inscribed in it, held as the disc of the same area
%      (sqrt(area) = sqrt(pi a c)) at the rectangle centre and the
%      area-weighted axial position. The merged crack takes the largest
%      member extension and the id of the member with the largest dK.
merged = false;
while numel(S.a) > 1
    A = interacting(S);
    if ~any(A(:))
        return
    end
    merged = true;
    comp = components(A);
    keep = true(numel(S.a), 1);
    for g = unique(comp(:)).'
        m = find(comp == g);
        if numel(m) < 2, continue, end
        [~, k] = max(S.dK(m));
        rep = m(k);
        r = S.a(m) / sqrt(pi);
        s0 = min(S.s(m) - r); s1 = max(S.s(m) + r);
        z0 = min(S.z(m) - r); z1 = max(S.z(m) + r);
        c = (s1 - s0) / 2;  h = (z1 - z0) / 2;
        w = S.a(m).^2;
        S.s(rep) = (s0 + s1) / 2;
        S.z(rep) = max((z0 + z1) / 2, 0);
        S.y(rep) = sum(w .* S.y(m)) / sum(w);
        S.a(rep) = sqrt(pi * c * h);
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

function A = interacting(S)
% Pairwise BS 7910 alignment + combination test.
r = S.a / sqrt(pi);
rr = r + r.';
H = abs(S.y - S.y.');
dIn = hypot(S.s - S.s.', S.z - S.z.');
g = max(dIn - rr, 0);
aligned = H == 0 | hypot(H, g) <= rr;               % D <= 0.5 (l1 + l2)
sx = abs(S.s - S.s.') - rr;                        % rectangle gaps
sz = abs(S.z - S.z.') - rr;
lmin = 2 * min(r, r.');
A = aligned & ((sx <= 0 & sz <= 0) | (sz <= 0 & sx <= lmin) | (sx <= 0 & sz <= rr));
A(1:size(A,1)+1:end) = false;
end

function comp = components(A)
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

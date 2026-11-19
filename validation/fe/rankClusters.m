function T = rankClusters(xyz, sqrtArea)
%RANKCLUSTERS  Kachanov interaction factor of every pore in a large XCT set.
%   T = rankClusters(xyz, sqrtArea)
%   xyz      N x 3 centres [x y z] in mm, with z along the load axis
%   sqrtArea N x 1 sqrt of the area projected on the plane normal to z, mm
%   Pores are grouped into connected components of pairs closer than
%   4 (r1 + r2) between edges (the cut-off of the interaction factor),
%   and flpm.internal.interactionFactor is evaluated per component, so
%   that sets of 1e4-1e5 pores do not need an N x N distance matrix.
%   T: table with gamma, the nearest-neighbour edge gap/(r1 + r2), the
%   axial offset/r to it and the component size, sorted by gamma.
n = numel(sqrtArea);
r = sqrtArea(:) / sqrt(pi);
rmax = max(r);
[idx, d] = rangesearch(xyz, xyz, 10 * rmax);
I = []; J = []; gap = inf(n, 1); nn = zeros(n, 1);
for i = 1:n
    j = idx{i}(2:end); dj = d{i}(2:end);
    if isempty(j), continue, end
    e = dj(:) - r(i) - r(j);
    rel = e ./ (r(i) + r(j));
    [m, k] = min(rel);
    gap(i) = m; nn(i) = j(k);
    keep = e < 4 * (r(i) + r(j));
    I = [I; i * ones(nnz(keep), 1)]; J = [J; j(keep)']; %#ok<AGROW>
end
G = graph(I, J, [], n);
comp = conncomp(G)';
gam = ones(n, 1);
cs = accumarray(comp, 1);
for c = find(cs > 1)'
    k = find(comp == c);
    if numel(k) < 2, continue, end
    S = struct('a', sqrtArea(k) / 1e3, 's', xyz(k, 1) / 1e3, 'z', xyz(k, 2) / 1e3, 'y', xyz(k, 3) / 1e3);
    gam(k) = flpm.internal.interactionFactor(S);
end
axialOffset = NaN(n, 1);
h = nn > 0;
axialOffset(h) = abs(xyz(h, 3) - xyz(nn(h), 3)) ./ r(h);
T = table((1:n)', sqrtArea(:), gam, gap, nn, axialOffset, cs(comp), ...
    'VariableNames', {'id', 'sqrtArea', 'gamma', 'gapRel', 'neighbour', 'axialOffsetRel', 'componentSize'});
T = sortrows(T, 'gamma', 'descend');
end

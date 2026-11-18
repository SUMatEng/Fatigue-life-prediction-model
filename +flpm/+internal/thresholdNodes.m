function [mats, w] = thresholdNodes(mat, n)
%THRESHOLDNODES  Gauss-Hermite nodes of a lognormal dKth_lc scatter.
%   mats{k} is mat with dKth_lc at node k (median = mat.dKth_lc), w the
%   weights; a single node if mat.dKth_cov is zero.
if nargin < 2, n = 7; end
if ~isfield(mat, 'dKth_cov') || mat.dKth_cov <= 0
    mats = {mat}; w = 1; return
end
s = sqrt(log(1 + mat.dKth_cov^2));
% Golub-Welsch for probabilists' Hermite polynomials
J = diag(sqrt(1:n-1), 1) + diag(sqrt(1:n-1), -1);
[V, D] = eig(J);
z = diag(D);
w = V(1, :).'.^2;
mats = cell(n, 1);
for k = 1:n
    mats{k} = drawn(mat, exp(s * z(k)));
end
end

function m = drawn(mat, f)
m = mat;
m.dKth_lc = max(mat.dKth_lc * f, mat.dKth_eff * (1 + 1e-6));
m.dKth_cov = 0;
end

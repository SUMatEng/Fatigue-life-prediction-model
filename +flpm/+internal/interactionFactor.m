function [gam, info] = interactionFactor(S, opts)
%INTERACTIONFACTOR  Continuous K amplification of interacting parallel cracks.
%   gam = interactionFactor(S) gives, for every defect of the solver state
%   S (sqrt(area) a, in-plane centre (s, z), axial position y), the ratio
%   of its largest stress intensity along the front to that of the same
%   crack alone. Each defect is the equivalent penny crack of radius
%   r = a/sqrt(pi) in the plane normal to the load.
%
%   Method (Kachanov, Int. J. Solids Struct. 23 (1987) 23; Kachanov &
%   Laures, Int. J. Fract. 41 (1989) 289): every crack i carries the
%   remote stress plus the stress induced by the others, each loaded by
%   its average traction p_j:
%       p_i = 1 + sum_j Lambda_ij p_j,   Lambda_ij = <sigma_j>_i
%   where <sigma_j>_i is the average over crack i of the normal stress of
%   crack j under unit pressure (flpm.internal.pennyStress). The resulting
%   non-uniform traction t(P) on crack i gives K along its front by the
%   exact weight function of the penny crack (Galin; Kassir & Sih),
%       K(Q) = 1/(pi sqrt(pi r)) int t(P) sqrt(r^2 - rho_P^2)/|PQ|^2 dA,
%   evaluated at 24 front points; gam is max_Q K(Q)/K_uniform(Q).
%   Pairs farther apart than 4 (r_i + r_j) between edges are taken as
%   non-interacting (amplification below about 1 %). Parallel cracks
%   offset along the load axis shield or amplify according to the sign
%   of the induced stress. The free surface is not included in the
%   interaction (it enters through Murakami's Y and A).
n = numel(S.a);
gam = ones(n, 1);
info = struct('pairs', 0);
if n < 2, return, end
if nargin < 2 || isempty(opts), opts = struct(); end
r = S.a / sqrt(pi);
D = sqrt((S.s - S.s.').^2 + (S.z - S.z.').^2 + (S.y - S.y.').^2);
near = D - (r + r.') < 4 * (r + r.') & ~eye(n);
if ~any(near(:)), return, end
info.pairs = nnz(near) / 2;

[xq, wq, xf] = quadrature();
idx = find(any(near, 2));
% Lambda over the cracks that have a neighbour (others: p = 1)
m = numel(idx);
Lam = zeros(m);
P = cell(m, 1);
for a = 1:m
    i = idx(a);
    P{a} = [S.s(i) + r(i) * xq(:, 1), S.z(i) + r(i) * xq(:, 2)];   % points on crack i
    for b = 1:m
        j = idx(b);
        if a == b || ~near(i, j), continue, end
        sig = induced(P{a}, S, j, r(j), S.y(i));
        Lam(a, b) = wq.' * sig;
    end
end
p = (eye(m) - Lam) \ ones(m, 1);
for a = 1:m
    i = idx(a);
    t = ones(size(xq, 1), 1);
    for b = 1:m
        j = idx(b);
        if a == b || ~near(i, j), continue, end
        t = t + p(b) * induced(P{a}, S, j, r(j), S.y(i));
    end
    % weight function at the front points Q (unit circle)
    dx = xq(:, 1) - xf(:, 1).';  dy = xq(:, 2) - xf(:, 2).';
    Wt = sqrt(max(1 - sum(xq.^2, 2), 0)) .* wq ./ (dx.^2 + dy.^2);
    gam(i) = max((t.' * Wt) ./ sum(Wt, 1));
end
end

function s = induced(Pts, S, j, rj, yi)
% normal stress on points Pts (plane y = yi) of crack j under unit pressure
rho = hypot(Pts(:, 1) - S.s(j), Pts(:, 2) - S.z(j)) / rj;
h = abs(yi - S.y(j)) / rj;
s = flpm.internal.pennyStress(rho, h * ones(size(rho)));
s(~isfinite(s)) = 0;
end

function [xq, wq, xf] = quadrature()
persistent X W F
if isempty(X)
    nr = 16; nt = 32;
    [u, wu] = gaussLegendre(nr);              % on [0, 1]
    rr = sqrt(u);                              % area-uniform in radius
    th = (0:nt-1)' / nt * 2 * pi;
    [R, T] = ndgrid(rr, th);
    X = [R(:) .* cos(T(:)), R(:) .* sin(T(:))];
    Wm = repmat(wu, 1, nt) / nt;               % area fractions, sum = 1
    W = Wm(:);
    tf = (0:23)' / 24 * 2 * pi;
    F = [cos(tf), sin(tf)];
end
xq = X; wq = W; xf = F;
end

function [x, w] = gaussLegendre(n)
k = 1:n-1;
b = k ./ sqrt(4 * k.^2 - 1);
[V, Dg] = eig(diag(b, 1) + diag(b, -1));
x = (diag(Dg) + 1) / 2;
w = (2 * V(1, :).'.^2) / 2;
end

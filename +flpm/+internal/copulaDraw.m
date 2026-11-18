function U = copulaDraw(cop, n, rs)
%COPULADRAW  n draws (n x 2) of the size-depth copula from stream rs.
R = [1 cop.rho; cop.rho 1];
Z = randn(rs, n, 2) * chol(R);
switch cop.family
    case 'gaussian'
        U = 0.5 * erfc(-Z / sqrt(2));
    case 't'
        W = 2 * randg(rs, cop.nu / 2, n, 1) / cop.nu;      % chi2(nu)/nu
        U = tcdf(Z ./ sqrt(W), cop.nu);
end
U = min(max(U, eps), 1 - eps);
end

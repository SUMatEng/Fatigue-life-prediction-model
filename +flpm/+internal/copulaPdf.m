function c = copulaPdf(cop, U)
%COPULAPDF  Density of the size-depth copula at U (n x 2).
switch cop.family
    case 'gaussian'
        c = copulapdf('Gaussian', U, cop.rho);
    case 't'
        c = copulapdf('t', U, cop.rho, cop.nu);
end
end

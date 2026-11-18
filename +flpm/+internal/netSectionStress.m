function S = netSectionStress(a, S0, D, e)
%NETSECTIONSTRESS  Maximum net-section stress of a bar with a surface crack.
%   S = flpm.internal.netSectionStress(a, S0, D, e) for crack depth a, remote
%   stress S0, bar diameter D (same length unit as a) and load-line
%   eccentricity factor e. Ported from net_stress.m (published code): the
%   cracked segment is bounded by a circular arc, and the net-section stress
%   combines the axial term T/A and the bending term from the shift of the
%   neutral axis. Valid for a < D/2. Scalar a.
c = (D/2)*atan((a*(D - a))/(D*((D/2) - a)));
R = D/2;
A0 = pi*R^2;
alpha = c/R;
R2 = R*tan(alpha);
beta = pi/2 - alpha;
A1 = (alpha - sin(2*alpha)/2)*R^2;
A2 = (beta - sin(2*beta)/2)*R2^2;
A = A0 - A1 - A2;
if alpha < 1e-6           % uncracked limit
    S = S0; return
end
l1 = (4*R*sin(alpha)^3)/(3*(2*alpha - sin(2*alpha)));
l2 = (4*R2*sin(beta)^3)/(3*(2*beta - sin(2*beta)));
x1 = -l1;
x2 = -(R*sec(alpha) - l2);
qx = -((A1*x1) + (A2*x2))/A;
I0 = (pi/4)*R^4 + A0*qx^2;
I1 = (R^4/4)*(alpha - sin(4*alpha)/4 - (16*sin(alpha)^6)/(9*(alpha - sin(alpha)*cos(alpha)))) + A1*(qx - x1)^2;
I2 = (R2^4/4)*(beta - sin(4*beta)/4 - (16*sin(beta)^6)/(9*(beta - sin(beta)*cos(beta)))) + A2*(qx - x2)^2;
Iy = I0 - I1 - I2;
T = S0*A0;
My = S0*A0*qx;
xa = -(R*cos(alpha) + qx) + e*R*(1 + cos(alpha));
xb = (R - qx) - e*R*(1 + cos(alpha));
SnA = T/A - (My/Iy)*xa;
SnB = T/A - (My/Iy)*xb;
S = max(abs(SnA), abs(SnB));
end

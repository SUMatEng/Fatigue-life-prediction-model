classdef FlpmFiniteLifeTest < matlab.unittest.TestCase
    %FLPMFINITELIFETEST  Growth-law calibration and the sealed-crack environment.

    methods (Test)
        function fitRecoversLongCrackLaw(tc)
            mat = flpm.material('Ti64-AB');
            U = flpm.internal.newmanU(0.1, mat.alpha);
            dK = logspace(log10(2.9), log10(20), 30)';
            q = [2e-10 2.9 0.4 2.6];
            y = q(1) * U^q(2) * dK.^(q(2) - q(3)) .* (dK - q(4)).^q(3);
            rs = RandStream('mt19937ar', 'Seed', 3);
            y = y .* exp(0.05 * randn(rs, size(y)));
            G = flpm.fitGrowthLaw(dK, y, 'LoadRatio', 0.1, 'Material', mat);
            tc.verifyEqual(G.dKth, q(4), 'RelTol', 0.05);
            tc.verifyEqual(G.m, q(2), 'RelTol', 0.1);
            tc.verifyLessThan(G.rmsLog, 0.05);
            tc.verifyEqual(G.rate(10, 0.1), q(1)*U^q(2)*10^(q(2)-q(3))*(10-q(4))^q(3), 'RelTol', 0.15);
        end

        function fitWithoutClosureForSealedCracks(tc)
            dK = linspace(2.2, 8, 25)';
            y = 5e-12 * dK.^(2.5) .* (dK - 2).^0.35;
            G = flpm.fitGrowthLaw(dK, y, 'Closure', false, 'FixExponents', true, ...
                'Material', flpm.material('Ti64-AB', 'm', 2.85, 'p', 0.35));
            tc.verifyEqual(G.dKth, 2, 'RelTol', 1e-3);
            tc.verifyEqual(G.C1, 5e-12, 'RelTol', 1e-3);
            tc.verifyFalse(G.sealed.closure);
            tc.verifyEqual(G.sealed.dKth_lc, G.sealed.dKth_eff);
        end

        function sealedEnvironmentActsOnInternalCracksOnly(tc)
            mat = flpm.material('Ti64-AB');
            slow = mat; slow.sealed = struct('C1', mat.C1 / 10);
            same = mat; same.sealed = struct('C1', mat.C1);
            dInt = flpm.defectSet(0.08, 0.4, 0);          % internal
            dSur = flpm.defectSet(0.08, 0.0, 0);          % surface-breaking
            s = 1.3 * flpm.fatigueStrength(dInt, mat).stress;
            n0 = flpm.fatigueLife(s, dInt, mat).N;
            tc.verifyEqual(flpm.fatigueLife(s, dInt, same).N, n0, 'RelTol', 1e-12);
            n1 = flpm.fatigueLife(s, dInt, slow).N;
            tc.verifyGreaterThan(n1, 2 * n0);             % slower until breakthrough
            ss = 1.3 * flpm.fatigueStrength(dSur, mat).stress;
            tc.verifyEqual(flpm.fatigueLife(ss, dSur, slow).N, flpm.fatigueLife(ss, dSur, mat).N, 'RelTol', 1e-12);
        end
    end
end

classdef FlpmSolverTest < matlab.unittest.TestCase
    %FLPMSOLVERTEST  Crack-growth solver and fatigue-strength search.

    properties
        mat
        single      % one near-surface pore
        multi       % five pores, two of which touch
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.mat = flpm.material('Ti64-AB');
            tc.single = flpm.defectSet(0.056, 0.256, 0);
            tc.multi = flpm.defectSet([0.06 0.045 0.07 0.05 0.055], ...
                [0.25 0.20 1.65 0.15 0.23], [0 0.05 5 9 0.12]);
        end
    end

    methods (Test)
        function stepSizeConvergence(tc)
            % halving Eta changes the life by < 0.5 %
            for s = [300 400]
                a = flpm.fatigueLife(s, tc.multi, tc.mat, flpm.options('Eta', 0.005));
                b = flpm.fatigueLife(s, tc.multi, tc.mat, flpm.options('Eta', 0.0025));
                if isfinite(b.N)
                    tc.verifyEqual(a.N, b.N, 'RelTol', 5e-3);
                end
            end
        end

        function noGrowthBelowOnset(tc)
            r = flpm.fatigueLife(50, tc.single, tc.mat);
            tc.verifyEqual(r.N, Inf);
            tc.verifyEqual(r.status, 'noGrowth');
        end

        function lifeDecreasesWithStress(tc)
            N = arrayfun(@(s) flpm.fatigueLife(s, tc.multi, tc.mat).N, [350 400 500 600]);
            tc.verifyTrue(all(diff(N(isfinite(N))) < 0));
        end

        function strengthIsMonotonicInSize(tc)
            sizes = [0.04 0.06 0.09];
            s = arrayfun(@(a) flpm.fatigueStrength(flpm.defectSet(a, 0.4, 0), tc.mat).stress, sizes);
            tc.verifyTrue(all(diff(s) < 0));
        end

        function nearSurfaceDefectIsMostSevere(tc)
            % Murakami: a pore whose ligament is within the proximity limit
            % acts as a larger surface defect, so it is more severe than the
            % same pore deep inside or cut by the surface (less area inside).
            % sqrt(area) = 0.05 mm: r = 0.028 mm, proximity limit z = 0.063 mm
            depths = [0.01 0.06 0.5 2];
            s = arrayfun(@(z) flpm.fatigueStrength(flpm.defectSet(0.05, z, 0), tc.mat).stress, depths);
            tc.verifyLessThan(s(2), s(1));
            tc.verifyLessThan(s(2), s(3));
            tc.verifyLessThanOrEqual(s(3), s(4) + 1e-9);
        end

        function bracketSolutionHitsTarget(tc)
            r = flpm.fatigueStrength(tc.multi, tc.mat);
            tc.verifyTrue(r.converged);
            lo = flpm.fatigueLife(r.stress - 0.2, tc.multi, tc.mat);
            if r.stepInSN
                tc.verifyGreaterThan(lo.N, 1e7);
            elseif ~r.onsetLimited
                tc.verifyEqual(r.N, 1e7, 'RelTol', 0.02);
                tc.verifyGreaterThan(lo.N, r.N);
            end
        end

        function failureByToughness(tc)
            r = flpm.fatigueLife(600, tc.multi, tc.mat);
            tc.verifyEqual(r.status, 'Kc');
            tc.verifyTrue(isfinite(r.N));
        end

        function monteCarloSmoke(tc)
            pop = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Weibull', 'a', 0.06, 'b', 5), ...
                'Depth', makedist('Exponential', 'mu', 0.4), 'Count', 6);
            T = flpm.monteCarlo(pop, tc.mat, 'SurfaceLength', 22.9, ...
                'NumSimulations', 4, 'Seed', 2);
            tc.verifySize(T, [4 14]);
            tc.verifyTrue(all(T.stress > 50 & T.stress < 600));
            T2 = flpm.monteCarlo(pop, tc.mat, 'SurfaceLength', 22.9, ...
                'NumSimulations', 4, 'Seed', 2);
            tc.verifyEqual(T2.stress, T.stress);
        end

        function emptyThresholdedSetGivesDefectFreeStrength(tc)
            pop = flpm.defectPopulation('volumetric', 'Size', 0.01, ...
                'Depth', 1, 'Count', 3, 'Threshold', 0.05);
            T = flpm.monteCarlo(pop, tc.mat, 'SurfaceLength', 10, 'NumSimulations', 2);
            tc.verifyEqual(T.stress, [450; 450]);
            tc.verifyEqual(T.status, ["noDefects"; "noDefects"]);
        end

        function inputValidation(tc)
            bad = flpm.defectSet(0.05, -0.1, 0);
            tc.verifyError(@() flpm.fatigueLife(300, bad, tc.mat), 'flpm:defects:depth');
            tc.verifyError(@() flpm.fatigueLife(-1, tc.single, tc.mat), ...
                'flpm:fatigueLife:stress');
        end
    end
end

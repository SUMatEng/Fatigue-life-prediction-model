classdef FlpmSolverTest < matlab.unittest.TestCase
    %FLPMSOLVERTEST  Crack-growth solver and fatigue-strength search.
    %
    %   Reference values marked "original code" were produced on 2026-11-18 by
    %   running N. Macallister's Nasgro.m / multi_coalition_version3.m
    %   (Feb 2025, 'Nic''s Model') with dKth,lc = 2.7 MPa*sqrt(m) and
    %   flow stress 1150 MPa (Mater. Des. 2026, Table 2), i.e. the parameters
    %   of flpm.material('Ti64-AB').

    properties
        mat
        geom
        single      % defect of the sanbox_multi single-defect case
        multi       % five defects, two of which coalesce
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.mat = flpm.material('Ti64-AB');
            tc.geom = flpm.geometry();
            tc.single = flpm.defectSet(0.055997925, 3.394247128, 2.544400909, 1);
            tc.multi = flpm.defectSet([0.06 0.045 0.07 0.05 0.055], ...
                [3.40 3.45 2.0 3.50 3.42], [0 0.03 1.5 -2 0.05], ...
                [0.6 0.8 0.5 1 0.7]);
        end
    end

    methods (Test)
        function legacyReproducesOriginalCycleCounts(tc)
            o = flpm.modelOptions('legacy');
            ref = [260 61200; 300 36700; 350 22000];          % original code
            for k = 1:size(ref, 1)
                r = flpm.fatigueLife(ref(k,1), tc.multi, tc.mat, tc.geom, o);
                tc.verifyEqual(r.N, ref(k,2), sprintf('multi, %g MPa', ref(k,1)));
                tc.verifyEqual(r.killer, 5);
            end
            r = flpm.fatigueLife(350, tc.single, tc.mat, tc.geom, o);
            tc.verifyEqual(r.N, 89900);                         % original code
            r = flpm.fatigueLife(300, tc.single, tc.mat, tc.geom, o);
            tc.verifyEqual(r.N, Inf);                           % run-out
        end

        function legacyStrengthMatchesOriginalFminsearch(tc)
            % original code: fminsearch from 400 MPa -> 315.703125 MPa
            r = flpm.fatigueStrength(tc.single, tc.mat, tc.geom, flpm.modelOptions('legacy'));
            tc.verifyEqual(r.stress, 315.703125, 'AbsTol', 1e-9);
            tc.verifyTrue(r.converged);
        end

        function bracketAgreesWithLegacySolver(tc)
            o = flpm.modelOptions('legacy', 'Stepping', 'adaptive', 'Eta', 0.005, ...
                'Solver', 'bracket', 'MaxLife', 1.02e7);
            r = flpm.fatigueStrength(tc.single, tc.mat, tc.geom, o);
            tc.verifyTrue(r.converged);
            tc.verifyEqual(r.stress, 315.703125, 'AbsTol', 0.5);
            tc.verifyEqual(r.N, 1e7, 'RelTol', 0.02);
        end

        function adaptiveStepConverges(tc)
            % Away from the threshold the adaptive scheme matches fine fixed
            % stepping to better than 0.5 %.
            of = flpm.modelOptions('paper', 'Stepping', 'fixed', 'FixedStep', 25);
            oa = flpm.modelOptions('paper', 'Eta', 0.005);
            for s = [300 400]
                rf = flpm.fatigueLife(s, tc.multi, tc.mat, tc.geom, of);
                ra = flpm.fatigueLife(s, tc.multi, tc.mat, tc.geom, oa);
                tc.verifyEqual(ra.N, rf.N, 'RelTol', 5e-3);
            end
            % long life: adaptive needs far fewer steps than fixed-100
            o = flpm.modelOptions('paper', 'Eta', 0.005);
            rf = flpm.fatigueLife(330, tc.single, tc.mat, tc.geom, ...
                flpm.modelOptions('paper', 'Stepping', 'fixed', 'FixedStep', 100));
            ra = flpm.fatigueLife(330, tc.single, tc.mat, tc.geom, o);
            if isfinite(rf.N)
                tc.verifyEqual(ra.N, rf.N, 'RelTol', 0.02);
                tc.verifyLessThan(ra.nSteps, rf.nSteps / 5);
            end
        end

        function noGrowthBelowOnset(tc)
            r = flpm.fatigueLife(50, tc.single, tc.mat, tc.geom, flpm.modelOptions());
            tc.verifyEqual(r.N, Inf);
            tc.verifyEqual(r.status, 'noGrowth');
        end

        function strengthIsMonotonicInSize(tc)
            o = flpm.modelOptions();
            s = zeros(1, 3);
            sizes = [0.04 0.06 0.09];
            for k = 1:3
                d = flpm.defectSet(sizes(k), 3.65 - sizes(k), 0, 1);   % surface
                s(k) = flpm.fatigueStrength(d, tc.mat, tc.geom, o).stress;
            end
            tc.verifyTrue(all(diff(s) < 0));
        end

        function surfaceDefectIsMoreSevereThanCentre(tc)
            % Acta Mater. 2022, Fig. 9: surface initiation lowers the
            % strength relative to an equal defect at the centre.
            o = flpm.modelOptions();
            sc = flpm.fatigueStrength(flpm.defectSet(0.05, 0, 0, 1), tc.mat, tc.geom, o);
            ss = flpm.fatigueStrength(flpm.defectSet(0.05, 3.60, 0, 1), tc.mat, tc.geom, o);
            tc.verifyLessThan(ss.stress, sc.stress);
        end

        function bracketSolutionHitsTarget(tc)
            r = flpm.fatigueStrength(tc.multi, tc.mat, tc.geom, flpm.modelOptions());
            tc.verifyTrue(r.converged);
            lo = flpm.fatigueLife(r.stress - 0.2, tc.multi, tc.mat, tc.geom, flpm.modelOptions());
            if r.stepInSN
                % N_f jumps at the strength: arrest just below it
                tc.verifyGreaterThan(lo.N, 1e7);
            elseif ~r.onsetLimited
                tc.verifyEqual(r.N, 1e7, 'RelTol', 0.02);
                tc.verifyGreaterThan(lo.N, r.N);
            end
        end

        function plasticZoneFactorMatters(tc)
            % pi/6 (code) vs 1/(6 pi) (Acta Eq. 12): larger rp -> shorter life
            a = flpm.fatigueLife(300, tc.multi, tc.mat, tc.geom, ...
                flpm.modelOptions('paper', 'PlasticZoneFactor', pi/6));
            b = flpm.fatigueLife(300, tc.multi, tc.mat, tc.geom, flpm.modelOptions('paper'));
            tc.verifyLessThan(a.N, b.N);
        end

        function monteCarloSmoke(tc)
            pop = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Weibull', 'a', 0.06, 'b', 5), ...
                'Depth', makedist('Exponential', 'mu', 0.4), ...
                'Aspect', makedist('Uniform', 'lower', 0.4, 'upper', 1), ...
                'Count', 6);
            T = flpm.monteCarlo(pop, tc.mat, tc.geom, 'NumSimulations', 4, 'Seed', 2);
            tc.verifySize(T, [4 16]);
            tc.verifyTrue(all(isfinite(T.stress)));
            tc.verifyTrue(all(T.stress > 50 & T.stress < 600));
            T2 = flpm.monteCarlo(pop, tc.mat, tc.geom, 'NumSimulations', 4, 'Seed', 2);
            tc.verifyEqual(T2.stress, T.stress);
        end

        function emptyThresholdedSetGivesDefectFreeStrength(tc)
            pop = flpm.defectPopulation('volumetric', 'Size', 0.01, ...
                'Radius', 1, 'Count', 3, 'Threshold', 0.05);
            T = flpm.monteCarlo(pop, tc.mat, tc.geom, 'NumSimulations', 2);
            tc.verifyEqual(T.stress, [450; 450]);
            tc.verifyEqual(T.status, ["noDefects"; "noDefects"]);
        end

        function inputValidation(tc)
            bad = flpm.defectSet(0.05, 4.0, 0, 1);              % outside the bar
            tc.verifyError(@() flpm.fatigueLife(300, bad, tc.mat, tc.geom), ...
                'flpm:defects:radius');
            tc.verifyError(@() flpm.fatigueLife(-1, tc.single, tc.mat, tc.geom), ...
                'flpm:fatigueLife:stress');
        end
    end
end

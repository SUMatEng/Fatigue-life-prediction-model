classdef FlpmSamplingTest < matlab.unittest.TestCase
    %FLPMSAMPLINGTEST  Population descriptors and Monte Carlo sampling.

    properties
        vol
        surf
        L = 22.9
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.vol = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Weibull', 'a', 0.055, 'b', 5), ...
                'Depth', makedist('Exponential', 'mu', 0.3), ...
                'Count', 12, 'Name', 'pores');
            tc.surf = flpm.defectPopulation('surface', ...
                'Size', makedist('Exponential', 'mu', 0.008), ...
                'SizeOffset', 0.020, 'CountPerLength', 0.5, ...
                'ShieldingFactor', 1.02, 'Name', 'roughness');
        end
    end

    methods (Test)
        function reproducibleAndOrderIndependent(tc)
            a = flpm.samplePopulation(tc.vol, 5, 'SurfaceLength', tc.L, 'Seed', 7);
            b = flpm.samplePopulation(tc.vol, 5, 'SurfaceLength', tc.L, 'Seed', 7);
            tc.verifyEqual(a, b);
            c = flpm.samplePopulation(tc.vol, 'simulation', 4, 'SurfaceLength', tc.L, 'Seed', 7);
            tc.verifyEqual(c, a{4});
            d = flpm.samplePopulation(tc.vol, 5, 'SurfaceLength', tc.L, 'Seed', 8);
            tc.verifyNotEqual(a{1}.size, d{1}.size);
        end

        function boundsAndPositions(tc)
            sets = flpm.samplePopulation(tc.vol, 200, 'SurfaceLength', tc.L, 'Seed', 2);
            Z = cell2mat(cellfun(@(x) x.depth, sets, 'UniformOutput', false));
            P = cell2mat(cellfun(@(x) x.position, sets, 'UniformOutput', false));
            S = cell2mat(cellfun(@(x) x.size, sets, 'UniformOutput', false));
            tc.verifyNumElements(sets{1}.size, 12);
            tc.verifyTrue(all(Z >= 0));
            tc.verifyTrue(all(P >= 0 & P < tc.L));
            tc.verifyTrue(all(S > 0));
        end

        function surfaceValleysGiveShieldedK(tc)
            s = flpm.samplePopulation(tc.surf, 'simulation', 1, 'SurfaceLength', tc.L, 'Seed', 3);
            tc.verifyNumElements(s.size, round(0.5*tc.L));
            tc.verifyEqual(s.depth, zeros(numel(s.size), 1));
            % K of the sampled disc equals Y_R s sqrt(pi Rv)
            S = flpm.internal.prepareState(s);
            [Y, A] = flpm.internal.stressIntensity(S);
            Rv = (s.size/sqrt(2)) / (1.02/0.65)^2 / 1000;
            tc.verifyEqual(Y .* sqrt(A), 1.02 * sqrt(Rv), 'RelTol', 1e-12);
            tc.verifyTrue(all(Rv*1000 >= 0.020 - 1e-12));
        end

        function appendedPopulations(tc)
            s = flpm.samplePopulation({tc.vol, tc.surf}, 'simulation', 1, ...
                'SurfaceLength', tc.L, 'Seed', 3);
            tc.verifyEqual(sum(s.population == 1), 12);
            tc.verifyEqual(sum(s.population == 2), round(0.5*tc.L));
        end

        function thresholdThinning(tc)
            p = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Lognormal', 'mu', log(0.05), 'sigma', 0.4), ...
                'Depth', makedist('Uniform', 'lower', 0, 'upper', 3), ...
                'Count', 100, 'Threshold', 0.08);
            s = flpm.samplePopulation(p, 100, 'SurfaceLength', tc.L, 'Seed', 5);
            n = cellfun(@(x) numel(x.size), s);
            Z = cell2mat(cellfun(@(x) x.size, s, 'UniformOutput', false));
            tc.verifyTrue(all(Z > 0.08));
            pExc = 1 - logncdf(0.08, log(0.05), 0.4);
            tc.verifyEqual(mean(n), 100*pExc, 'RelTol', 0.1);
        end

        function dataVectorGivesKernelFit(tc)
            rng(0);
            data = 0.04 + 0.02*rand(300, 1);
            p = flpm.defectPopulation('volumetric', 'Size', data, 'Depth', 0.3, 'Count', 9);
            tc.verifyClass(p.size.value, 'prob.KernelDistribution');
            s = flpm.samplePopulation(p, 'simulation', 1, 'SurfaceLength', tc.L);
            tc.verifyEqual(s.depth, repmat(0.3, 9, 1));
        end

        function countPerLengthScalesWithLength(tc)
            p = flpm.defectPopulation('volumetric', 'Size', 0.05, 'Depth', 1, ...
                'CountPerLength', 0.4);
            s1 = flpm.samplePopulation(p, 'simulation', 1, 'SurfaceLength', 10);
            s2 = flpm.samplePopulation(p, 'simulation', 1, 'SurfaceLength', 40);
            tc.verifyNumElements(s1.size, 4);
            tc.verifyNumElements(s2.size, 16);
        end

        function poissonCount(tc)
            p = flpm.defectPopulation('volumetric', 'Size', 0.05, 'Depth', 1, ...
                'Count', 10, 'CountModel', 'poisson');
            s = flpm.samplePopulation(p, 400, 'SurfaceLength', tc.L, 'Seed', 11);
            n = cellfun(@(x) numel(x.size), s);
            tc.verifyEqual(mean(n), 10, 'AbsTol', 0.6);
            tc.verifyGreaterThan(var(n), 5);
        end

        function invalidInputs(tc)
            tc.verifyError(@() flpm.defectPopulation('volumetric', 'Size', 0.05, 'Count', 3), ...
                'flpm:population:depth');
            tc.verifyError(@() flpm.defectPopulation('surface', 'Size', 0.01), ...
                'flpm:population:count');
            tc.verifyError(@() flpm.defectPopulation('volumetric', 'Depth', 1, 'Count', 3), ...
                'flpm:population:size');
            tc.verifyError(@() flpm.samplePopulation(tc.vol, 2), 'flpm:sample:length');
        end
    end
end

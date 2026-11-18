classdef FlpmSamplingTest < matlab.unittest.TestCase
    %FLPMSAMPLINGTEST  Population descriptors and Monte Carlo sampling.

    properties
        geom
        vol
        surf
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.geom = flpm.geometry('Diameter', 7.3);
            tc.vol = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Weibull', 'a', 0.055, 'b', 5), ...
                'Depth', makedist('Exponential', 'mu', 0.3), ...
                'Aspect', makedist('Beta', 'a', 6, 'b', 2), ...
                'Count', 12, 'Name', 'pores');
            tc.surf = flpm.defectPopulation('surface', ...
                'Size', makedist('Exponential', 'mu', 0.008), ...
                'SizeOffset', 0.020, 'CountPerLength', 0.5, ...
                'ShieldingFactor', 1.02, 'Name', 'roughness');
        end
    end

    methods (Test)
        function reproducibleAndOrderIndependent(tc)
            a = flpm.samplePopulation(tc.vol, tc.geom, 5, 'Seed', 7);
            b = flpm.samplePopulation(tc.vol, tc.geom, 5, 'Seed', 7);
            tc.verifyEqual(a, b);
            % simulation k does not depend on how many are drawn
            c = flpm.samplePopulation(tc.vol, tc.geom, 'simulation', 4, 'Seed', 7);
            tc.verifyEqual(c, a{4});
            d = flpm.samplePopulation(tc.vol, tc.geom, 5, 'Seed', 8);
            tc.verifyNotEqual(a{1}.size, d{1}.size);
        end

        function boundsAreRespected(tc)
            s = flpm.samplePopulation(tc.vol, tc.geom, 50, 'Seed', 1);
            s = s{1};
            sets = flpm.samplePopulation(tc.vol, tc.geom, 200, 'Seed', 2);
            R = cell2mat(cellfun(@(x) x.radius, sets, 'UniformOutput', false));
            A = cell2mat(cellfun(@(x) x.aspect, sets, 'UniformOutput', false));
            Z = cell2mat(cellfun(@(x) x.size, sets, 'UniformOutput', false));
            tc.verifyNumElements(s.size, 12);
            tc.verifyTrue(all(R >= 0 & R <= 3.65));
            tc.verifyTrue(all(A > 0 & A <= 1));
            tc.verifyTrue(all(Z > 0));
        end

        function surfaceCountAndPlacement(tc)
            s = flpm.samplePopulation(tc.surf, tc.geom, 'simulation', 1, 'Seed', 3);
            tc.verifyNumElements(s.size, round(0.5*pi*7.3));
            tc.verifyEqual(s.radius, repmat(3.65, numel(s.size), 1));
            tc.verifyTrue(all(s.size >= 0.020));
            tc.verifyEqual(unique(s.yFloor), 1.02);
            tc.verifyEqual(unique(s.aspect), 1);
        end

        function appendedPopulations(tc)
            s = flpm.samplePopulation({tc.vol, tc.surf}, tc.geom, 'simulation', 1, 'Seed', 3);
            tc.verifyEqual(sum(s.population == 1), 12);
            tc.verifyEqual(sum(s.population == 2), round(0.5*pi*7.3));
            tc.verifyEqual(unique(s.yFloor(s.population == 1)), 1.0335);
        end

        function thresholdThinning(tc)
            p = flpm.defectPopulation('volumetric', ...
                'Size', makedist('Lognormal', 'mu', log(0.05), 'sigma', 0.4), ...
                'Radius', makedist('Uniform', 'lower', 0, 'upper', 3.6), ...
                'Count', 100, 'Threshold', 0.08);
            s = flpm.samplePopulation(p, tc.geom, 100, 'Seed', 5);
            n = cellfun(@(x) numel(x.size), s);
            Z = cell2mat(cellfun(@(x) x.size, s, 'UniformOutput', false));
            tc.verifyTrue(all(Z > 0.08));
            % expected number of exceedances = 100*(1 - F(u))
            pExc = 1 - logncdf(0.08, log(0.05), 0.4);
            tc.verifyEqual(mean(n), 100*pExc, 'RelTol', 0.1);
        end

        function dataVectorGivesKernelFit(tc)
            rng(0);
            data = 0.04 + 0.02*rand(300, 1);
            p = flpm.defectPopulation('volumetric', 'Size', data, ...
                'Radius', 3.0, 'Aspect', 0.8, 'Count', 9);
            tc.verifyClass(p.size.value, 'prob.KernelDistribution');
            s = flpm.samplePopulation(p, tc.geom, 'simulation', 1);
            tc.verifyEqual(s.radius, repmat(3.0, 9, 1));
            tc.verifyEqual(s.aspect, repmat(0.8, 9, 1));
        end

        function poissonCount(tc)
            p = flpm.defectPopulation('volumetric', 'Size', 0.05, ...
                'Radius', makedist('Uniform', 'lower', 0, 'upper', 3.5), ...
                'Count', 10, 'CountModel', 'poisson');
            s = flpm.samplePopulation(p, tc.geom, 400, 'Seed', 11);
            n = cellfun(@(x) numel(x.size), s);
            tc.verifyEqual(mean(n), 10, 'AbsTol', 0.6);
            tc.verifyGreaterThan(var(n), 5);
        end

        function invalidInputs(tc)
            tc.verifyError(@() flpm.defectPopulation('volumetric', 'Size', 0.05, 'Count', 3), ...
                'flpm:population:location');
            tc.verifyError(@() flpm.defectPopulation('surface', 'Size', 0.01), ...
                'flpm:population:count');
            tc.verifyError(@() flpm.defectPopulation('volumetric', 'Radius', 1, 'Count', 3), ...
                'flpm:population:size');
        end
    end
end

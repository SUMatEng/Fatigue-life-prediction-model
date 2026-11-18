classdef FlpmClusteringTest < matlab.unittest.TestCase
    %FLPMCLUSTERINGTEST  Thomas cluster process: sampling and fitting.

    methods (Test)
        function defaultIsUnchanged(tc)
            % no 'Clustering': the random stream use is that of the
            % Poisson sampler, so earlier results are reproduced
            p = flpm.defectPopulation('volumetric', 'Size', 0.05, 'Depth', 1, 'Count', 5);
            tc.verifyEmpty(p.clustering);
            s = flpm.samplePopulation(p, 'simulation', 1, 'SurfaceLength', 10, 'Seed', 3);
            rs = RandStream('mrg32k3a', 'Seed', 3); rs.Substream = 1;
            tc.verifyEqual(s.position, 10 * rand(rs, 5, 1));
        end

        function clusteredCountAndBounds(tc)
            c = flpm.clustering('MeanClusterSize', 5, 'Spread', [0.05 0.2]);
            p = flpm.defectPopulation('volumetric', 'Size', 0.05, ...
                'Depth', makedist('Exponential', 'mu', 0.5), ...
                'CountPerArea', 2, 'CountModel', 'poisson', 'Clustering', c);
            s = flpm.samplePopulation(p, 300, 'SurfaceLength', 10, 'AxialLength', 3, 'Seed', 1);
            n = cellfun(@(x) numel(x.size), s);
            tc.verifyEqual(mean(n), 60, 'RelTol', 0.1);          % kappa mu |A|
            tc.verifyGreaterThan(var(n), 3 * mean(n));           % over-dispersed
            P = cell2mat(cellfun(@(x) [x.depth x.position x.axial], s, 'UniformOutput', false));
            tc.verifyTrue(all(P(:,1) >= 0 & P(:,2) >= 0 & P(:,2) < 10 & P(:,3) >= 0 & P(:,3) < 3));
            % fixed count: exact total
            p2 = flpm.defectPopulation('volumetric', 'Size', 0.05, 'Depth', 1, ...
                'Count', 12, 'Clustering', c);
            s2 = flpm.samplePopulation(p2, 'simulation', 1, 'SurfaceLength', 10);
            tc.verifyNumElements(s2.size, 12);
            % one cluster holding every defect
            p3 = flpm.defectPopulation('volumetric', 'Size', 0.05, 'Depth', 1, ...
                'Count', 9, 'Clustering', flpm.clustering('MeanClusterSize', 9, 'Spread', 0.1));
            s3 = flpm.samplePopulation(p3, 'simulation', 1, 'SurfaceLength', 10, 'AxialLength', 2);
            tc.verifySize(s3.position, [9 1]);
        end

        function clusteringShortensNeighbourDistance(tc)
            base = {'Size', 0.05, 'Depth', 1, 'Count', 200};
            p0 = flpm.defectPopulation('volumetric', base{:});
            p1 = flpm.defectPopulation('volumetric', base{:}, 'Clustering', ...
                flpm.clustering('MeanClusterSize', 10, 'Spread', 0.05));
            d = zeros(1, 2); P = {p0, p1};
            for k = 1:2
                s = flpm.samplePopulation(P{k}, 'simulation', 1, 'SurfaceLength', 50, ...
                    'AxialLength', 50, 'Seed', 2);
                [~, dd] = knnsearch([s.position s.axial], [s.position s.axial], 'K', 2);
                d(k) = median(dd(:, 2));
            end
            tc.verifyLessThan(d(2), 0.5 * d(1));
        end

        function fitRecoversThomas2D(tc)
            % simulate a Thomas process in a box and recover its parameters
            rs = RandStream('mrg32k3a', 'Seed', 4);
            kappa = 2; mu = 8; sig = 0.08; Lb = 10;
            np = poissinv(rand(rs), kappa * (Lb + 1)^2);
            P = -0.5 + (Lb + 1) * rand(rs, np, 2);
            k = poissinv(rand(rs, np, 1), mu);
            X = repelem(P, k, 1) + sig * randn(rs, sum(k), 2);
            X = X(all(X >= 0 & X <= Lb, 2), :);
            f = flpm.fitClustering(X, 'Box', [0 0; Lb Lb], 'RMax', 1, 'NumSimulations', 99);
            tc.verifyEqual(f.kappa, kappa, 'RelTol', 0.3);
            tc.verifyEqual(f.mu, mu, 'RelTol', 0.3);
            tc.verifyEqual(f.sigma, sig, 'RelTol', 0.2);
            tc.verifyGreaterThan(f.pThomas, 0.05);
            tc.verifyLessThan(f.pPoisson, 0.05);
        end

        function poissonDataGivesWeakClustering(tc)
            rs = RandStream('mrg32k3a', 'Seed', 5);
            X = 10 * rand(rs, 400, 2);
            f = flpm.fitClustering(X, 'Box', [0 0; 10 10], 'RMax', 1.5, 'NumSimulations', 19);
            % K of CSR: the cluster term is negligible over the fitted range
            big = f.K.r >= 0.5;
            tc.verifyLessThan(max(abs(f.K.Khat(big) - pi * f.K.r(big).^2) ./ (pi * f.K.r(big).^2)), 0.15);
            tc.verifyGreaterThan(f.pPoisson, 0.05);
        end

        function fitAnisotropic3D(tc)
            rs = RandStream('mrg32k3a', 'Seed', 6);
            kappa = 1; mu = 10; sp = 0.05; sa = 0.2; R = 3; h = 6;
            lo = [-R-1 -R-1 -1]; hi = [R+1 R+1 h+1];
            np = poissinv(rand(rs), kappa * prod(hi - lo));
            P = lo + (hi - lo) .* rand(rs, np, 3);
            k = poissinv(rand(rs, np, 1), mu);
            X = repelem(P, k, 1) + [sp sp sa] .* randn(rs, sum(k), 3);
            X = X(hypot(X(:,1), X(:,2)) <= R & X(:,3) >= 0 & X(:,3) <= h, :);
            f = flpm.fitClustering(X, 'Cylinder', [R 0 h], 'RMax', 0.8, 'NumSimulations', 9);
            tc.verifyEqual(f.sigma(1), sp, 'RelTol', 0.25);
            tc.verifyEqual(f.sigma(2), sa, 'RelTol', 0.25);
            tc.verifyEqual(f.mu, mu, 'RelTol', 0.35);
            % edge correction: disc window, Poisson -> K ~ pi r^2 near the edge too
            Y = []; while size(Y,1) < 600
                y = 2*R*rand(rs, 800, 2) - R; Y = [Y; y(hypot(y(:,1), y(:,2)) <= R, :)]; %#ok<AGROW>
            end
            g = flpm.fitClustering(Y(1:600, :), 'Disc', R, 'RMax', 1.5, 'NumSimulations', 9);
            tc.verifyEqual(g.K.Khat(end), pi * 1.5^2, 'RelTol', 0.15);
        end

        function inhomogeneousIntensityIsNotClustering(tc)
            % Poisson points concentrated near the free surface of a bar
            rs = RandStream('mrg32k3a', 'Seed', 9);
            R = 3; h = 6; X = zeros(0, 3);
            while size(X, 1) < 500
                y = [2*R*rand(rs, 4000, 2) - R, h*rand(rs, 4000, 1)];
                r = hypot(y(:,1), y(:,2));
                keep = r <= R & rand(rs, 4000, 1) < exp(-(R - r) / 0.4);
                X = [X; y(keep, :)]; %#ok<AGROW>
            end
            X = X(1:500, :);
            fh = flpm.fitClustering(X, 'Cylinder', [R 0 h], 'RMax', 0.8, 'NumSimulations', 39);
            fi = flpm.fitClustering(X, 'Cylinder', [R 0 h], 'RMax', 0.8, 'NumSimulations', 39, ...
                'Intensity', 'radial');
            tc.verifyLessThan(fh.pPoisson, 0.05);          % homogeneous null rejected
            tc.verifyGreaterThan(fi.pPoisson, 0.05);       % inhomogeneous null not rejected
            big = fi.K.r >= 0.4;
            tc.verifyLessThan(max(abs(fi.K.Khat(big) ./ fi.K.Kpoisson(big) - 1)), 0.3);
        end

        function markTests(tc)
            rs = RandStream('mrg32k3a', 'Seed', 7);
            X = 10 * rand(rs, 300, 2);
            s = X(:, 1) + 0.5 * randn(rs, 300, 1);      % size trends with x
            f = flpm.fitClustering(X, 'Box', [0 0; 10 10], 'Sizes', s, 'Depth', X(:, 1), ...
                'NumSimulations', 9, 'NumPermutations', 199);
            tc.verifyLessThan(f.marks.pNN, 0.01);
            tc.verifyGreaterThan(f.marks.depthRho, 0.8);
        end

        function invalidClustering(tc)
            tc.verifyError(@() flpm.clustering('Spread', 0.1), 'flpm:clustering:input');
            tc.verifyError(@() flpm.defectPopulation('volumetric', 'Size', 0.05, ...
                'Depth', 1, 'Count', 3, 'Clustering', struct('a', 1)), ...
                'flpm:population:clustering');
        end
    end
end

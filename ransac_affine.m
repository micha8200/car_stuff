function [bestT, inliers] = ransac_affine(points1, points2, threshold, maxIter)
    % RANSAC to estimate affine transformation between two sets of 2D points
    % points1, points2: Nx2 matrices of corresponding points
    % threshold: distance threshold to consider a point an inlier
    % maxIter: number of RANSAC iterations

    numPoints = size(points1, 1);
    bestInlierCount = 0;
    bestT = eye(3);
    inliers = [];

    for i = 1:maxIter
        % Randomly select 3 point pairs
        idx = randperm(numPoints, 3);
        p1 = points1(idx, :);
        p2 = points2(idx, :);

        % Estimate affine transformation
        T = estimateAffine(p1, p2);

        % Apply transformation to all points
        points1_h = [points1, ones(numPoints, 1)]';
        transformed = T * points1_h;
        transformed = transformed(1:2, :)';

        % Compute distances
        dists = sqrt(sum((transformed - points2).^2, 2));
        currentInliers = find(dists < threshold);
        inlierCount = length(currentInliers);

        % Update best model
        if inlierCount > bestInlierCount
            bestInlierCount = inlierCount;
            bestT = T;
            inliers = currentInliers;
        end
    end
end

function T = estimateAffine(p1, p2)
    % Estimate affine transformation from 3 point pairs
    A = [];
    b = [];
    for i = 1:3
        A = [A;
             p1(i,1), p1(i,2), 1, 0, 0, 0;
             0, 0, 0, p1(i,1), p1(i,2), 1];
        b = [b; p2(i,1); p2(i,2)];
    end
    x = A \ b;
    T = [x(1), x(2), x(3);
         x(4), x(5), x(6);
         0,    0,    1];
end
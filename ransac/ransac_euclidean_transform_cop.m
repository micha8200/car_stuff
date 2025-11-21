function [R_best, t_best, inliers_best] = ransac_euclidean_transform_cop(P1, P2, threshold, max_iter)
% Estimate Euclidean transformation (rotation + translation) using RANSAC
% P1, P2: Nx2 and Mx2 matrices of 2D points (can be different sizes)
% threshold: inlier distance threshold
% max_iter: number of RANSAC iterations

    num_P1 = size(P1, 1);
    num_P2 = size(P2, 1);
    best_inlier_count = 0;
    R_best = eye(2);
    t_best = zeros(2,1);
    inliers_best = [];

    for i = 1:max_iter
        % Randomly sample 2 points from each set
        idx1 = randperm(num_P1, 2);
        idx2 = randperm(num_P2, 2);
        A = P1(idx1, :);
        B = P2(idx2, :);

        % Estimate transformation from A to B
        [R, t] = estimate_euclidean(A, B);

        % Transform all P1 points
        P1_transformed = (R * P1')' + t';

        % For each transformed point, find nearest neighbor in P2
        inliers = [];
        for j = 1:num_P1
            dists = vecnorm(P2 - P1_transformed(j,:), 2, 2);
            [min_dist, idx] = min(dists);
            if min_dist < threshold
                inliers = [inliers; j, idx];
            end
        end

        % Update best model if more inliers found
        if size(inliers,1) > best_inlier_count
            best_inlier_count = size(inliers,1);
            R_best = R;
            t_best = t;
            inliers_best = inliers;
        end
    end
end

function [R, t] = estimate_euclidean(A, B)
% Estimate Euclidean transform (rotation + translation) from A to B
    centroid_A = mean(A);
    centroid_B = mean(B);

    A_centered = A - centroid_A;
    B_centered = B - centroid_B;

    H = A_centered' * B_centered;
    [U, ~, V] = svd(H);
    R = V * U';

    if det(R) < 0
        V(:,2) = -V(:,2);
        R = V * U';
    end

    t = centroid_B' - R * centroid_A';
end
% MATLAB code to estimate 2D Euclidean transformation using RANSAC
% between two point sets: pointsA (2x5) and pointsB (2x8)
% Assumes pointsA is the moving set, pointsB is the fixed set with possible outliers

function [R, t, ixBtoA, ixAtoB] = ransac_euclidean_transform_grok(pointsA, Na, pointsB, Nb, dist_thresh, dist_euclid)
    % pointsA: 2xNa matrix, each column is a point [x; y]
    % pointsB: 2xNb matrix, each column is a point [x; y]
    % dist_thresh: distance threshold for inliers

    R           = eye(2);
    t           = zeros(2, 1);
    ixBtoA      = zeros(1, Nb);
    ixAtoB      = zeros(1, Na);    
    inlier_pa   = zeros(2, min(Nb, Na));
    inlier_pb   = zeros(2, min(Nb, Na));
    if Na==0 || Nb==0
        return
    end
    best_inliers = 0;
    best_R = eye(2);
    best_t = zeros(2, 1);

    is_optim = 1;
    if is_optim
        dA          = sqdist(pointsA(:, 1:Na)');
        dB          = sqdist(pointsB(:, 1:Nb)');
        dABabs      = abs(dB-dA');
        [dAB, ind]  = stalinsort(dABabs(:), 2*dist_thresh); % upon reaching list where pair distances above 2THRESH it is pointless to test them
        % [dAB, ind]=sort(dABabs(:));
        % dAB<dist_thresh % relevant points
        listA       = find(tril(ones(Na), -1));
        listB       = find(tril(ones(Nb), -1));
    end
    
    for iter = 1:length(ind)

        if is_optim
            if dAB(iter)>2*dist_thresh % protection. max(dAB) should not exceed 2*dist_thresh due to stalinsort
                break
            end
            [iB, iA] = ind2sub(size(dABabs), ind(iter));
            % index iA in self-dist group A will point to the relevant pair
            [idxA(1), idxA(2)]=ind2sub([Na Na], listA(iA));
            [idxB(1), idxB(2)]=ind2sub([Nb Nb], listB(iB));
        else % basic: random permutation of 2 points from pool of Na/Nb
            idxA = randperm(Na, 2);
            idxB = randperm(Nb, 2);
        end
        % Sample 2 distinct points from A
        p1 = pointsA(:, idxA(1));
        p2 = pointsA(:, idxA(2));
        
        % Sample 2 distinct points from B
        q1 = pointsB(:, idxB(1));
        q2 = pointsB(:, idxB(2));
        
        % Compute distances
        v_p = p2 - p1;
        v_q = q2 - q1;
        dA  = norm(v_p);
        dB  = norm(v_q);
        
        % Skip if distances don't match (Euclidean preserves distances)
        if abs(dA - dB) > dist_thresh
            continue;
        end
        
        if dA < 1e-6 || dB < 1e-6
            continue; % Avoid division by zero
        end
        
        % Unit vectors
        u_p = v_p / dA;
        u_q = v_q / dB;
        
        % Compute rotation: cos_theta = dot(u_p, u_q), sin_theta = det([u_p, u_q])
        cos_theta = dot(u_p, u_q);
        sin_theta = u_p(1)*u_q(2) - u_p(2)*u_q(1);
        
        R = [cos_theta, -sin_theta; sin_theta, cos_theta];
        
        % Translation
        t = q1 - R * p1;
        
        % Verify second point
        if norm(q2 - (R * p2 + t)) > dist_thresh
            continue;
        end
        
        % Transform all points in A
        pointsA_trans(:, 1:Na) = bsxfun(@plus, R * pointsA(:, 1:Na),  t);
        
        % Count inliers: for each transformed point, check min distance to B
        inliers = 0;
        for i = 1:Na
            diff = pointsA_trans(:, i) - pointsB(:, 1:Nb); % 2xM
            if dist_euclid
                dists = sqrt(sum(diff.^2, 1)); % 1xM
            else
                dists = sum(abs(diff), 1); % manhatttan
            end
            min_d = min(dists);
            if min_d < dist_thresh
                inliers = inliers + 1;
            end
        end
        
        if inliers > best_inliers
            best_inliers = inliers;
            best_R = R;
            best_t = t;
        end
    end
    
    fprintf('ran for %1.0f iterations\n', iter)
    % Refine using all inliers if possible
    if best_inliers >= 2
        % Transform A with best model
        pointsA_trans(:, 1:Na) = bsxfun(@plus, best_R * pointsA(:, 1:Na),  best_t); 
        
        % Find inlier correspondences
        k = int32(0);
        for i = 1:Na
            diff    = pointsA_trans(:, i) - pointsB(:, 1:Nb);
            if dist_euclid
                dists = sqrt(sum(diff.^2, 1)); % 1xM
            else
                dists = sum(abs(diff), 1); % manhatttan
            end
            [min_d, min_idx] = min(dists);
            if min_d < dist_thresh
                ixAtoB(i)       = min_idx; 
                ixBtoA(min_idx) = i; 
                k               = k + 1;
                inlier_pa(:,k)  = pointsA(:, i);
                inlier_pb(:,k)  = pointsB(:, min_idx);
            end
        end

        if k >= 2
            % Compute centroids
            cent_a = mean(inlier_pa(:, 1:k), 2);
            cent_b = mean(inlier_pb(:, 1:k), 2);
            
            % Centered points
            pa_c = bsxfun(@minus, inlier_pa, cent_a);
            pb_c = bsxfun(@minus, inlier_pb, cent_b);
            
            % SVD for rotation
            H = pa_c * pb_c';
            [U, ~, V] = svd(H);
            R_refined = V * U';
            
            % Ensure proper rotation (det=1)
            if det(R_refined) < 0
                V(:, end) = -V(:, end);
                R_refined = V * U';
            end
            
            t_refined = cent_b - R_refined * cent_a;
            
            best_R = R_refined;
            best_t = t_refined;
        end
    end
    
    R = best_R;
    t = best_t;
end

function [valss, ix]=stalinsort(vals, thresh)
valss   = zeros(length(vals), 1);
ix      = zeros(length(vals), 1);
j       = int32(0);
for i=1:length(vals)
    if vals(i)<thresh
        j(1)        = j + 1;
        valss(j)    = vals(i);
        ix(j)       = i;
    end
end
valss = valss(1:j);
ix = ix(1:j);
end

function test
%%
% Example usage:
% Generate sample data
rng(42); % For reproducibility
pointsA = rand(2, 5) * 10; % 5 random points

% Known transformation
theta_true = pi/6; % 30 degrees
R_true = [cos(theta_true), -sin(theta_true); sin(theta_true), cos(theta_true)];
t_true = [2; 3];
t_rep = repmat(t_true, 1, 5);
pointsA_trans = R_true * pointsA + t_rep;

% Add 3 outliers
outliers = rand(2, 3) * 20 + [5; 5];
pointsB = [pointsA_trans, outliers];

% Shuffle B
perm = randperm(8);
pointsB = pointsB(:, perm);

% Run RANSAC
max_iters = 1000;
dist_thresh = 0.5;
[R_est, t_est] = ransac_euclidean_transform(pointsA, pointsB, max_iters, dist_thresh);

% Display results
disp('Estimated Rotation Matrix:');
disp(R_est);
disp('Estimated Translation:');
disp(t_est);
disp('True Rotation Matrix:');
disp(R_true);
disp('True Translation:');
disp(t_true);

end
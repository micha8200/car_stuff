% function script_ransac(seed)

% if nargin<1
%     seed = 2;
% end

% RANSAC: ok with up to 50% outliers
% p - rrequired success probability (99.9%)
% e - outlier rate              (10%, 20%)
% s - minimum number of points to define the function model (line,
% transformation, etc)
% transformations to test:
% translation 
% euclidean (translation + rotation) * preserves length + orientatios 3DOF
% affine (euclidean + non-coherent stretch/squeeze) preserves parallelism
p = 0.99999;
e = 0.5;
s = 3;
nRansacIter = log(1-p)/log(1-(1-e)^s); % nIter>= that number
% nRansacIter = 150;

% [files, products] = matlab.codetools.requiredFilesAndProducts('script_ransac_example_test.m');
% disp(files')
%%
% BAD results!
% 6; 8
% 3; 10

numPoints   = 7; % 14
Nfab        = 9;
degNoiseStd = 7;
r0          = 85; % distance to observers
numObserv   = 2;
tscale      = 10; % object feature displacement scale
oscale      = 30; % observer displacement scale
obsFOVdeg   = [6 8];


YPRE        = [1.5 2 2; -1.5 -0.77 4]*pi/180;

seed        = 1; % 2 8
rng(seed)
points = rand(numPoints, 3) * tscale + rand(1, 3)*r0 ;
obs = rand(numObserv, 3) * oscale ;

set(groot, 'defaultLegendItemHitFcn', @LegendItemHitFnc);


%% plot figures
iseucl          = 0;
isplot          = 1;
isconvert2uv    = 0;
ix3r            = randsample(1:numPoints, numPoints); % shuffle points
for i=1
    % complexity for eachobserver:
    % Na(Na-1) point pairs in A vs Nb(Nb-1) point pairs in B
    % O(N^4) pair-vs-pair full transformation check (norm, dot, cross)
    % option to reduce by:
    % calc point-to-point dist (euclid/manhattan) on group A, and group B
    % Na(Na-1), Nb(Nb-1)
    % calculate difference between point pairs in A and in B total of
    % NaNb(Na-1)(Nb-1) distance differences and sort them (can be stalinsort!) (also O(N^4))
    % procede to apply (norm, dot, cross)  process on a subset of lines (A
    % vs B) that passed the sort (are below a specific threshold - tied to the RANSAC process threshold)



    obj                     = clsObserver(obs(i,:), [], YPRE(i,:), obsFOVdeg(i));
    ae                      = obj.getPOV(points, 0); % object features, as seen by observer
    % LOS                     = obj.getLOS(points);
    ae3                     = obj.getPOV(points(ix3r,:), 1); % object features, in relation to observer (cast 3D data on observer)
    [ae, isns]              = obj.addNoiseFeatures(ae, Nfab, degNoiseStd); % N, std for noise
    [rectWin, isvis]        = obj.getLOS(ae, 'mean');



    tic
    pA          = zeros(2, 50);
    pB          = zeros(2, 50);
    Na          = size(ae3, 1);
    Nb          = sum(isvis);
    pA(:, 1:Na) = ae3';
    pB(:, 1:Nb) = ae(isvis,:)';
    
    if isconvert2uv
        pA = ae2uv_local(pA);
        pB = ae2uv_local(pB);
    end
    fprintf('\n\n----\n')

    [R, t, iBtoA, iAtoB]    = ransac_euclidean_transform_grok(pA,Na, pB, Nb, 0.5, iseucl);
    toc
    ae3r                    = (R*ae3' + t)';
    
    Rangle      = atan2d(R(2, 1), R(1,1));
    R3          = eye(3);
    R3(1:2,1:2) = R;
    ypr         = rad2deg(rotm2eulerZYX(R3));

    fprintf('chosen rotation: %1.0f hits   %1.1f degrees\n', sum(iBtoA>0), Rangle)
    if isplot
        ff      = figure();
        sb(1)   = subplot(2,3,1);hold on;
        sb(2)   = subplot(2,3,4);hold on;
        sb(3)   = subplot(2,3,[2 3 5 6]);hold on;
        clrs    = lines(7);
    
        set(ff, 'CurrentAxes', sb(1))
        title('visible by observer')
        rectangle('Position', rectWin, 'EdgeColor', 'k'); % , 'DisplayName', 'FOV'
        plot(ae(isvis&~isns, 1), ae(isvis&~isns, 2), 'o', 'DisplayName','2D visible', 'MarkerSize',12, 'Color',clrs(1,:));
        plot(ae(~isvis, 1), ae(~isvis, 2), 'x', 'DisplayName','2D missing', 'MarkerSize',12, 'Visible','off', 'Color',clrs(1,:));
        plot(ae(isvis&isns, 1), ae(isvis&isns, 2), 'd', 'DisplayName','2D noise', 'MarkerSize',12, 'Color',clrs(1,:));
        plot(ae3(:, 1), ae3(:, 2), '*', 'DisplayName','3D in2DCoords', 'Color',clrs(4,:));
        % plot(ae3r(:, 1), ae3r(:, 2), '+', 'DisplayName','3D in2DCoords corrected', 'Visible','off', 'Color',clrs(3,:));
        legend('show', 'Location','best')
    
        set(ff, 'CurrentAxes', sb(2))
        title('fixed features for observer + display others')
        rectangle('Position', rectWin, 'EdgeColor', 'k'); % , 'DisplayName', 'FOV'
        plot(ae(isvis&~isns, 1), ae(isvis&~isns, 2), 'o', 'DisplayName','2D visible', 'MarkerSize',12, 'Color',clrs(1,:));
        plot(ae(~isvis, 1), ae(~isvis, 2), 'x', 'DisplayName','2D invisible', 'MarkerSize',12, 'Color',clrs(1,:));
        plot(ae(isns, 1), ae(isns, 2), 'd', 'DisplayName','2D noise', 'MarkerSize',12, 'Color',clrs(1,:));     
        plot(ae3r(:, 1), ae3r(:, 2), '*', 'DisplayName','3D in2DCoords corrected', 'Color',clrs(4,:));
        legend('show', 'Location','best')
    
    
        set(ff, 'CurrentAxes', sb(3))
        N = 3*numPoints;
        pline1 = zeros(N, 3);
        pline1(1:3:N, 1) = points(:,1);
        pline1(1:3:N, 2) = points(:,2);
        pline1(1:3:N, 3) = points(:,3);
        pline1(2:3:N, 1) = obs(i,1);
        pline1(2:3:N, 2) = obs(i,2);
        pline1(2:3:N, 3) = obs(i,3);
        pline1(3:3:N, :) = nan;
    
        scatter3(points(:,1), points(:,2), points(:,3), 100, 'o', 'DisplayName','points');
        scatter3(obs(:,1), obs(:,2), obs(:,3), 100, '*k', 'DisplayName','observers');
        plot3(pline1(:, 1),pline1(:, 2),pline1(:, 3), 'r-', 'DisplayName',sprintf('observer%1.0f\\_LOS', i))
        % axis equal;
        legend('show', 'Location','best')
    end
end



%%


function uv = ae2uv_local(ae)
sinel   = sin(ae(2,:));
cosel   = cos(ae(2,:)); 
sinaz   = sin(ae(1,:));
uv      = [cosel.*sinaz; sinel];
end

% Initial state
SDp = 0.5;

P = eye(2)*SDp; % initial uncertainty
dt = 0.1;   % time step
N = 200;
ts = (0:N-1)*dt;
x0 = 6;
v0 = 8.5;
a0 = 0.1;
x = [7 ; 9];
q0 = 0.1;

% Simulated measurements
measurements = x0 + SDp*randn(1, N) + v0 * ts + 0.5 * a0 * ts.^2;
% ix = floor(N*0.33):N;
% measurements(ix) = measurements(ix) + 0.9;
% ix = floor(N*0.45):N;
% measurements(ix) = measurements(ix) - 1.2;

v_real = v0 + a0 * ts;

xs  = zeros(size(ts));
xe  = zeros(size(ts));
vs  = zeros(size(ts));
ve  = zeros(size(ts));
for k = 1:length(measurements)
    [x, P] = position_velocity_estimator(measurements(k), SDp, x, P, dt, q0);
    % fprintf('Time %.1f: Position = %.2f, Velocity = %.2f\n', k*dt, x(1), x(2));
    xs(k) = x(1);
    vs(k) = x(2);
    xe(k) = sqrt(P(1));
    ve(k) = sqrt(P(4));
end

close all
figure(1)
hold on
plot(ts, measurements,'o', 'DisplayName','meas');
errorbar(ts, xs, xe, 'DisplayName','est');
legend show

figure(2)
hold on
plot(ts, v_real,'o', 'DisplayName','real\_vel');
errorbar(ts, vs, ve, 'DisplayName','est');
legend show

function [x_est, P_est] = position_velocity_estimator(z, R, x_prev, P_prev, dt, q0)
% POSITION_VELOCITY_ESTIMATOR - Estimates position and velocity from position measurements
% using a simple Kalman filter in 1D.
%
% Inputs:
%   z      - Measured position at current time step
%   R      - Measurement noise covariance
%   x_prev - Previous state estimate [position; velocity]
%   P_prev - Previous covariance matrix
%   dt     - Time step duration
%   q0     - Process noise covariance (lower means closer to model, higher means closer to measurements)
%
% Outputs:
%   x_est  - Updated state estimate [position; velocity]
%   P_est  - Updated covariance matrix

% State transition matrix
A = [1 dt; 0 1];

% Measurement matrix (we only measure position)
H = [1 0];

% Process noise covariance
Q = [1 0; 0 1] .* q0; % 0.01

% Predict (propogate)
x_pred = A * x_prev;
P_pred = A * P_prev * A' + Q;

% Kalman gain
K = P_pred * H' / (H * P_pred * H' + R);

% Update
x_est = x_pred + K * (z - H * x_pred);
P_est = (eye(2) - K * H) * P_pred;
end
function ypr = rotationMatrixToEulerAngles(R)
    % Validate input
    if ~isequal(size(R), [3, 3])
        error('Input must be a 3x3 rotation matrix');
    end

    % Check for gimbal lock
    if abs(R(3,1)) < 1 - 1e-6
        pitch = -asin(R(3,1));
        roll  = atan2(R(3,2)/cos(pitch), R(3,3)/cos(pitch));
        yaw   = atan2(R(2,1)/cos(pitch), R(1,1)/cos(pitch));
    else
        % Gimbal lock case
        pitch = pi/2 * sign(-R(3,1));
        roll  = atan2(-R(1,2), R(2,2));
        yaw   = 0;
    end
    ypr = [yaw pitch roll];
end
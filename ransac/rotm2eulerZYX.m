function ypr = rotm2eulerZYX(R)
    % ROTM2EULERZYX Convert rotation matrix to Euler angles (ZYX convention)
    % Input:
    %   R - 3x3 rotation matrix
    % Output:
    %   yaw   - rotation about Z axis
    %   pitch - rotation about Y axis
    %   roll  - rotation about X axis
    
    % Check matrix size
    if ~all(size(R) == [3 3])
        error('Input must be a 3x3 matrix');
    end
    
    % Extract angles
    yaw   = atan2(R(2,1), R(1,1));  % ψ
    pitch = atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2)); % θ
    roll  = atan2(R(3,2), R(3,3));  % φ
    ypr = [yaw pitch roll];
end
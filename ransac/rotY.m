function R = rotY(theta)
% ROTY returns a 3x3 rotation matrix about the Y-axis
% theta: rotation angle in radians

R = [cos(theta)  0  sin(theta);
          0      1      0;
     -sin(theta) 0  cos(theta)];
end
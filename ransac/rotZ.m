function R = rotZ(theta)
% ROTZ returns a 3x3 rotation matrix about the Z-axis
% theta: rotation angle in radians

R = [cos(theta) -sin(theta)  0;
     sin(theta)  cos(theta)  0;
          0           0      1];
end
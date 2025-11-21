function R = rotX(theta)
% ROTX returns a 3x3 rotation matrix about the X-axis
% theta: rotation angle in radians

R = [1      0           0;
     0  cos(theta) -sin(theta);
     0  sin(theta)  cos(theta)];
end
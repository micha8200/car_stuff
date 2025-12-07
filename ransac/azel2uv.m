function [u,v,w] = azel2uv(az, el, Bx, By, Bz)
% AZEL2UV  Transform azimuth/elevation to u-v-w coordinates
%
%   [u,v,w] = azel2uv(az, el, Bx, By, Bz)
%
%   Inputs:
%       az  - Azimuth angle (radians)
%       el  - Elevation angle (radians)
%       Bx, By, Bz - Baseline vector components (in wavelengths)
%
%   Outputs:
%       u,v - Projected baseline coordinates in the u-v plane
%       w   - Projection along source direction

    % Rotation matrix from az/el to uvw
    R = [ sin(az),              cos(az),             0;
         -sin(el)*cos(az),  sin(el)*sin(az),   cos(el);
          cos(el)*cos(az), -cos(el)*sin(az),   sin(el) ];

    % Baseline vector
    B = [Bx; By; Bz];

    % Transform
    uvw = R * B;

    % Extract components
    u = uvw(1);
    v = uvw(2);
    w = uvw(3);
end
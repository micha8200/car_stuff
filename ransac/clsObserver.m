classdef clsObserver < handle
    properties
        name = ''
        pos = [0;0;0]
        OrientError = [0 0 0]  % yaw pitch roll
        Orient      = [0 0 0]  % yaw pitch roll
        R = eye(3)
        RE = eye(3)
        FOV = [55 45]
    end
    methods
        function obj = clsObserver(pos, ypr, yprE, FOVAngle)
            obj.pos(:)          = pos;
            if nargin>1 && ~isempty(ypr)
                obj.Orient(:)       = ypr;
                obj.R               = rotX(ypr(3))*rotY(ypr(2))*rotZ(ypr(1));
            end
            if nargin>2 && ~isempty(yprE)
                obj.OrientError(:)  = yprE;
                obj.RE              = rotX(yprE(3))*rotY(yprE(2))*rotZ(yprE(1));
            end    
            if nargin>3 && ~isempty(FOVAngle)
                obj.FOV(:)          = FOVAngle;
            end
        end

        function [ae, ix] = addNoiseFeatures(obj, ae, N, std)
            % N - number of fabicated points to add
            if nargin<4
                std = [1 1];
            end
            LOS = mean(ae, 1);
            ae1 = LOS + std.*randn(N, 2);
            ix  = [false(size(ae, 1), 1); true(N, 1)];
            ae = [ae;ae1];
            
        end

        function [losRect, ix, ixnon] = getLOS(obj, ae, type)
            % LOS to object features (select if mean, lower-left, etc
            % ae = getPOV(obj, pos, 0); % during calculations of LOS, observer orientation error should be included
            np  = size(ae, 1);
            if nargin<3
                type = 'mean';
            end
            LOS = [0 0];
            switch type
                case 'mean'
                    LOS = mean(ae, 1);
                case 'left'
                    [~, ix]     = min(ae(:, 1), 1);
                    LOS         = ae(ix, :);
                case 'lowerleft'
                    [~, ix]     = min(ae.^2, 1);
                    LOS         = ae(ix, :);
            end
            dAng    = ae - LOS;
            ix      = all(abs(dAng) < 0.5*obj.FOV, 2);
            ixnon   = ~ix;
            losRect = [LOS - 0.5*obj.FOV  obj.FOV(1) obj.FOV(2)];
        end

        function ae = getPOV(obj, posOther, castOnly)
            if nargin<3
                castOnly = true;
            end
            p_rel           = (posOther' - obj.pos)';                    % features position relative to observer
            if ~castOnly
                p_rel           = (obj.RE*obj.R*p_rel')';   % features position relative to observer in his LOS coordinates
            end
            [az, el, rn2]   = cart2sphd(p_rel(:, 1),p_rel(:, 2),p_rel(:, 3)); % in LGC
            ae              = [az, el]; % in deg

        end
        function plotVisible()
        end
    end
end

function [az, el, r]=cart2sphd(x, y, z)
[az, el, r]=cart2sph(x,y,z);
az = -az * 180/pi;
el = el * 180/pi;
end
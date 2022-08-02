classdef WDObject < handle
    % WDObject    class of real world objects
    %
    % The WDObject class is the base of all kinds of real world objects 
    % used in computational cases. By class design any real world object can
    % have 
    %   * a name,
    %   * a unique identifier,
    %   * one or more nested real world objects,
    %   * one or more positions and orientations, and
    %   * one of more time values.
    % Typically, there are as many time values as position and/or orientation
    % values.
    %
    %
    % Properties (private)
    % --------------------
    % R                         position [m] and/or [rad]
    % E__rad                    Euler angles [rad]
    % name                      name or specifier
    % t__s                      time axis [s]
    % uid                       unique identifier or lot number
    % wd_object                 real world objects as components
    %
    %
    % Methods (public)
    % ----------------
    % WDObject                  construct this instance of WDObject class
    % getName                   get name
    % getPosition               get position structure
    % getTime                   get time values
    % getUID                    get unique identifier
    % getVelocity               get velocity structure
    % getWDObject               get nested components
    % setPosition1D             set one or more 1D position values
    % setPosition2Dct           set 2D cartesian position values
    % setPosition3Dct           set 3D cartesian position values
    % setTime                   set time values
    % setVelocity2Dct           set 2D cartesian velocity values
    % computeArcLength          compute arc length of position sequence
    % computeNormOfVelocity     compute norm of velocity vector
    % computeTangentAngle       compute tangential angle along track
    % differentiatePosition     compute position differences
    % differentiateTime         compute time differences
    %
    % ________________________________________________________________________
    %
    % todo:
    % 1. complete get and set methods
    % 2. insert tag of current coordinate system
    % 3. handle unit conversion
    % 4. ...
    %
    % author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    
%% properties (private) %%
    properties (Access = private)
        E__rad                  % Euler angles [rad]
        R                       % position [m] and/or [rad]
        V                       % velocity [m/s] and/or [rad/s]
        name                    % name or specifier
        t__s                    % time axis [s]
        uid                     % unique identifier or lot number
        wd_object               % real world objects as components
    end%properties
    
    
%% methods (public) %%
    methods
    %% constructor %%
        function wdo = WDObject(name, uid)
            % WDObject    construct this instance of WDObject class
            %
            % Construct a new real world object and sets its name and unique
            % identifier.
            % 
            % Arguments:
            % name              name of this real world object
            % uid               unique identifier of this real world object
            % 
            % Returns:
            % wdo               the new real world object
            arguments
                name (1, 1) string = "rwo"
                uid  (1, 1) string = "000"
            end%arguments
            wdo.name = name;
            wdo.uid = uid;
        end
        
    %% get methods %%
        function name = getName(wdo)
            % getName    get name
            name = wdo.name;
        end%function
        function R = getPosition(wdo)
            % getPosition    get position structure
            R = wdo.R;
        end%function
        function t__s = getTime(wdo)
            % getTime    get time values
            t__s = wdo.t__s;
        end%function
        function uid = getUID(wdo)
            % getUID    get unique identifier
            uid = wdo.uid;
        end%function
        function V = getVelocity(wdo)
            % getVelocity    get velocity structure
            V = wdo.V;
        end%function
        function wd_object = getWDObject(wdo)
            % getWDObject    get nested components
            wd_object = wdo.wd_object;
        end%function
        
    %% set methods %%
        function R = setPosition1D(wdo, X__m)
            % setPosition1D    set one or more 1D position values
            % 
            % Set one or more position values of this real world object in 
            % metric units. The method sets the spatial dimension of the 
            % positions in this real world object to one.
            %
            % Arguments:
            % X__m              position values along X axis in [m]
            %
            % Returns:
            % R                 new one dimensional position values [m]   
            arguments
                wdo WDObject
                X__m {mustBeNumeric}
            end%arguments
            R = struct('X__m', X__m);
            wdo.R = R;
        end%function
        function R = setPosition2Dct(wdo, X__m, Y__m)
            % setPosition2Dct    set 2D cartesian position values
            % 
            % Set one or more position values of this real world object in 
            % metric units. The method sets the spatial dimension of the 
            % positions in this real world object to two. The position values
            % are cartesian coordinates, and must have the same size along 
            % each axis.
            %
            % Arguments:
            % X__m              position values along X axis in [m]
            % Y__m              position values along Y axis in [m]
            %
            % Returns:
            % R                 new two dimensional position values [m]   
            arguments
                wdo WDObject
                X__m {mustBeNumeric}
                Y__m {mustBeNumeric}
            end%arguments
            if ~all(size(X__m) == size(Y__m))
                throw(MException( ...
                    'WDObject:SizeError', ...
                    'Values must have equal size along each axis.'));
            else
                R = struct('X__m', X__m, 'Y__m', Y__m);
                wdo.R = R;
            end%if
        end%function
        function R = setPosition3Dct(wdo, X__m, Y__m, Z__m)
            % setPosition3Dct    set 3D cartesian position values
            % 
            % Set one or more position values of this real world object in 
            % metric units. The method sets the spatial dimension of the 
            % positions in this real world object to three. The position 
            % values are cartesian coordinates, and must have the same size 
            % along each axis.
            %
            % Arguments:
            % X__m              position values along X axis in [m]
            % Y__m              position values along Y axis in [m]
            % Z__m              position values along Z axis in [m]
            %
            % Returns:
            % R                 new three dimensional position values [m]   
            arguments
                wdo WDObject
                X__m {mustBeNumeric}
                Y__m {mustBeNumeric}
                Z__m {mustBeNumeric}
            end%arguments
            if ~all((size(X__m) == size(Y__m)) & (size(X__m) == size(Z__m)))
                throw(MException( ...
                    'WDObject:SizeError', ...
                    'Values must have equal size along each axis.'));
            else
                R = struct('X__m', X__m, 'Y__m', Y__m, 'Z__m', Z__m);
                wdo.R = R;
            end%if
        end%function
        function t__s = setTime(wdo, t__s)
            % setTime    set values of time axis
            arguments
                wdo WDObject
                t__s (:, 1) {mustBeNumeric}
            end%arguments
            wdo.t__s = t__s;
        end%function
        function V = setVelocity2Dct(wdo, V_X__m_s, V_Y__m_s)
            % setVelocity2Dct    set 2D cartesian velocity values
            %
            % Set one or more velocity values of this real world object in
            % SI units. The method sets the velocity values in two spatial
            % dimensions of the cartesian system. The velocity values must
            % have the same size along each axis.
            %
            % Arguments:
            % V_X__m_s         velocity values along X axis in [m/s]
            % V_Y__m_s         velocity values along Y axis in [m/s]
            % 
            % Returns:
            % V                new two dimensional velocity values [m/s]
            arguments
                wdo WDObject
                V_X__m_s {mustBeNumeric}
                V_Y__m_s {mustBeNumeric}
            end%arguments
            if ~all(size(V_X__m_s) == size(V_Y__m_s))
                throw(MException( ...
                    'WDObject:SizeError', ...
                    'Values must have equal size along each axis.'));
            else
                V = struct('V_X__m_s', V_X__m_s, 'V_Y__m_s', V_Y__m_s);
                wdo.V = V;
            end%if
        end%function
        
        %% computational methods %%
        function s__m = computeArcLength(wdo)
            % computeArcLength    compute arc length along track
            dR = wdo.differentiatePosition('forward');
            field = fieldnames(dR);
            s__m = dR.(field{1, 1}) .^ 2.0;
            for f = 2:numel(field)
                s__m = s__m + dR.(field{f, 1}) .^ 2.0;
            end%for
            s__m = [0.0 * s__m(1, :); cumsum(sqrt(s__m), 1)];
        end%function
        function v__m_s = computeNormOfVelocity(wdo)
            % computeNormOfVelocity    compute norm of velocity vector
            field = fieldnames(wdo.V);
            v__m_s = wdo.V.(field{1, 1}) .^ 2.0;
            for f = 2:numel(field)
                v__m_s = v__m_s + wdo.V.(field{f, 1}) .^ 2.0;
            end%for
            v__m_s = sqrt(v__m_s);
        end%function
        function Phi__rad = computeTangentAngle(wdo)
            % computeTangentAngle    compute tangential angle along track
            dR = wdo.differentiatePosition('forward');
            Phi__rad = unwrap(atan2(dR.Y__m, dR.X__m), [], 1);
        end%function
        function dR = differentiatePosition(wdo, direction)
            % differentiatePosition    compute position differences
%             arguments
%                 wdo WDObject
%                 direction {mustBeA(direction, "string")} = 'forward';
%             end%arguments
            dR = struct();
            field = fieldnames(wdo.R);
            for f = 1:numel(field)
                dR.(field{f, 1}) = diff(wdo.R.(field{f, 1}), 1, 1);
                if strcmp(direction, 'backward')
                    dR.(field{f, 1}) = [ ...
                        0.0 * dR.(field{f, 1})(1, :); dR.(field{f, 1})];
                elseif strcmp(direction, 'central')
                    dR.(field{f, 1}) = 0.5 * ([ ...
                        0.0 * dR.(field{f, 1})(1, :); dR.(field{f, 1})] + ...
                        [dR.(field{f, 1}); 0.0 * dR.(field{f, 1})(end, :)]);
                else
                    dR.(field{f, 1}) = [ ...
                        dR.(field{f, 1}); 0.0 * dR.(field{f, 1})(end, :)];
                end%if
            end%for
        end%function
        function dt__s = differentiateTime(wdo, direction)
            % differentiateTime    compute time differences
%             arguments
%                 wdo WDObject
%                 direction {mustBeA(direction, "string")} = 'forward';
%             end%arguments
            dt__s = diff(wdo.t__s);
            if strcmp(direction, 'backward')
                dt__s = [0.0; dt__s];
            elseif strcmp(direction, 'central')
                dt__s = 0.5 * ([0.0; dt__s] + [dt__s; 0.0]);
            else
                dt__s = [dt__s; 0.0];
            end%if
        end%function
    end%methods
end%classdef


%% end of module WDObject

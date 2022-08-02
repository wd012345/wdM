classdef WDriveSet < WDObject
    % WDriveSet   set of machine axes drives
    %
    % Set of drives for axis motion on machines, typically realizing
    % cartesion or polar coordinates.
    %
    %
    % Properties (private)
    % --------------------
    % Phi__rad                  tangential angle along track
    % S__m                      metric scales of increments
    % a_n_max__m_s2             max. limit of drives acceleration
    % dR                        nominal and actual increments
    % dt__s                     time increment
    % s__m                      arc length along trajectory
    % v__m_s                    nominal and actual velocity norm
    % v_n_max__m_s              max. limit of drives speed
    %
    %
    % Methods (public)
    % ----------------
    % WDriveSet                 construct this WDriveSet object
    % getAccelerationMax        return nominal acceleration maximum value
    % getArcLength              return arc lengths along track
    % getIncrementOfPosition    return position increments along track
    % getTangentAngle           return tangent angles along track
    % getTime                   return time along track
    % getVelocity               return velocity vectors and norms along track
    % getVelocityMax            return nominal velocity maximum value
    % setAccelerationMax        assign nominal acceleration maximum value
    % setPosition3Dct           scale and assign drives positions
    % setVelocity2Dct           scale and assigne drives velocity
    % setVelocityMax            assign nominal velocity maximum value
    % computeArcLength          compute nominal and actual arc length
    % computeNormOfVelocity     compute norm of velocity vector
    % computeTangentAngle       compute tangential angle along track
    % differentiatePosition     compute position increments
    % differentiateTime         compute time increments
    %    
    %
    % ?setTimeProblemIndex       assign time problems index vector
    % ?synchronizeToFrames       synchronize frame captures to time axis
    %
    % ________________________________________________________________________
    %
    % todo:
    % 1. ...
    %
    % author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    
%% properties (private) %%
    properties (Access = private)
        Phi__rad                % tangential angle along track
        S__m                    % metric scales of increments
        a_n_max__m_s2           % max. limit of drives acceleration
        dR                      % nominal and actual increments
        dt__s                   % time increment
        s__m                    % arc length along trajectory
        v__m_s                  % nominal and actual velocity norm
        v_n_max__m_s            % max. limit of drives speed
%         ?idxTimeProblemVec       % time jumps index
    end
    
    
%% methods %%
    methods
    %% constructor %%
        function wd_drives = WDriveSet(name, uid)
            % WDriveSet  construct this WDriveSet object
            wd_drives = wd_drives@WDObject(name, uid);
            % set SI unit scale factor of position increments and feeds
            wd_drives.S__m.XY = 24.0e-3 / 2^20;
            wd_drives.S__m.Z = 10.0e-3 / 2^20;
            wd_drives.S__m.V = 24.0 * 2^-30;
        end%function
        
    %% get and set methods %%
        function a_n_max__m_s2 = getAccelerationMax(wd_drives)
            % getAccelerationMax    return nominal acceleration maximum value
            a_n_max__m_s2 = wd_drives.a_n_max__m_s2;
        end%function
        function s__m = getArcLength(wd_drives)
            % getArcLength    return arc length along track
            s__m = wd_drives.s__m;
        end%function
        function dR = getIncrementOfPosition(drives)
            % getPositionsIncrements    return position increments along track
            dR = drives.dR;
        end%function
        function Phi__rad = getTangentAngle(wd_drives)
            % getTangentAngles    return tangent angles along track
            Phi__rad = wd_drives.Phi__rad;
        end%function
        function [V, v__m_s] = getVelocity(wd_drives)
            % getVelocities   return velocity vectors and norms along track
            % Return the vector components and their norm value of the
            % feedback, nominal and frame capture velocities.
            V = wd_drives.getVelocity@WDObject();
            v__m_s = wd_drives.v__m_s;
        end%function
        function v_n_max__m_s = getVelocityMax(drives)
            % getVelocityMax   return nominal velocity maximum value
            v_n_max__m_s = drives.v_n_max__m_s;
        end%function
        function setAccelerationMax(drives, a_n_max__m_s2)
            % setAccelerationMax   assign nominal acceleration maximum value
            drives.a_n_max__m_s2 = a_n_max__m_s2;
        end%function
        function setPosition3Dct(drives, X, Y, Z)
            % setPosition3Dct    scale and assign drives positions
            %
            % Set the position in three cartesian dimensions. 
            drives.setPosition3Dct@WDObject( ...
                drives.S__m.XY * X, ...
                drives.S__m.XY * Y, ...
                drives.S__m.Z * Z);
        end%function
        function setVelocity2Dct(wd_drives, V_X, V_Y)
            % setVelocity2Dct    scale and assigne drives velocity
            wd_drives.setVelocity2Dct@WDObject( ...
                wd_drives.S__m.V * V_X, ...
                wd_drives.S__m.V * V_Y);
        end%function
        function setVelocityMax(wd_drives, v_n_max__m_s)
            % setVelocityMax   nominal velocity maximum value
            wd_drives.v_n_max__m_s = v_n_max__m_s;
        end%function

    %% computation methods %%
        function s__m = computeArcLength(wd_drives)
            % computeArcLength    compute nominal and actual arc length
            s__m = wd_drives.computeArcLength@WDObject();
            wd_drives.s__m = s__m;
        end
        function [V, v__m_s] = computeVelocity2D(wd_drives)
            % computeVelocity    compute 2D track velocity and norm
            for s = 1:size(wd_drives.dR.X__m, 2)
                V.V_X__m_s(:, s) = wd_drives.dR.X__m(:, s) ./ wd_drives.dt__s;
                V.V_Y__m_s(:, s) = wd_drives.dR.Y__m(:, s) ./ wd_drives.dt__s;
            end%for
            v__m_s = sqrt(V.V_X__m_s .^ 2.0 + V.V_Y__m_s .^ 2.0);
        end%function
        function v__m_s = computeNormOfVelocity(wd_drives)
            % computeNormOfVelocity    compute norm of velocity vector
            v__m_s = wd_drives.computeNormOfVelocity@WDObject();
            wd_drives.v__m_s = v__m_s;
%             drives.v__m_s = sqrt( ...
%                 gradient(drives.R__m.X, t__s.ebh) .^ 2.0 + ...
%                 gradient(drives.R__m.Y, t__s.ebh) .^ 2.0 + ...
%                 gradient(drives.R__m.Z, t__s.ebh) .^ 2.0);
%             drives.v_n__m_s = sqrt( ...
%                 gradient(drives.R_n__m.X, t__s.ebh) .^ 2.0 + ...
%                 gradient(drives.R_n__m.Y, t__s.ebh) .^ 2.0 + ...
%                 gradient(drives.R_n__m.Z, t__s.ebh) .^ 2.0);
        end%function
        function Phi__rad = computeTangentAngle(wd_drives)
            % computeTangentAngle    compute tangential angle along track
            Phi__rad = wd_drives.computeTangentAngle@WDObject();
            wd_drives.Phi__rad = Phi__rad;
        end%function
        function dR = differentiatePosition(wd_drives, direction)
            % differentiatePosition    compute position increments
            dR = wd_drives.differentiatePosition@WDObject(direction);
            wd_drives.dR = dR;
        end%function
        function dt__s = differentiateTime(wd_drives, direction)
            % differentiateTime    compute time increments
            dt__s = wd_drives.differentiateTime@WDObject(direction);
            wd_drives.dt__s = dt__s;
        end%function
%         function i_tp = findIndexOfTimeStepError(wd_drive)
%                     % find index of time step problem
%             wd_machine.t__s.nok = find(wd_machine.t__s.vbh == 0);
% %             drives.idxTimeProblemVec = cell(length(idxTimeProblem), 1);
% %             idxTimeProblemVec = cell(length(idxTimeProblem), 1);
%             if ~isempty(wd_machine.t__s.nok)
%                 throw(MException("WDMachine:TimeAxisError", ...
%                     "A time axis problem occurred."));
% %             if ~isempty(idxTimeProblem)
% %                 for ijk = 1:length(idxTimeProblem) - 1
% %                     exitflag = 0;
% %                     ijktmp = 1;
% %                     while ~exitflag
% %                         if (idxTimeProblem(ijk + ijktmp) - ...
% %                                 idxTimeProblem(ijk)) > ijktmp
% %                             exitflag = 1;
% %                             drives.idxTimeProblemVec{ijk} = [ ...
% %                                 idxTimeProblem(ijk), ...
% %                                 idxTimeProblem(ijk + ijktmp - 1)];
% %                         else
% %                             ijktmp = ijktmp + 1;
% %                         end%if
% %                     end%while
% %                 end%for
% %                 drives.idxTimeProblemVec{ijk + 1} = [ ...
% %                     idxTimeProblem(end), idxTimeProblem(end)];
%             end%if
%         function df_s = synchronizeToFrames(drives, wd_video)
%             % synchronizeToFrames   synchronize frame captures to time axis
%             % Synchronize the capture of the camera images with the
%             % EtherCAT time axis.
%             
%             % prepare local scope data
%             N_s = length(wd_video);
%             df_s = cell(N_s, 1);
%             % loop over all streams
%             for s = 1:N_s
%                 % compute valid image frames
%                 % (valid: f_v = 1, invalid: f_v = 0)
%                 df_s{s, 1} = zeros(wd_video{s, 1}.getNumberOfFrames(), 1);
% %                 a_f = wd_video{s, 1}.getIndexOfFrameSampling();
%                 f_s = wd_video{s, 1}.getIndexOfFrames();
%                 df_s{s, 1} = [diff(f_s); 0];
% %                 j = 0;
% %                 if f_s(1) == a_f(1)
% %                     for i = 1 : length(f_s)
% %                         if f_s(i) == a_f(i - j)
% %                             df_s{s}(i) = 1;
% %                         else
% %                             j = j + 1;
% %                         end%if
% %                     end%for
% %                 elseif f_s(1) == a_f(1) + 1
% %                     for i = 1 : length(f_s)
% %                         if f_s(i) == a_f(i - j) + 1
% %                             df_s{s}(i) = 1;
% %                         else
% %                             j = j + 1;
% %                         end%if
% %                     end%for
% %                 end%if
%                 % compute compensation of sampling interpolation
%                 dX = zeros(size(f_s));
%                 dY = zeros(size(f_s));
%                 dZ = zeros(size(f_s));
%                 f_red = wd_video{s, 1}.getIndexOfFrameTimes();
%                 t_f__s = wd_video{s, 1}.getTimeOfFrames();
%                 [dt_f__s, dtau_f] = ...
%                     wd_video{s, 1}.getTimeErrorOfFrameSampling();
%                 for f = 1 : numel(t_f__s)
%                     if dt_f__s(f) < 0.0
%                         dX(f) = dtau_f(f) * ( ...
%                             drives.R_n__m.X(f_red(f)) - ...
%                             drives.R_n__m.X(f_red(f) - 1));
%                         dY(f) = dtau_f(f) * ( ...
%                             drives.R_n__m.Y(f_red(f)) - ...
%                             drives.R_n__m.Y(f_red(f) - 1));
%                         dZ(f) = dtau_f(f) * ( ...
%                             drives.R_n__m.Z(f_red(f)) - ...
%                             drives.R_n__m.Z(f_red(f) - 1));
%                     else % 0.0 <= wd_video{s, 1}.dt_f__s(f)
%                         dX(f) = dtau_f(f) * ( ...
%                             drives.R_n__m.X(f_red(f) + 1) - ...
%                             drives.R_n__m.X(f_red(f)));
%                         dY(f) = dtau_f(f) * ( ...
%                             drives.R_n__m.Y(f_red(f) + 1) - ...
%                             drives.R_n__m.Y(f_red(f)));
%                         dZ(f) = dtau_f(f) * ( ...
%                             drives.R_n__m.Z(f_red(f) + 1) - ...
%                             drives.R_n__m.Z(f_red(f)));
%                     end%if
%                 end%for
%                 % compute compensated nominal frame positions
%                 drives.R_nf__m(s, 1).X = ...
%                     drives.R_n__m.X(f_red) - drives.R_n__m.X(f_red(1)) + dX;
%                 drives.R_nf__m(s, 1).Y = ...
%                     drives.R_n__m.Y(f_red) - drives.R_n__m.Y(f_red(1)) + dY;
%                 drives.R_nf__m(s, 1).Z = ...
%                     drives.R_n__m.Z(f_red) - drives.R_n__m.Z(f_red(1)) + dZ;
%                 % compute nominal frame position increments
%                 drives.dR_nf__m(s, 1).X = [diff(drives.R_nf__m(s, 1).X); 0.0];
%                 drives.dR_nf__m(s, 1).Y = [diff(drives.R_nf__m(s, 1).Y); 0.0];
%                 drives.dR_nf__m(s, 1).Z = [diff(drives.R_nf__m(s, 1).Z); 0.0];
%                 % set contour arc length and velocity norm
%                 drives.s_nf__m{1, s} = [0.0; cumsum(sqrt( ...
%                     drives.dR_nf__m(s, 1).X .^ 2.0 + ...
%                     drives.dR_nf__m(s, 1).Y .^ 2.0 + ...
%                     drives.dR_nf__m(s, 1).Z .^ 2.0))];
%                 drives.v_nf__m_s{1, s} = drives.v_n__m_s(f_red);
%                 % compute tangential angles at frame positions
%                 drives.Phi_nf__rad{1, s} = drives.Phi_n__rad(f_red);
%             end%for
%         end%function        

        
        
            
%             try % old code
%                 % nominal (command) position data
%                 drives.R_n__m.X = drives.S__m.XY * double( ...
%                     eth.Drive_X1_AX5112_0000_0203_MDT_Position_command_value);
%                 drives.R_n__m.Y = drives.S__m.XY * double( ...
%                     eth.Drive_Y_Z_AX5206_0000_0203_MDT_1_Position_command_value);
%                 drives.R_n__m.Z = drives.S__m.Z * double( ...
%                     eth.Drive_Y_Z_AX5206_0000_0203_MDT_2_Position_command_value);
%                 % actual (feedback) position data
%                 drives.R__m.X = drives.S__m.XY * (double( ...
%                     eth.Drive_X1_AX5112_0000_0203_AT_Position_feedback_1_value) + ...
%                     double( ...
%                     eth.Drive_X2_AX5112_0000_0203_AT_Position_feedback_1_value)) / 2;
%                 drives.R__m.Y = drives.S__m.XY * double( ...
%                     eth.Drive_Y_Z_AX5206_0000_0203_AT_1_Position_feedback_1_value);
%                 drives.R__m.Z = drives.S__m.Z * double( ...
%                     eth.Drive_Y_Z_AX5206_0000_0203_AT_2_Position_feedback_1_value);
%                 % compute nominal and actual velocities
%                 drives.V_n__m_s.X = drives.S__m.V * double( ...
%                     eth.Drive_X1_AX5112_0000_0203_AT_Effective_velocity_command_value);
%                 drives.V_n__m_s.Y = drives.S__m.V * double( ...
%                     eth.Drive_Y_Z_AX5206_0000_0203_AT_1_Effective_velocity_command_valu);
%                 drives.V__m_s.X = drives.S__m.Z * double( ...
%                     eth.Drive_X1_AX5112_0000_0203_AT_Velocity_feedback_value_1);
%                 drives.V__m_s.Y = drives.S__m.Z * double( ...
%                     eth.Drive_Y_Z_AX5206_0000_0203_AT_1_Velocity_feedback_value_1);
%             catch % then new code
%                 % nominal (command) position data
%                 drives.R_n__m.X = drives.S__m.XY * ...
%                     double(eth.DriveX1_Rx_PositionCommand);
%                 drives.R_n__m.Y = drives.S__m.XY * ...
%                     double(eth.DriveY_Rx_PositionCommand);
%                 drives.R_n__m.Z = drives.S__m.Z * ...
%                     double(eth.DriveZ_Rx_PositionCommand);
%                 % actual (feedback) position data
%                 drives.R__m.X = drives.S__m.XY * ( ...
%                     double(eth.DriveX1_Tx_PositionFeedback1) + ...
%                     double(eth.DriveX2_Tx_PositionFeedback1)) / 2;
%                 drives.R__m.Y = drives.S__m.XY * ...
%                     double(eth.DriveY_Tx_PositionFeedback1);
%                 drives.R__m.Z = drives.S__m.Z * ...
%                     double(eth.DriveZ_Tx_PositionFeedback1);
%                 % compute nominal and actual velocities
%                 drives.V_n__m_s.X = drives.S__m.V * ...
%                     double(eth.DriveX2_Rx_VeloCommand);
%                 drives.V_n__m_s.Y = drives.S__m.V * ...
%                     double(eth.DriveY_Rx_VeloCommand);
%                 drives.V__m_s.X = drives.S__m.V * ...
%                     double(eth.DriveX2_Tx_VeloFeedback1);
%                 drives.V__m_s.Y = drives.S__m.V * ...
%                     double(eth.DriveY_Tx_VeloFeedback1);
%             end
%         end%function
        
    end%methods
    
    methods (Access = private)
        function setPositions(drives, R__m, R_n__m)
            % setPositions   assign non-empty positions vectors
            % Assign given position vectors as far as they are not empty and
            % recompute position increments and dependant properties. If you
            % have to change the values of the drives properties, make sure
            % to call setTimes() first, then call set Positions(), and in a
            % third step call setVelocities().
            if ~isempty(R__m)
                drives.R__m = R__m;
                drives.dR__m = WDFieldSet('Actual Position Increments', struct( ...
                    'X', diff(R__m.X), ...
                    'Y', diff(R__m.Y), ...
                    'Z', diff(R__m.Z)));
                drives.s__m = [0.0; cumsum(sqrt( ...
                    drives.dR__m.X .^ 2.0 + ...
                    drives.dR__m.Y .^ 2.0 + ...
                    drives.dR__m.Z .^ 2.0))];
                drives.Phi__rad = unwrap(atan2(drives.dR__m.Y, drives.dR__m.X));
            end%if
            if ~isempty(R_n__m)
                drives.R_n__m = R_n__m;
                drives.dR_n__m = WDFieldSet('Nominal Position Increments', struct( ...
                    'X', diff(R_n__m.X), ...
                    'Y', diff(R_n__m.Y), ...
                    'Z', diff(R_n__m.Z)));
                drives.s_n__m = [0.0; cumsum(sqrt( ...
                    drives.dR_n__m.X .^ 2.0 + ...
                    drives.dR_n__m.Y .^ 2.0 + ...
                    drives.dR_n__m.Z .^ 2.0))];
                drives.Phi_n__rad = unwrap(atan2(drives.dR_n__m.Y, drives.dR_n__m.X));
            end%if
        end%function
        function setTimes(drives, t__s)
            % setTimes   assign various times sampled along track
            drives.t__s = t__s;
        end%function
        
    end%methods (private)
    
    
end%classdef %DriveSet

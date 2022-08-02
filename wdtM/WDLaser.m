classdef WDLaser < WDObject
    % WDLaser   class of generic laser objects
    %
    % Define properties and methods common to generic laser objects.
    %
    %
    % Properties
    % ----------
    % P_n__W            nominal optical output power
    % P_nf__W           nominal optical output power at frame sampling
    % name              name of laser model
    % uid               uid (e.g. serial id) of this laser object
    %
    %
    % Methods
    % -------
    % WDLaser                   construct generic laser object
    % getPower                  return nominal power values
    % setPowerNominal           assign nominal power values
    % synchronizeToFrames       synchronize property value to frame sampling times
    %
    %
    % Author
    % Stefan Wittwer, info@wittwer-datatools.ch
    %
    % ToDo
    % - ...
    %
    
    properties (Access = private)
        P_n__W          % nominal optical output power
        P_nf__W         % nominal optical output power at frame sampling
    end
    
    methods
    % constructor %
        function wd_laser = WDLaser(name, uid)
            % WDLaser   construct generic laser object
            wd_laser = wd_laser@WDObject(name, uid);
        end
        
    % data I/O methods %
        function clone(laser, eth)
            try % old code
                laser.setPowerNominal( ...
                    eth.A2G_SSC_V2_SolidStateCuttinghead_RxPDO_Laserpower);
            catch % then new code
                laser.setPowerNominal(eth.Laser_Tx_InputPower);
            end%try
        end%function

    % get methods %
        function [P_n__W, P_nf__W] = getPower(wd_laser)
            % getPower   return nominal power values
            P_n__W = wd_laser.P_n__W;
            P_nf__W = wd_laser.P_nf__W;
        end
        
    % set methods %
        function setPowerNominal(wd_laser, P_n__W)
            % setPowerNominal   assign nominal power values
            wd_laser.P_n__W = max(0.0, double(P_n__W));
        end%function
        
        function synchronizeToFrames(wd_laser, wd_video)
            % synchronizeToFrames   synchronize to frame sampling times
            N_s = length(wd_video);
            % loop over all streams
            for s = 1:N_s
                f_red = wd_video{s, 1}.getIndexOfFrames();
                wd_laser.P_nf__W{1, s} = wd_laser.P_n__W(f_red);
            end%for
        end%function
    end
end


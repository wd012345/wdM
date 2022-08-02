classdef WDMachine < WDObject
    % WDMachine   class of generic machine object
    %
    % Define the class WDMachine which is the base class of all objects 
    % representing generic processing machines. The generic machine objects
    % contain subcomponents, typically drives, machining heads etc.
    %
    %
    % Properties (private)
    % --------------------
    % wd_drives                 drive set containing drive axes
    % wd_head                   generic processing head object
    % wd_laser                  generic laser object
    %
    %
    % Methods
    % -------
    % WDMachine         construct this generic machine object
    % getDriveSet       get drive set component
    % getHead           get machining head component
    % getLaser          get generic machining laser
    % getName           get name and uid of this machine
    % getNozzle         get nozzle object of the head on this machine
    %
    % ________________________________________________________________________
    % 
    % todo:
    % 1. constructor name and uid must be cell arrays to contain components
    %    names and uids as well
    % 2. ...
    %
    % author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    properties (Access = private)
        wd_drives      WDriveSet
        wd_head        WDHead
        wd_laser       WDLaser
    end%properties
    
    methods
    % constructor %
        function wd_machine = WDMachine(name, uid)
            % WDMachine   construct this generic machine object
            %
            % Construct instance of this generic machine object and its
            % subcomponents. Initialize object data container properties.
            %
            % Arguments
            % * name        name of this generic machine object
            % * uid         uniqe identifier of this generic machine object
            %
            % Returns
            % wd_machine    this generic machine object
            wd_machine = wd_machine@WDObject(name, uid);
            % construct machine components
            wd_machine.wd_drives = WDriveSet('Drives', 'DrivesUID');
            wd_machine.wd_head = WDHead('Head', 'HeadUID');
            wd_machine.wd_laser = WDLaser('Laser', 'LaserUID');
        end%function
        
    % data I/O methods %
        function clone(wd_machine, wd_io, data)
            % clone  assign input data to machine and subcomponents
            %
            % Arguments
            % * wd_io           data I/O objects
            % * data            input data of machine and subcomponents
            % * parameter       case parameter
            
            % assign sampling times of Beckhoff EtherCAT and Vision
            wd_machine.t__s.vtc = 1.0e-6 * ( ...
                double(data.vis.tsCam) - double(data.vis.tsCam(1)));
            wd_machine.t__s.vux = 1.0e-6 * ( ...
                double(data.vis.unixtime) - double(data.vis.unixtime(1)));
            wd_machine.t__s.vbh = 1.0e-9 * (double(data.vis.beckhofftime));
            wd_machine.t__s.ebh = 1.0e-9 * (double(data.eth.beckhofftime));
            % assign input data to subcomponents
            wd_machine.drives.clone(data.cfg, data.eth);
            wd_machine.head.clone(wd_io, data);
            wd_machine.laser.clone(data.eth);
        end%function
        
    % pre- and postprocessing methods %
    function [t__s, dt__s, df_s] = preprocess(wd_machine, c_0)
            % preprocess  preprocess generic machine data
            
            % compute forward time increments
            wd_machine.dt__s.ebh = [diff(wd_machine.t__s.ebh); 0.0];
            wd_machine.dt__s.vtc = [diff(wd_machine.t__s.vtc); 0.0];
            wd_machine.dt__s.vbh = [diff(wd_machine.t__s.vbh); 0.0];
            wd_machine.dt__s.vux = [diff(wd_machine.t__s.vux); 0.0];
            % preprocess drive set
            wd_machine.getDriveSet().preprocess(wd_machine.t__s);
            % preprocess camera
            wd_machine.getCamera().preprocess(wd_machine.t__s, c_0);
            % synchronize drives increments to camera frames
            df_s = wd_machine.getDriveSet().synchronizeToFrames( ...
                wd_machine.getCamera().getVideo());
            % return time structures
            t__s = wd_machine.t__s;
            dt__s = wd_machine.dt__s;
        end%function
        
    %% get methods %%        
        function wd_drives = getDriveSet(wd_machine)
            % getDriveSet   get drive set object of this machine
            wd_drives = wd_machine.wd_drives;
        end%function
        function wd_head = getHead(wd_machine)
            % getHead    get machining head component
            wd_head = wd_machine.wd_head;
        end%function
        function wd_laser = getLaser(wd_machine)
            % getLaser    get generic machining laser
            wd_laser = wd_machine.wd_laser;
        end%function        
        function t__s = getTime(wd_machine)
            % getTime   return various times sampled on trajectory
            t__s = wd_machine.t__s;
        end%function
        function dt__s = getTimeIncrement(wd_machine)
            % getTimeIncrement  return various time increments on trajectory
            dt__s = wd_machine.dt__s;
        end%function
        
%         function z_0__m = getActualFocusPosition(machine)
%             z_0__m = machine.head.getActualFocusPosition();
%         end%function
%         
%         function [Phi__rad, Phi_n__rad, Phi_nf__rad] = ...
%                 getAnglesOfTrajectory(machine)
%             [Phi__rad, Phi_n__rad, Phi_nf__rad] = machine.drives.getAngles();
%         end%function
%         
%         function [s__m, s_n__m, s_nf__m] = getArcLength(machine)
%             [s__m, s_n__m, s_nf__m] = machine.drives.getArcLength();
%         end%function
%         
%         function d_n__m = getDiameterOfNozzle(machine)
%             d_n__m = machine.head.getDiameterOfNozzle();
%         end%function
%         
%         function [R__m, R_f__m, R_n__m, R_nf__m] = getPositions(machine)
%             [R__m, R_f__m, R_n__m, R_nf__m] = machine.drives.getPositions();
%         end%function
%         
%         function [dR__m, dR_f__m, dR_n__m, dR_nf__m] = ...
%                 getPositionsIncrements(machine)
%             [dR__m, dR_f__m, dR_n__m, dR_nf__m] = ...
%                 machine.drives.getPositionsIncrements();
%         end%function
%         
%         function [ ...
%                 V__m_s, V_n__m_s, V_nf__m_s, ...
%                 v__m_s, v_n__m_s, v_nf__m_s] = getVelocities(machine)
%             [V__m_s, V_n__m_s, V_nf__m_s, ...
%                 v__m_s, v_n__m_s, v_nf__m_s] = machine.drives.getVelocities();
%         end%function
%         
%         function fig_camera = postprocessHead(machine, title)
%             fig_camera = machine.head.postprocessHead(title);
%         end
%         
%         function fig_times = postprocessMachine(machine, title)
%             % VIEW_TIME_SCALES Create Graphs of Time Scales
%             
%             % create array of graphics objects
%             fig_times = gobjects(2, 1);
%             axs = gobjects(2, 1);
%             lgd = gobjects(2, 1);
%             
%             % create figure with time scales
%             fig_times(1) = figure(11);
%             fig_times(1).Color = [1 1 1];
%             fig_times(1).Visible = 'off';
%             axs(1) = axes('Parent', fig_times(1));
%             plot(axs(1), machine.drives.t__s, '.');
%             axs(1).NextPlot = 'add';
%             plot(axs(1), ...
%                 machine.drives.t_camb__s - machine.drives.t_camb__s(1), ...
%                 ':*', 'MarkerSize', 1.5);
%             plot(axs(1), machine.drives.t_f__s, ':o', 'MarkerSize', 1.5);
%             axs(1).NextPlot = 'replace';
%             axs(1).Box = 'on';
%             axs(1).FontName = 'Verdana';
%             axs(1).Title.String = {'Time Scales'; title};
%             axs(1).Title.FontName = 'Verdana';
%             axs(1).Title.Interpreter = 'none';
%             axs(1).XGrid = 'on';
%             axs(1).XMinorGrid = 'on';
%             axs(1).XLabel.String = 'Sample No.';
%             axs(1).YGrid = 'on';
%             axs(1).YMinorGrid = 'on';
%             axs(1).YLabel.String = 'Time Scales [s]';
%             lgd(1) = legend(axs(1), ...
%                 'EtherCAT', 'Frame Beckhoff', 'Illuminated');
%             lgd(1).FontName = 'Verdana';
% 
%             % create figure with time increments
%             fig_times(2) = figure(12);
%             fig_times(2).Color = [1 1 1];
%             fig_times(2).Visible = 'off';
%             axs(2) = axes('Parent', fig_times(2));
%             plot(axs(2), diff(machine.drives.t__s), '.');
%             axs(2).NextPlot = 'add';
%             plot(axs(2), diff(machine.drives.t_camb__s), ...
%                 ':*', 'MarkerSize', 1.5);
%             plot(axs(2), diff(machine.drives.t_f__s), ':o', 'MarkerSize', 1.5);
%             axs(2).NextPlot = 'replace';
%             axs(2).Box = 'on';
%             axs(2).FontName = 'Verdana';
%             axs(2).Title.String = {'Time Increments'; title};
%             axs(2).Title.FontName = 'Verdana';
%             axs(2).Title.Interpreter = 'none';
%             axs(2).XGrid = 'on';
%             axs(2).XMinorGrid = 'on';
%             axs(2).XLabel.String = 'Sample No.';
%             axs(2).YGrid = 'on';
%             axs(2).YMinorGrid = 'on';
%             axs(2).YLabel.String = 'Time Increments [s]';
%             axs(2).YLim = [-1, 1] * max(diff(machine.drives.t_f__s));
%             lgd(2) = legend(axs(2), ...
%                 'EtherCAT', 'Frame Beckhoff', 'Illuminated');
%             lgd(2).FontName = 'Verdana';
%             
%             fig_times(1).Visible = 'on';
%             fig_times(2).Visible = 'on';
%         end%function
% 
%         
%         
%         
%         %%%%%%%
%         function [d_n, f_p__1_m, p_0, q_0, p_c, q_c, M_c, p_n, q_n] = ...
%                 getDiameterOfAperture(machine)
%             % PREPROCESS_DATA Prepare Input Data for Processing
%             %   Synchronize frame capture times with EtherCAT time. Read
%             %   the image files extracted from the video. Identify the
%             %   stationary form from pixel value distribution and determine
%             %   its diameter in pixel units. Compute the pixel scale from
%             %   known aperture diameter and mask the outside region.
%             
%             [d_n, f_p__1_m, p_0, q_0, p_c, q_c, M_c, p_n, q_n] = ...
%                 machine.head.camera.getDiameterOfAperture( ...
%                 machine.getDiameterOfNozzle());
%         end%function
    end%methods
end%classdef


%     %initialize machine components for app computation
%     %rev.20171031/swi
%     function initAppComputation(machine)
%       machine.head.optic.initAppComputation();
%       machine.head.initAppComputation();
%     end
%     
%     %initialize machine components for app simulation
%     %rev.20171024/swi
%     function initAppSimulation(machine, D__m)
%       machine.head.optic.initAppSimulation(D__m);
%       machine.head.initAppSimulation();
%     end
%     
%     %get point data of machine objects (PCC solver output)
%     %rev.20171010/swi
%     function [pd_nozzle, pd_optic, pd_beam_tool] = getPointDataPCC(machine)
%       [pd_nozzle, pd_optic, pd_beam_tool] = machine.head.getPointDataPCC();
%     end
%     
%     %get point data of machine objects (PDC solver output)
%     %rev.20171010/swi
%     function [pdc_head, pdc_nozzle, pdc_optic, pdc_beam_tool] = ...
%         getPointDataPDC(machine)
%       [pdc_head, pdc_nozzle, pdc_optic, pdc_beam_tool] = ...
%         machine.head.getPointDataPDC();
%     end
%     
%     %set point data of machine components (PCC solver input)
%     %rev.20171004/swi
%     function setPointDataPCC(machine, point_data)
%       machine.head.setPointDataPCC( ...
%         point_data(:, [1 3 4]), machine.v_c_max__m_s);
%       machine.laser.setPointDataPCC(point_data{1, 2});
%     end
%     
%     %set point data of machine components (PDC solver input)
%     %rev.20171010/swi
%     function setPointDataPDC(machine, point_data)
%       machine.head.setPointDataPDC(point_data(:, [2 3]));
%     end
%         
%     %set system data of machine objects
%     %rev.20171004/swi
%     function setSystemData(machine, system_data)
%       %set system data to machine components
%       machine.head.setSystemData( ...
%         system_data{2}, ...  SDHead
%         system_data{1}, ...  SDGas
%         system_data{5}, ...  SDNozzle
%         system_data{6});   % SDOptic
%       machine.laser.setSystemData(system_data{3});
%       %set system data of machine object
%       machine.a_c_max__m_s2 = system_data{4}{2, 1};
%       machine.v_c_max__m_s = system_data{4}{2, 2};
%       machine.uid = system_data{4}{2, 3};
%     end
%     
%     % view systems data table of machine
%     % rev.20170410/swi
%     function sd_machine = viewSystemData(machine)
%       sd_machine = table();
%       sd_machine.a_c_max__m_s2 = machine.a_c_max__m_s2;
%       sd_machine.v_c_max__m_s = machine.v_c_max__m_s;
%       sd_machine.uid = machine.uid;
%     end
% 
%     % vizualize unit part data of machine components
%     % rev.20170522/swi
%     function vizualizePartDataArray(machine, n_s, u_c)
%       % create figure
%       fig01PCC = figure();
%       axs010PCC = axes();
%       axs010PCC.Parent = fig01PCC;
%       axs010PCC.NextPlot = 'add';
%       lin010PCC = line(n_s, machine.laser.c_m);
%       lin011PCC = line(n_s, machine.laser.q_L);
%       lin012PCC = line(n_s, machine.laser.q_m);
%       lin013PCC = line(n_s, machine.head.q_c);
%       lin014PCC = line(n_s, machine.head.q_H);
%       lin015PCC = line(n_s, u_c);
%       axs010PCC.NextPlot = 'replace';
%       lin010PCC.Color = [0.00 0.00 1.00];
%       lin011PCC.Color = [0.00 0.20 0.80];
%       lin012PCC.Color = [0.00 0.20 0.60];
%       lin013PCC.Color = [0.00 1.00 0.00];
%       lin014PCC.Color = [0.00 0.65 0.00];
%       lin015PCC.Color = [0.75 0.75 0.75];
%       lin010PCC.LineWidth = 1.25;
%       lin011PCC.LineWidth = 2.00;
%       lin012PCC.LineWidth = 1.25;
%       lin013PCC.LineWidth = 2.00;
%       lin014PCC.LineWidth = 2.00;
%       lin015PCC.LineWidth = 1.25;
%       lin010PCC.Parent = axs010PCC;
%       lin011PCC.Parent = axs010PCC;
%       lin012PCC.Parent = axs010PCC;
%       lin013PCC.Parent = axs010PCC;
%       lin014PCC.Parent = axs010PCC;
%       lin015PCC.Parent = axs010PCC;
%       axs010PCC.XGrid = 'on';
%       axs010PCC.XLabel.String = 'Node Number n_s';
%       axs010PCC.XMinorGrid = 'on';
%       axs010PCC.YGrid = 'on';
%       axs010PCC.YLabel.String = 'Unit Part Data';
%       axs010PCC.YMinorGrid = 'on';
%       %axs010PCC.YScale = 'log';
%       legend(axs010PCC, ...
%         'laser.c_m', 'laser.q_L', 'laser.q_m', ...
%         'head.q_c', 'head.q_H', ...
%         'u_c', ...
%         'Location', 'northeast');
%     end
%     
% 
% end of Machine.m

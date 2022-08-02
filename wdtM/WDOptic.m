classdef WDOptic < WDObject
    % WDOptic   class of generic optic objects
    %
    % Define the class of generic optic objects (e.g. lenses, waveguides).
    %
    %
    % Properties (private)
    % --------------------
    % M_w               magnification of focal radius
    % dz_0__m           shift of focus position (absolute value)
    % dz0_zRdP__1_W     shift of focus position (relative to focus length)
    % z_0__m            actual focus position
    % z_0n__m           nominal focus position
    %
    %
    % Methods (public)
    % ----------------
    % getFocusPosition         get actual and nominal focus position
    % setFocusPosition         set nominal and actual focus position
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
%         M_w             % magnification of focal radius
%         dz_0__m         % shift of focus position (absolute value)
%         dz0_zRdP__1_W   % shift of focus position (relative to focus length)
        z_0__m          % actual focus position
        z_0n__m         % nominal focus position
    end
    
    
%% methods (public) %%
    methods (Access = public)
    %% constructor %%
        function wd_optic = WDOptic(name, uid)
            % WDOptic    construct this optic object
            wd_optic = wd_optic@WDObject(name, uid);
%             optic.beam_tool = SuperGauss();
        end%function
        
    %% get methods %%
        function [z_0__m, z_0n__m] = getFocusPosition(optic)
            % getFocusPosition    get actual and nominal focus position
            z_0__m = optic.z_0__m;
            z_0n__m = optic.z_0n__m;
        end%function
        
    %% set methods %%
        function setFocusPosition(optic, z_0__m, z_0n__m)
            % setFocusPositions   set nominal and actual focus positions
            if ~isempty(z_0n__m)
                optic.z_0n__m = z_0n__m;
            end%if
            if ~isempty(z_0__m)
                optic.z_0__m = z_0__m;
            end%if
        end%function
        
    %% other methods ...
    
%     %compute power induced focus shift and effective focus position
%     %           1      dP_L
%     %  dz_0 = ----- * ------ * P_L * z_R   ;   z_0_eff = z_0_nom + dz_0
%     %          z_R     dz_0
%     %rev.20170913/swi
%     function [z_0_eff__m, dz_0__m] = computeFocusPositionEff(optic, P_L__W)
%       %compute focus shift
%       dz_0__m = optic.dz0_zRdP__1_W * P_L__W .* optic.beam_tool.z_R__m;
%       %compute effective focus position
%       z_0_eff__m = optic.z_0_nom__m + dz_0__m;
%       %read results to optic object
%       optic.dz_0__m = dz_0__m;
%       optic.z_0_eff__m = z_0_eff__m;
%     end
%     
%     % compute nominal focus position by considering power induced focus shift
%     %                         1      dP_L
%     %   z_0,nom = z_0,eff - ----- * ------ * z_R
%     %                        z_R     dz_0
%     % rev.20170720/swi
%     function z_0_nom__m = computeFocusPositionNom(optic, P_L__W)
%       dz_0__m = optic.dz0_zRdP__1_W * P_L__W .* optic.beam_tool.z_R__m;
%       z_0_nom__m = optic.z_0_eff__m - dz_0__m;
%       optic.z_0_nom__m = z_0_nom__m;
%     end
%     
%     % compute focus or Rayleigh length in [m] from given magnification according
%     % rev.20161110/swi
%     function z_R__m = computeFocusLength(optic, d_c__m)
%       z_R__m = optic.beam_tool.computeFocusLength(optic.M_w, d_c__m);
%     end
%     
%     %initialize optic components/attributes for app computation
%     %rev.20171031/swi
%     function initAppComputation(optic)
%       %compute focus (Rayleigh) length
%       optic.beam_tool.computeFocusLength(optic.M_w);
% 
%     end
%     
%     %initialize optic components/attributes for app simulation
%     %rev.20171024/swi
%     function initAppSimulation(optic, D__m)
%       %compute focus (Rayleigh) length
%       optic.beam_tool.computeFocusLength(optic.M_w);
%       optic.beam_tool.computeWaistRadius();
%       %compute effective focus position
%       P_L_avg__W = optic.beam_tool.c_m .* optic.beam_tool.P_L__W;
%       optic.computeFocusPositionEff(P_L_avg__W);
%       % compute top and ejection beam radii
%       optic.beam_tool.computeRadiusAtTop(optic.z_0_eff__m);
%       optic.beam_tool.computeRadiusAtEjection(D__m, optic.z_0_eff__m);
%     end
%     
%     %get point data of optic object (PCC solver output)
%     %rev.20171004/swi
%     function [pd_optic, pd_beam_tool] = getPointDataPCC(optic)
%       %get point data of optic components
%       pd_beam_tool = optic.beam_tool.getPointDataPCC();
%       %get point data of optic object
%       pd_optic = [ ...
%         (0 : length(optic.M_w) - 1)', ...
%         optic.M_w, ...
%         1.0e3 * optic.dz_0__m, ...
%         1.0e3 * optic.z_0_eff__m];
%     end
% 
%     %get point data of optic object (PDC solver output)
%     %rev.20171010/swi
%     function [pdc_optic, pdc_beam_tool] = getPointDataPDC(optic)
%       %get point data of optic components
%       pdc_beam_tool = optic.beam_tool.getPointDataPDC();
%       %get point data of optic object
%       pdc_optic = [ ...
%         (0 : length(optic.M_w) - 1)', ...
%         optic.M_w, ...
%         1.0e3 * optic.z_0_nom__m];
%     end
% 
%     %set point data of optic object (PCC solver input)
%     %rev.20171004/swi
%     function setPointDataPCC(optic, pd_optic)
%       optic.M_w = pd_optic(:, 2);
%       optic.z_0_nom__m = 1.0e-3 * pd_optic(:, 3);
%     end
%     
%     %set system data of optic object
%     %rev.20171004/swi
%     function setSystemData(optic, sd_optic)
%       optic.dz0_zRdP__1_W = 1.0e-3 * sd_optic{2, 1};
%     end
%     
%     % view system data table of optic
%     % rev.20170413/swi
%     function sd_optic = viewSystemData(optic)
%       sd_optic = table();
%       sd_optic.dz0_zRdP__1_W = optic.dz0_zRdP__1_W;
%     end
%     
%     % write process characteristic
%     % rev.20170531/swi
%     function process_characteristic = ...
%         writeProcessCharacteristic(optic, process_characteristic, id)
%       % create matrix of zeros for process characteristic of workpiece
%       pcc_optic = zeros(length(id), 2);
%       % insert index colum
%       pcc_optic(:, 1) = id;
%       % insert process characteristics
%       pcc_optic(:, 2) = 1.0e3 * optic.z_0_eff__m;
%       % return results
%       process_characteristic{5} = pcc_optic;
%     end
    end%methods
end%classdef


% end of Optic.m

classdef WDHead < WDObject
    % WDHead   class of generic processing head
    %
    % The class WDHead contains properties and methods common to all
    % processing machine head components. The machining head contains its 
    % own components listed in the properties.
    %
    %
    % Properties (private)
    % --------------------
    % wd_camera                 generic camera component
    % wd_nozzle                 generic nozzle component
    % wd_optic                  generic optic component
    %
    %
    % Methods
    % -------
    % WDHead                    construct generic processing head object
    % getCamera                 get generic camera component
    % getNozzle                 get generic nozzle component
    % getOptic                  get generic optic component
    % preprocessData    preprocess head subcomponents
    %
    %
    % ToDos
    % - review commented code
    % - ...
    %
    % Authors
    % Stefan Wittwer, info@wittwer-datatools.ch    
    
    
%% properties (private) %%
    properties (Access = private)
        wd_camera  WDCamera
        wd_nozzle  WDNozzle
        wd_optic   WDOptic
    end%properties
    
    methods
    %% constructor %%
        function wd_head = WDHead(name, uid)
            % WDHead    construct generic processing head object
            %
            % Construct this generic head object and call subcomponent
            % constructors.
            wd_head = wd_head@WDObject(name, uid);
            wd_head.wd_camera = WDCamera('Camera', 'CameraUID');
            wd_head.wd_nozzle = WDNozzle('Nozzle', 'NozzleUID');
            wd_head.wd_optic = WDOptic('Optic', 'OpticUID');
        end%function
                
    %% get methods %%
        function wd_camera = getCamera(wd_head)
            % getCamera    get generic camera component
            wd_camera = wd_head.wd_camera;
        end
        function wd_nozzle = getNozzle(wd_head)
            % getNozzle    get generic nozzle component
            wd_nozzle = wd_head.wd_nozzle;
        end%function
        function wd_optic = getOptic(wd_head)
            % getOptic    get generic optic component
            wd_optic = wd_head.wd_optic;
        end%function
        
        
        
        %% other methods ...
        function frame_data = preprocessData(head, ...
                wd_preprocessor, fb_parameter, t__s, t_camb__s)
            % preprocessData   preprocess head subcomponents
            % Call preprocessing methods of each subcomponent of this
            % generic processing head object, thereby passing the
            % preprocessor instance.
            frame_data = head.camera.preprocessData( ...
                wd_preprocessor, fb_parameter, t__s, t_camb__s);
            head.nozzle.preprocessData(wd_preprocessor);
            head.optic.preprocessData(wd_preprocessor);
        end%function
        
%         function fig_camera = postprocessHead(head, title)
%             fig_camera = head.camera.postprocessCamera(title);
%         end%function
%                
    end%methods
end%classdef

    
%     %initialize head components for app simulation
%     %rev.20171024/swi
%     function initAppSimulation(head)
%       W_t__m = head.optic.beam_tool.W_t__m;
%       head.nozzle.computeCouplingEfficiency(W_t__m);
%     end
%     
%     %get point data of head objects (PCC solver output)
%     %rev.20171004/swi
%     function [pd_nozzle, pd_optic, pd_beam_tool] = getPointDataPCC(head)
%       pd_nozzle = head.nozzle.getPointDataPCC();
%       [pd_optic, pd_beam_tool] = head.optic.getPointDataPCC();
%     end
%     
%     %get point data of head objects (PDC solver output)
%     %rev.20171010/swi
%     function [pdc_head, pdc_nozzle, pdc_optic, pdc_beam_tool] = ...
%         getPointDataPDC(head)
%       %get point data of head components
%       pdc_nozzle = head.nozzle.getPointDataPDC();
%       [pdc_optic, pdc_beam_tool] = head.optic.getPointDataPDC();
%       %get point data of head object
%       pdc_head = [ ...
%         (0 : length(head.u_c) - 1)', ...
%         1.0e-5 * head.dp_H__Pa, ...
%         60.0 * head.v_c__m_s];
%     end
%     
%     %set point data of head object (PCC solver input)
%     %20171004/swi
%     function setPointDataPCC(head, point_data, v_c_max__m_s)
%       %set point data of head components
%       head.nozzle.setPointDataPCC(point_data{1, 2});
%       head.optic.setPointDataPCC(point_data{1, 3});
%       %set point data of head object
%       pd_head = point_data{1, 1};
%       head.dp_H__Pa = 1.0e5 * pd_head(:, 2);
%       head.u_c = pd_head(:, 3);
%       head.v_c__m_s = max(head.v_c_min__m_s, pd_head(:, 4) / 60.0);
%       %set relative point data
%       head.q_H = head.dp_H__Pa / head.dp_H_max__Pa;
%       head.q_c = head.v_c__m_s / v_c_max__m_s;
%     end
%     
%     %set point data of head objects (PDC solver input)
%     %rev.20171010/swi
%     function setPointDataPDC(head, point_data)
%       %set point data of head components
%       head.nozzle.setPointDataPDC(point_data{1, 2});
%       %set point data of head object
%       pc_head = point_data{1, 1};
%       head.u_c = pc_head(:, 2);
%     end
%     
%     %set system data of head objects
%     %rev.20171004/swi
%     function setSystemData(head, sd_head, sd_gas, sd_nozzle, sd_optic)
%       %SDHead: set system data of head object
%       head.dp_H_max__Pa = 1.0e5 * sd_head{2, 1};
%       %SDGas: set type of processing gas
%       switch sd_gas{2, 1}
%         case 'N2'
%           head.gas = N2();
%         case 'O2'
%           head.gas = O2();
%         otherwise
%           error('Unknown type of processing gas.');
%       end
%       %SDNozzle: set type of nozzle
%       switch sd_nozzle{2, 1}
%         case 'HK'
%           head.nozzle = HK();
%         case 'NK'
%           head.nozzle = NK();
%         case 'PC'
%           head.nozzle = PC();
%         otherwise
%           error('Unkown type of nozzle.');
%       end
%       %set system data of optic object
%       head.optic.setSystemData(sd_optic);
%     end
%     
%     % write part data from head object
%     % rev. 20170719/swi
%     function part_data = writePartData(head)
%       part_data = [ ...
%         (0 : length(head.dp_H__Pa) - 1)', ...
%         1.0e-5 * head.dp_H__Pa', ...
%         60.0 * head.v_c__m_s'];
%     end
%     
%     % view system data table of head
%     % rev.20170413/swi
%     function sd_head = viewSystemData(head)
%       sd_head = table();
%       sd_head.gas = class(head.gas);
%       sd_head.max_dp_H__bar = 1.0e-5 * head.dp_H_max__Pa;
%     end
%     
%     
%   end
%   
%   % public class operations: computation
%   methods (Access = public)
%     
%     % compute minimum constraint of feed
%     % rev.20170405/swi
%     function [inf_v_c__m_s, sup_v_c__m_s] = computeConstraintForFeed(head, ...
%         max_v_XO__m_s, ref_v_c__m_s)
%       % compute minimum constraint for feed
%       inf_v_c__m_s = max_v_XO__m_s ...
%         .* (head.nozzle.mu_j .* head.q_H) .^ 0.25 ...
%         / head.max_v_c__m_s;
%       head.inf_v_c__m_s = inf_v_c__m_s;
%       % compute maximum constraint for feed
%       sup_v_c__m_s = (ref_v_c__m_s ...
%         .* (head.nozzle.mu_j .* head.q_H) .^ 0.25 ...
%         + max_v_XO__m_s) ...
%         .* (head.nozzle.mu_j .* head.q_H) .^ 0.25 ...
%         / head.max_v_c__m_s;
%       head.sup_v_c__m_s = sup_v_c__m_s;
%     end
%     
%     % compute suggestion for part data values within constraint region
%     % rev.20170710/swi
%     function cstr_dp_H__Pa = computePartDataCompatibility(head, ...
%         min_eta_j, max_eta_j)
%       cstr_dp_H__Pa = head.dp_H__Pa;
%       % compute pressure constraint
%       inf_dp_H__Pa = min_eta_j .^ 2.0 ./ head.nozzle.mu_j * head.dp_H_max__Pa;
%       sup_dp_H__Pa = max_eta_j .^ 2.0 ./ head.nozzle.mu_j * head.dp_H_max__Pa;
%       % compute fraction
%       q = (head.dp_H__Pa - inf_dp_H__Pa) ...
%         ./ (sup_dp_H__Pa - inf_dp_H__Pa);
%       % check pressure constraint
%       is_constraint = ...
%         any(head.dp_H__Pa < inf_dp_H__Pa) ...
%         | any(sup_dp_H__Pa < head.dp_H__Pa);
%       % show error in case of constraint
%       if is_constraint
%         cstr_dp_H__Pa = ...
%           min(q) * sup_dp_H__Pa + (1.0 - min(q)) * inf_dp_H__Pa;
% %           max(inf_dp_H__Pa, ...
% %           min(sup_dp_H__Pa, head.dp_H__Pa));
%         disp(num2str(1.0e-5 * cstr_dp_H__Pa', '%10.5e\n'));
%         err_msg = ['Pressure data input out of constraint.' ...
%           'Use displayed values in [bar] instead.'];
%         error(err_msg);
%       end      
%     end
% 
%     
%     % make part data arrays over outline
%     % rev.20170522/swi
%     function makePartDataArray(head, u_c, max_v_c__m_s)
%       assert(all(head.dp_H__Pa < head.dp_H_max__Pa));
%       assert(all(head.v_c__m_s < max_v_c__m_s));
%       head.dp_H__Pa(head.dp_H__Pa < 0.0) = 0.0;
%       head.v_c__m_s(head.v_c__m_s < head.v_c_min__m_s) = head.v_c_min__m_s;
%       % make unit overpressure array
%       unit_offset = min(head.dp_H__Pa) / head.dp_H_max__Pa;
%       unit_amplitude = ...
%         (max(head.dp_H__Pa) - min(head.dp_H__Pa)) / head.dp_H_max__Pa;
%       head.q_H = unit_offset + unit_amplitude * u_c;
%       % make unit feed array
%       unit_offset = min(head.v_c__m_s) / max_v_c__m_s;
%       unit_amplitude = ...
%         (max(head.v_c__m_s) - min(head.v_c__m_s)) / max_v_c__m_s;
%       head.q_c = unit_offset + unit_amplitude * u_c;
%     end
%     

% end of Head.m

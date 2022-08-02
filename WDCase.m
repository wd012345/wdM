classdef WDCase < handle
    % WDCase   interface of computation cases
    % 
    % The WDCase class is the interface of all application cases of
    % scientific computation. Examples of computation cases are statistics 
    % reports, risk analysis, finite element computations etc.
    %
    %
    % Properties (protected)
    % ----------------------
    % data              case data in structures
    % wd_io             zero or more data I/O interfaces
    % wd_object         real world objects of this case
    % wd_postprocessor  campaign-specific postprocessor object
    % wd_preprocessor   campaign-specific preprocessor object
    %
    %
    % Abstract Methods
    % ----------------
    % clone             call clone() of data I/O objects
    % commit            call commit() of data I/O objects
    % postprocess       postprocessing method of this case
    % preprocess        preprocessing method of this case
    % runDatatools      run specific datatools to compute case data
    % visualize         call visualize() of WDGraphic objects
    %
    % ________________________________________________________________________
    %
    % ToDos:
    % 1. ...
    %
    % Author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    properties (Access = protected)
        data                    % case data in structured tables
        wd_io                   % zero or more data I/O interfaces
        wd_object               % real world objects of this case
        wd_postprocessor        % case-specific postprocessor object
        wd_preprocessor         % case-specific preprocessor object
        wd_graphic              % case-specific graphic objects
    end
    
    methods (Abstract)
        clone()                 % call clone() of data I/O objects
        commit()                % call commit() of data I/O objects
        postprocess()           % postprocessor object of campaign
        preprocess()            % preprocessor object of campaign
        runDatatools()          % run campaign-specific datatools
        visualize()             % call visualize() of WDGraphic objects
    end%methods
end%classdef


% end of module wdM.wdtM.WDCase

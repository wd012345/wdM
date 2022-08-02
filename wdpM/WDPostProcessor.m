classdef WDPostProcessor < handle
    % WDPostProcessor   generic postprocessor interface
    % The postprocessing interface inherits the postprocess method which must
    % be implemented by the various postprocessors. The postprocess method 
    % does the postprocess of case data before saving and exiting the case.
    %
    %
    %   Properties
    %   ----------
    %
    %
    %   Methods
    %   -------
    %
    %
    %   Methods (abstract)
    %   ------------------
    %   postprocess             postprocess campaign data
    %
    % ________________________________________________________________________
    %
    % ToDo
    % ----
    % 1. ...
    %
    % Author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
        
    methods (Abstract)
        postprocess()           % postprocess case data
    end%methods
    
    
end%classdef


% end of module WDPostProcessor

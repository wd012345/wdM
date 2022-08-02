classdef (Abstract) WDPreProcessor < handle
    % WDPreProcessor   preprocessor interface
    % The abstract class WDPreProcessor is an interface to various
    % preprocessor implementations which implement the preprocessing of
    % user specified data.
    %
    % Properties (protected)
    % ----------------------
    % raw               input raw data to preprocess
    %
    % Methods
    % -------
    % WDPreProcessor    construct this preprocessor interface
    %
    %
    % ToDos:            ...
    %
    % Authors:          Stefan Wittwer, info@wittwer-datatools.ch
    
% properties (protected) %
%     properties (Access = protected)
%         raw        % input raw data to preprocess
%     end%properties
%     
% % constructor *
%     methods
%         function wd_prp = WDPreProcessor(raw)
%             % WDPreProcessor   construct this preprocessor interface
%             wd_prp.raw = raw;
%         end%function
%     end%methods
    
% abstract methods %
    methods (Abstract = true)
        preprocess(wd_prp);        
    end%methods
end%classdef


% end of module WDPrePreprocessor

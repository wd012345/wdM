classdef WDGraphic < handle
    % WDGraphic  data visualizing interface
    % 
    %
    % Properties
    % ----------
    % N_g               length of graphic array
    % ui                reference to case UI (app)
    % graph             graphic array with figure handles
    % video             case video
    %
    %
    % Methods
    % -------
    % WDGraphic         construct WDGraphic object
    %
    %
    % Methods (abstract)
    % ------------------
    % visualize         inherited visualize method declaration
    %
    % ________________________________________________________________________
    %
    % ToDo
    % 1. ...
    %
    % Author: Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    properties (Access = protected)
        N_g
        graph
        ui
        video
    end
    
    methods
    % constructor %
        function wd_graphic = WDGraphic(N_g)
            % WDGraphic  construct WDGraphic object
            %
            % Construct this WDGraphic object and initializes it.
            %
            % Arguments:
            % N_g               length of graphic array
            %
            % Returns:
            % wd_graphic        newly constructed WDGraphic object
            arguments
                N_g (1, 1) int32 = 0
            end%arguments
            wd_graphic.N_g = max(0, N_g);
            if 0 < wd_graphic.N_g
                wd_graphic.graph = gobjects(wd_graphic.N_g, 1);
            end%if
        end%function
        
        % get methods %
        function graph = getGraph(wd_graphic)
            graph = wd_graphic.graph;
        end%function
        function ui = getUI(wd_graphic)
            ui = wd_graphic.ui;
        end%function
        function video = getVideo(wd_graphic)
            video = wd_graphic.video;
        end%function
    end%methods
    
    methods (Abstract = true)
        visualize(data, n_g)    % inherited visualize method declaration
    end
end


% end of module WDGraphic

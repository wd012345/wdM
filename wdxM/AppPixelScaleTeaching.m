classdef AppPixelScaleTeaching < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PixelScaleTeachingAppUIFigure  matlab.ui.Figure
        DoneButton                     matlab.ui.control.Button
        TextAreaLabel                  matlab.ui.control.Label
        TextArea                       matlab.ui.control.TextArea
        UIAxes                         matlab.ui.control.UIAxes
        UIAxes2                        matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        d_n % nozzle diameter
        f_n % selected frame index
        S_n__au % pixel scale [a.u.]
        p_0 % p coordinate of tool center
        q_0 % q coordinate of tool center
    end
    
    properties (Access = private)
        x % Description
        y % Description
        R                       % nominal frame positions on part contour
        c                       % image sequence
        hline % Description
        himage % Description
        d_ref__au % Description
    end
    
    methods (Access = public)
        
        function plotContour(app, R)
            app.hline = plot(app.UIAxes, R.X__m(:, 1), R.Y__m(:, 1), 'o');
            app.hline.ButtonDownFcn = @app.selectPosition;
            app.x = [];
            app.y = [];
            app.R = R;
            app.TextArea.Value = ...
                'Proceed by selecting a position on the left graph';
        end
        
        function selectPosition(app, ~, event)
            if isempty(app.x) && isempty(app.y)
                X = event.IntersectionPoint(1);
                Y = event.IntersectionPoint(2);
                % evaluate and append frame number
                f_x = find( ...
                    round(1.0e6 * X) == ...
                    round(1.0e6 * app.R.X__m(:, 1)));
                f_y = find( ...
                    round(1.0e6 * Y) == ...
                    round(1.0e6 * app.R.Y__m(:, 1)));
                f_xy = (f_x == f_y);
                if numel(f_xy) == numel(f_x)
                    if numel(app.f_n) == numel(app.d_n)
                        app.f_n = [app.f_n; f_x(f_xy)];
                    elseif numel(app.f_n) == numel(app.d_n) + 1
                        app.f_n(end) = f_x(f_xy);
                    end%if
                else
                    if numel(app.f_n) == numel(app.d_n)
                        app.f_n = [app.f_n; f_y(f_xy)];
                    elseif numel(app.f_n) == numel(app.d_n) + 1
                        app.f_n(end) = f_y(f_xy);
                    end%if
                end%if
                app.x = [];
                app.y = [];                
                % display selected frame
                app.himage = imagesc(app.UIAxes2, app.c{app.f_n(end), 1});
                app.himage.ButtonDownFcn = @app.selectPositions;
                app.UIAxes2.XLabel.String = 'q [pixel]';
                app.UIAxes2.XLim = [1 size(app.c{app.f_n(end), 1}, 2)];
                app.UIAxes2.YLabel.String = 'p [pixel]';
                app.UIAxes2.YLim = [1 size(app.c{app.f_n(end), 1}, 1)];                
                app.TextArea.Value = { ...
                    'Selected position:'; ...
                    ['X = ' num2str(1.0e6 * X) ' {\mu}m']; ...
                    ['Y = ' num2str(1.0e6 * Y) ' {\mu}m']; ...
                    ['at frame ' num2str(app.f_n(end))]; ...
                    'Proceed with selecting first point on nozzle diagonal';
                    'or select another position on the contour';
                    'or press Done to exit.'};
            else
                app.x = [];
                app.y = [];
                app.TextArea.Value = { ...
                    'Missing second position on nozzle diagonal.'; ...
                    'Reset to proceed with selecting contour position.'};
            end
        end
        
        function selectPositions(app, ~, event)
            if isempty(app.x) && isempty(app.y)
                app.x = event.IntersectionPoint(1);
                app.y = event.IntersectionPoint(2);
                app.TextArea.Value = { ...
                    'First point on nozzle diagonal has been selected'; ...
                    ['x = ' num2str(app.x)]; ...
                    ['y = ' num2str(app.y)]; ...
                    'Proceed with selecting opposite point on '; ...
                    'nozzle diagonal'};
            else
               if numel(app.f_n) == numel(app.d_n)
                   app.d_n(end) = sqrt( ...
                    (app.x - event.IntersectionPoint(1)) ^ 2.0 + ...
                    (app.y - event.IntersectionPoint(2)) ^ 2.0);
               elseif numel(app.f_n) == numel(app.d_n) + 1
                   app.d_n = [app.d_n; sqrt( ...
                       (app.x - event.IntersectionPoint(1)) ^ 2.0 + ...
                       (app.y - event.IntersectionPoint(2)) ^ 2.0)];
               end%if
                if ~isempty(app.d_ref__au)
                    app.S_n__au = [app.S_n__au; app.d_ref__au / app.d_n(end)];
                end%if
                app.q_0 = [app.q_0; ...
                    0.5 * (app.x + event.IntersectionPoint(1))];
                app.p_0 = [app.p_0; ...
                    0.5 * (app.y + event.IntersectionPoint(2))];
                app.x = [];
                app.y = [];
                app.TextArea.Value = { ...
                    'Done with selection of nozzle diagonal.'; ...
                    ['d_n = ' num2str(app.d_n(end))]; ...
                    ['S_n__au = ' num2str(app.S_n__au(end))]; ...
                    ['p_0 = ' num2str(app.p_0(end))]; ...
                    ['q_0 = ' num2str(app.q_0(end))]; ...
                    'Proceed with selecting another position ';
                    'or press Done to exit.'};
            end            
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, c, d_ref__au)
            app.d_ref__au = d_ref__au;  % set reference nozzle diameter
            app.c = c;                  % load image sequence
            app.d_n = [];               % initialize properties
            app.f_n = [];
            app.S_n__au = [];
            app.p_0 = [];
            app.q_0 = [];
        end

        % Button pushed function: DoneButton
        function DoneButtonPushed(app, ~)
            if numel(app.f_n) == numel(app.d_n) + 1
                app.f_n(end) = [];
            end%ir
            if numel(app.d_n) == numel(app.f_n)
                s = struct( ...
                    'd_n', app.d_n, ...
                    'f_n', app.f_n, ...
                    'S_n__au', app.S_n__au, ...
                    'p_0', app.p_0, ...
                    'q_0', app.q_0);
                save(['.' filesep 'tmp' filesep 'PixelScalesTmp.mat'], ...
                    '-struct', 's');
                app.delete();
            else
                app.TextArea.Value = ...
                    'Select missing nozzle diagonal before closing.';
            end%if
        end

        % Close request function: PixelScaleTeachingAppUIFigure
        function PixelScaleTeachingAppUIFigureCloseRequest(app, ~)
            if numel(app.f_n) == numel(app.d_n) + 1
                app.f_n(end) = [];
            end%ir
            if numel(app.d_n) == numel(app.f_n)
                s = struct( ...
                    'd_n', app.d_n, ...
                    'f_n', app.f_n, ...
                    'S_n__au', app.S_n__au, ...
                    'p_0', app.p_0, ...
                    'q_0', app.q_0);
                save(['.' filesep 'tmp' filesep 'PixelScalesTmp.mat'], ...
                    '-struct', 's');
                app.delete();
            else
                app.TextArea.Value = ...
                    'Select missing nozzle diagonal before closing.';
            end%if
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PixelScaleTeachingAppUIFigure and hide until all components are created
            app.PixelScaleTeachingAppUIFigure = uifigure('Visible', 'off');
            app.PixelScaleTeachingAppUIFigure.Position = [100 100 640 480];
            app.PixelScaleTeachingAppUIFigure.Name = 'Pixel Scale Teaching App';
            app.PixelScaleTeachingAppUIFigure.CloseRequestFcn = createCallbackFcn(app, @PixelScaleTeachingAppUIFigureCloseRequest, true);

            % Create DoneButton
            app.DoneButton = uibutton(app.PixelScaleTeachingAppUIFigure, 'push');
            app.DoneButton.ButtonPushedFcn = createCallbackFcn(app, @DoneButtonPushed, true);
            app.DoneButton.FontName = 'Verdana';
            app.DoneButton.Position = [99 143 111 45];
            app.DoneButton.Text = 'Done';

            % Create TextAreaLabel
            app.TextAreaLabel = uilabel(app.PixelScaleTeachingAppUIFigure);
            app.TextAreaLabel.HorizontalAlignment = 'right';
            app.TextAreaLabel.FontName = 'Verdana';
            app.TextAreaLabel.Position = [10 84 62 22];
            app.TextAreaLabel.Text = 'Text Area';

            % Create TextArea
            app.TextArea = uitextarea(app.PixelScaleTeachingAppUIFigure);
            app.TextArea.Editable = 'off';
            app.TextArea.FontName = 'Verdana';
            app.TextArea.Position = [87 13 541 95];
            app.TextArea.Value = {'''Select position on trajectory graph.'''};

            % Create UIAxes
            app.UIAxes = uiaxes(app.PixelScaleTeachingAppUIFigure);
            title(app.UIAxes, 'Select Position')
            xlabel(app.UIAxes, 'X [a.u.]')
            ylabel(app.UIAxes, 'Y [a.u.]')
            app.UIAxes.PlotBoxAspectRatio = [1.32432432432432 1 1];
            app.UIAxes.FontName = 'Verdana';
            app.UIAxes.Position = [16 264 243 202];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.PixelScaleTeachingAppUIFigure);
            title(app.UIAxes2, 'Select Nozzle')
            xlabel(app.UIAxes2, 'q')
            ylabel(app.UIAxes2, 'p')
            app.UIAxes2.PlotBoxAspectRatio = [1.03040540540541 1 1];
            app.UIAxes2.FontName = 'Verdana';
            app.UIAxes2.ClippingStyle = 'rectangle';
            app.UIAxes2.Position = [276 116 352 350];

            % Show the figure after all components are created
            app.PixelScaleTeachingAppUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AppPixelScaleTeaching(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.PixelScaleTeachingAppUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.PixelScaleTeachingAppUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PixelScaleTeachingAppUIFigure)
        end
    end
end
classdef WDCamera < WDObject
    % WDCamera    generic camera object
    %
    % Definition of generic camera objects used to capture single images
    % or image series (videos, series of frames).    
    %
    %
    % Properties (private)
    % --------------------
    % f_B                       basic frame index over all video streams
    %                           (f_B == idxFrameList)
    % wd_video                  video streams
    %
    %
    % Methods (public)
    % ----------------
    % WDCamera                  construct generic camera object
    %
    % getHeightOfFrames         get number or pixel rows of video streams
    % getIndexOfFrames          get basic frame index over video streams
    % getNumberOfFrames         get number of stream frames
    % getNumberOfStreams        get number of streams
    % getNumberOfValuesPerPixel get number of values per pixel
    % getStream                 get specified or all video streams
    % getVideo                  get specified or all video objects
    % getWidthOfFrames          get number or pixel colums
    % setIndexOfFrames          set base index of frames
    % setVideo                  set one or more video streams
    %
    % computeDistribution       compute probabilities of pixel values
    % computeThreshold          compute threshold of image sequence
    % generateBinary            generate binary images from threshold
    % normalizeToUnity          normalize image sequence to unity
    % squeezeToChannel          average to single channel image streams
    % trimPixelValues           trim pixel values to given range
    % ________________________________________________________________________
    %
    % todo:
    % 1. set number of values per pixel correctly (line 201)
    % 2. ...
    %
    % author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    
%% properties (private) %%
    properties (Access = private)
        f_B                     % basic frame index over all video streams
                                % (f_B == idxFrameList)
        wd_video                % video streams
    end%properties
    
    
%% methods (public)
    methods
    %% constructor %%
        function wd_camera = WDCamera(name, uid)
            % WDCamera    construct generic camera object
            %
            % Construct this generic camera object and initialize the 
            % property of the video streams as empty cell array.
            % 
            % Returns:
            % wd_camera         this WDCamera object
            wd_camera = wd_camera@WDObject(name, uid);
            wd_camera.wd_video = {};
        end%function
        
    %% get and set methods %%
        function N_p = getHeightOfFrames(wd_camera)
            % getHeightOfFrames    get number or pixel rows of video streams
            %
            % Return the number of pixel rows of the video streams.
            %
            % Returns:
            % N_p               frame height or number of pixel rows
            N_s = wd_camera.getNumberOfStreams();
            N_p = cell(N_s, 1);
            for s = 1:N_s
                N_p{s, 1} = wd_camera.wd_video{s, 1}.getHeight();
            end%for
        end%function
        function f_B = getIndexOfFrames(wd_camera)
            % getIndexOfFrames    get base index of frames
            f_B = wd_camera.f_B;
        end%function
        function [N_f, N_s] = getNumberOfFrames(wd_camera)
            % getNumberOfFrames    get length of image series
            N_s = wd_camera.getNumberOfStreams();
            N_f = zeros(N_s, 1);
            for s = 1:N_s
                N_f(s, 1) = wd_camera.wd_video{s, 1}.getNumberOfFrames();
            end%for
        end%function
        function N_s = getNumberOfStreams(wd_camera)
            % getNumberOfStreams    get number of streams
            N_s = numel(wd_camera.wd_video);
        end%function
        function [N_b, N_s] = getNumberOfValuesPerPixel(wd_camera)
            % getNumberOfValuesPerPixel    get number of values per pixel
            N_s = wd_camera.getNumberOfStreams();
            N_b = zeros(N_s, 1);
            for s = 1:N_s
                N_b(s, 1) = wd_camera ...
                    .wd_video{s, 1} ...
                    .getNumberOfValuesPerPixel();
            end%for
        end%function
        function wd_image = getStream(wd_camera, varargin)
            % getStream    get specified or all video streams
            %
            % Return one or more video streams as sequence of WDImage objects.
            % 
            % Arguments:
            % varargin          index of stream or empty
            %
            % Returns:
            % wd_images         specified or all video streams
            if nargin == 1
                [N_f, N_s] = wd_camera.getNumberOfFrames();
                wd_image = cell(max(N_f), N_s);
                for s = 1:N_s
                    wd_image(:, s) = wd_camera.wd_video{s, 1}.getStream();
                end%for
            else
                wd_image = wd_camera.wd_video{varargin{1, 1}, 1}.getStream();
            end%if
        end%function
        function wd_video = getVideo(wd_camera, varargin)
            % getVideo    get specified or all video objects
            %
            % Arguments:
            % varargin          index of video object or empty
            %
            % Returns:
            % wd_video          specified or all video objects
            if nargin == 1
                wd_video = wd_camera.wd_video;
            else
                s = max(1, min(numel(wd_camera.wd_video), varargin{1, 1}));
                wd_video = wd_camera.wd_video{s, 1};
            end%if
        end%function
        function N_q = getWidthOfFrames(wd_camera)
            % getWidthOfFrames    get number or pixel colums
            %
            % Return the number of pixel colums of the video streams.
            % image or image series.
            %
            % Returns:
            % N_q               frame width or number of pixel colums
            N_s = wd_camera.getNumberOfStreams();
            N_q = cell(N_s, 1);
            for s = 1:N_s
                N_q{s, 1} = wd_camera.wd_video{s, 1}.getWidth();
            end%for
        end%function
        function setIndexOfFrames(wd_camera, f_B)
            % setIndexOfFrames    set base index of frames
            %
            % Arguments:
            % f_B               specified basic frame index
            wd_camera.f_B = f_B;
        end%function
        function setVideo(wd_camera, video, name)
            % setVideo     set one or more video streams
            %
            % Set one or more named video objects.
            %
            % Arguments:
            % video             video stream(s)
            % name              name of each given video stream
            N_s = length(name);
            wd_camera.wd_video = cell(N_s, 1);
            for s = 1:N_s
               wd_camera.wd_video{s, 1} = WDVideo("", name(s));
               wd_camera ...
                   .wd_video{s, 1} ...
                   .setFrameRate(video.info(s, 1).f__Hz);
%                wd_camera ...
%                    .wd_video{s, 1} ...
%                    .setHeight(video.info(s, 1).N_p);
               wd_camera ...
                   .wd_video{s, 1} ...
                   .setNumberOfValuesPerPixel(2^10);
               wd_camera.wd_video{s, 1}.setStream(video.stream(:, s));
%                wd_camera ...
%                    .wd_video{s, 1} ...
%                    .setWidth(video.info(s, 1).N_q);
            end%for
        end%function
        
    %% computational methods %%
        function [P_h, N_b, N_s] = computeDistribution(wd_camera)
            % computeDistribution    compute probabilities of pixel values
            %
            % Compute probability distribution of pixel values of the
            % entire image series.
            %
            % Returns:
            % P_h       probabilities of pixel values
            % N_b       number of pixel values
            % N_s       number of video streams
            N_s = wd_camera.getNumberOfStreams();
            P_h = cell(N_s, 1);
            N_b = zeros(N_s, 1);
            for s = 1:N_s
                [P_h{s, 1}, N_b(s, 1)] = wd_camera ...
                    .wd_video{s, 1} ...
                    .computeDistribution();
                P_h{s, 1} = P_h{s, 1} ./ sum(P_h{s, 1}, 2); % normalize
            end%for
        end%function
        function [C_t, P_h] = computeThreshold(wd_camera, T_C)
            % computeThreshold    compute threshold of image sequence
            %
            % Compute image threshold of each image frame in the image
            % sequence using the method given in the argument. If needed 
            % compute the probability distribution of the pixel values.
            %
            % Arguments:
            % T_C               generating method of image threshold
            %
            % Returns:
            % C_0               image thresholds
            % P_h               probability distribution
            N_s = wd_camera.getNumberOfStreams();
            C_t = cell(N_s, 1);
            P_h = cell(N_s, 1);
%             [P_h, N_b, ~] = wd_camera.computeDistribution();
            for s = 1:N_s
                [C_t{s, 1}, P_h{s, 1}] = wd_camera ...
                    .getVideo(s) ...
                    .computeThreshold(T_C(s, 1));
%                 switch T_C(s, 1)
%                     case "0<C"
%                         C_0{s, 1} = zeros(N_f(s, 1), 1);
%                     case "1<C"
%                         C_0{s, 1} = ones(N_f(s, 1), 1);
%                     case "(2^8-1)<C"
%                         C_0{s, 1} = ones(N_f(s, 1), 1) * (2^8 - 1);
%                     case "(N_b-2)<C"
%                         C_0{s, 1} = ones(N_f(s, 1), 1) * (N_b(s, 1) - 2);
%                     case "1.moment"
%                         C_0{s, 1} = sum((0:N_b(s, 1) - 1) .* P_h{s, 1}, 2);
%                     case "2.moment"
%                         mu = sum((0:N_b(s, 1) - 1) .* P_h{s, 1}, 2);
%                         C_0{s, 1} = sqrt(sum( ...
%                             ((0:(N_b(s, 1) - 1)) - mu) .^ 2.0 .* ...
%                             P_h{s, 1}, 2));
%                     case "otsu"
%                         C_0{s, 1} = zeros(N_f(s, 1), 1);
%                         for f = 1:N_f(s, 1)
%                             C_0{s, 1}(f, 1) = wd_camera ...
%                                 .wd_video{s, 1} ...
%                                 .getFrame(f) ...
%                                 .computeOtsu(N_b(s, 1));
%                         end%for
%                     case "C(P=5%)"
%                         C_0{s, 1} = zeros(N_f(s, 1), 1);
%                         for f = 1:N_f(s, 1)
%                             C_0{s, 1}(f, 1) = find( ...
%                                 0.05 < cumsum(P_h{s, 1}(f, :)), 1, 'first');
%                         end%for
%                     case "C(P=10%)"
%                         C_0{s, 1} = zeros(N_f(s, 1), 1);
%                         for f = 1:N_f(s, 1)
%                             C_0{s, 1}(f) = find( ...
%                                 0.10 < cumsum(P_h{s, 1}(f, :)), 1, 'first');
%                         end%for
%                     case "C(P=25%)"
%                         C_0{s, 1} = zeros(N_f(s, 1), 1);
%                         for f = 1:N_f(s, 1)
%                             C_0{s, 1}(f, 1) = find( ...
%                                 0.25 < cumsum(P_h{s, 1}(f, :)), 1, 'first');
%                         end%for
%                     case "F<=1/C"
%                         C_0{s, 1} = zeros(N_f(s, 1), 1);
%                         for f = 1:N_f(s, 1)
%                             c_0 = find( ...
%                                 (1.0 ./ double(2:size(P_h{s, 1}, 2))) < ...
%                                 P_h{s, 1}(f, 2:end), 1, 'first');                            
%                             if ~isempty(c_0)
%                                 C_0{s, 1}(f, 1) = c_0;
%                             end%if                           
%                         end%for
%                     otherwise
%                         throw(MException( ...
%                             "WDCamera:UnsupportedImageThreshold", ...
%                             "Image threshold method " + ...
%                             num2str(T_C(s, 1)) + " not supported."));
%                 end%switch
            end%for
        end%function
%         function [dC, dC_df] = computeVariation(wd_camera)
%             % computeVariation    compute variation of stream frames
%             %
%             % Compute the standard deviation and the first derivative of
%             % the pixel values of the image frames in the image sequence.
%             [N_f, N_s] = wd_camera.getNumberOfFrames();
%             C = 
%             wd_image = wd_camera.getStream();
%             
%             % compute standard deviation of frames
%             dC = std(C, 0, 3);
%             
%             % compute difference of frames
%             dC_df = diff(C, 1, 3);
%         end%function
        function binary = generateBinary(wd_camera, C_t, isaPositive)
            % generateBinary    generate binary images from threshold
            %
            % Generate binary images of each frame of the image series by
            % by setting pixel value below and above threshold to
            % zeros and ones, respectively.
            %
            % Arguments:
            % C_t               image threshold
            % isaPositive       if true generate positive binary image
            %
            % Returns:
            % binaries          binary image series
            [N_f, N_s] = wd_camera.getNumberOfFrames();
            binary = cell(max(N_f), N_s);
            for s = 1:N_s
                b = wd_camera ...
                    .wd_video{s, 1} ...
                    .generateBinary(C_t{s, 1}, isaPositive);
                binary(1:numel(b), s) = b;
            end%for
        end%function
        function c = normalizeToUnity(wd_camera)
            % normalizeToUnity    normalize image sequence to unity
            %
            % Normalize the pixel values of image sequences to the unit 
            % interval [0; 1].
            %
            % Returns:
            % c                 channels array with unit pixel values
            N_s = wd_camera.getNumberOfStreams();
            c = cell(N_s, 1);
            for s = 1:N_s
                c{s, 1} = wd_camera ...
                    .getVideo(s) ...
                    .normalizeToUnity();
            end%for
        end%function
        function G = squeezeToChannel(wd_camera)
            % squeezeToChannel    average to single channel image streams
            %
            % Compute the average value of the channels of each image of all
            % camera streams.
            %
            % Returns:
            % G                 single channel image streams
            N_s = wd_camera.getNumberOfStreams();
            G = cell(N_s, 1);
            for s = 1:N_s
                G{s, 1} = wd_camera.wd_video{s, 1}.squeezeToChannel();
            end%for
        end%function
        function C = trimPixelValues(wd_camera, C_min, C_max)
            % trimPixelValues    trim pixel values to given range
            %
            % Restrict pixel values in each channel of each image frame in
            % the image series to the minimum and maximum value given in
            % the argument list. The image series are updated with those bounded
            % to interval.
            %
            % Arguments:
            % C_min     minimum pixel value
            % C_max     maximum pixel value
            %
            % Returm:
            % C         channels trimmed to interval
            N_s = wd_camera.getNumberOfStreams();
            C = cell(N_s, 1);
            for s = 1:N_s
                c = wd_camera ...
                    .wd_video{s, 1} ...
                    .trimPixelValues(C_min, C_max{s, 1});
                C(1:numel(c), s) = c;
            end%for
        end%function
        

        
        
    %% other methods ...
%         function [N_p, N_q, N_s] = getSizeOfImages(wd_camera)
%             % getSizeOfImages    get number or pixel rows and colums
%             %
%             % Return the number of pixel rows and pixel colums of the
%             % image or image series.
%             %
%             % Returns:
%             % N_p       number of pixel rows
%             % N_q       number of pixel colums
%             % N_s       number of video streams
%             N_s = wd_camera.getNumberOfStreams();
%             N_p = zeros(N_s, 1);
%             N_q = zeros(N_s, 1);
%             for s = 1:N_s
%                 wd_image = wd_camera.wd_video{s, 1}.getFrame(1);
%                 [N_p(s, 1), N_q(s, 1), ~] = wd_image.getSize();
%             end%for
%         end%function
%         function url = getURL(wd_cam)
%             % getURL   return url of stream and image series
%             url = cell(wd_cam.N_s, 1);
%             for s = 1:wd_cam.N_s
%                 url{s, 1} = wd_cam.wd_video{s, 1}.getURL();
%             end%for
%         end%function
%         function setImageSeries(camera, images)
%             % setImageSeries   assign image series of each stream
%             % Assign image series of each stream, organized in cell colums.
%             camera.images = images;
%         end%function
%         function setNumberOfImages(camera, N_f)
%             % setNumberOfImages   assign number of images of each stream
%             camera.N_f = N_f;
%         end%function
%         function setNumberOfStreams(camera, N_s)
%             % setNumberOfStreams   assign numbr of streams
%             camera.N_s = N_s;
%         end%function
%         function setSizeOfImages(camera, N_p, N_q)
%             % setSizeOfImages   assign number of pixel rows and colums
%             % Assign the number of pixel rows (image height) and the number of
%             % pixel colums (image width) for each stream.
%             if ~isempty(N_p)
%                 camera.N_p = N_p;
%             end%if
%             if ~isempty(N_q)
%                 camera.N_q = N_q;
%             end%if
%         end%function
%         function setStreamSeries(camera, streams)
%             % setStreamSeries   assign stream series properties
%             camera.streams = streams;
%         end%function
%         function setURL(camera, url)
%             % setURL   assign url of stream and image series
%             camera.url = url;
%         end%function
%         function setValuesPerPixel(camera, N_b)
%             % setValuesPerPixel   assign number of values per pixel
%             camera.N_b = max(1,N_b);
%         end%function  
%         function frame_data = preprocessData(camera, wd_preprocessor, ...
%                 fb_parameter, t__s, t_camb__s)
%             % preprocessData   preprocess input data of this camera object
%             % Preprocess input data of this camera object by using the
%             % given preprocessor instance and additional parameters.
%             
%             % load and preprocess input data
%             [camera.images, wdfs] = wd_preprocessor.preprocessCamera( ...
%                 fb_parameter, t__s, t_camb__s);
%             
%             % assign preprocessed input data
%             camera.N_b = wdfs.N_b;
%             camera.N_f = wdfs.N_f;
%             camera.N_p = wdfs.N_p;
%             camera.N_q = wdfs.N_q;
%             camera.f__Hz = wdfs.f__Hz;
%             camera.frameListArray = wdfs.frameListArray;
%             camera.frameListIllu = wdfs.frameListIllu;
%             camera.idxFrameList = wdfs.idxFrameList;
%             camera.idxReducedArray = wdfs.idxReducedArray;
% %             camera.url = wdfs.url;
%             
%             % return frame data
%             frame_data = wdfs;
%         end%function
%         function aperture = computeDiametersOfAperture(camera)
%             %
%             %
%             aperture = WDFieldSet('Aperture', struct( ...
%                 'OK', false(camera.N_f, 1), ...
%                 'd_n', zeros(camera.N_f, 1), ...
%                 'p_0', zeros(camera.N_f, 1), ...
%                 'q_0', zeros(camera.N_f, 1)));
%             for f = 1 : camera.N_f
%                 % scale to [0; 255]
%                 C = camera.images{f}.scaleTo(0.0, 255.0);
%                 % apply 5% low pass
%                 C = (C <= 0.05 * 255.0) .* C;
%                 % apply negative incircle mask
%                 [M_c, aperture.p_0(f), aperture.q_0(f), aperture.d_n(f)] = ...
%                     camera.images{f}.maskIncircle(0);
%                 C = M_c .* C;
%                 % compute gradient
%                 [gradC_p, gradC_q] = camera.images{f}.computeGradient( ...
%                     ~(C == 0), 'prewitt');
%                 % filter gradients that make sense
%                 C_tl = (gradC_p > 0.0) & (gradC_q < 0.0);
%                 C_tr = (gradC_p > 0.0) & (gradC_q > 0.0);
%                 C_bl = (gradC_p < 0.0) & (gradC_q < 0.0);
%                 C_br = (gradC_p < 0.0) & (gradC_q > 0.0);
%                 C = (C_tl | C_tr | C_bl | C_br) .* C;
%                 % fit to circle
%                 dr = Inf;
%                 [p_n, q_n] = find(C);
%                 p_0 = aperture.p_0(f);
%                 q_0 = aperture.q_0(f);
%                 while 0.1 < dr
%                     cp = cos(atan2(q_n - q_0, p_n - p_0));
%                     sp = sin(atan2(q_n - q_0, p_n - p_0));
%                     P = mean(p_n); % a
%                     Q = mean(q_n); % a
%                     CP = mean(cp); % b
%                     SP = mean(sp); % b
%                     PCP = mean(p_n .* cp); % c
%                     QSQ = mean(q_n .* sp); % c
%                     CP2 = mean(cp .^ 2.0); % d
%                     SP2 = mean(sp .^ 2.0); % d
%                     p_0new = (CP .* PCP - P .* CP2) ./ (CP .^ 2.0 - CP2);
%                     q_0new = (SP .* QSQ - Q .* SP2) ./ (SP .^ 2.0 - SP2);
%                     dr = abs(p_0 - p_0new) + abs(q_0 - q_0new);
%                     p_0 = p_0new;
%                     q_0 = q_0new;
%                 end%while
%                 d_n = 2.0 * mean(P .* CP - PCP) ./ (CP .^ 2.0 - CP2);
% %                 d_n = 2,0 * mean(P .* SP - QSQ) ./ (SP .^ 2.0 - SP2);
%                 % check for valid fit
%                 OK = (0.9 * min(size(C)) < d_n);
%                 OK = OK & (0.45 * min(size(C)) < p_0);
%                 OK = OK & (0.45 * min(size(C)) < q_0);
%                 aperture.OK(f) = OK;
%                 if OK
%                     aperture.d_n(f) = d_n;
%                     aperture.p_0(f) = p_0;
%                     aperture.q_0(f) = q_0;
%                 end%if
%             end%for
%         end%function
%         
%         function [d_n, f_p__1_m, p_0, q_0, p_c, q_c, M_c, p_n, q_n] = ...
%                 getDiameterOfAperture(camera, d_n__m)
%             [d_n, p_0, q_0, p_c, q_c, M_c, p_n, q_n] = ...
%                 camera.capture.getDiameterOfAperture();
%             f_p__1_m = d_n / d_n__m;
%         end%function
%         
%         function fig_camera = postprocessCamera(camera, title)
%             fig_camera = gobjects(2, 1);
%             axs = gobjects(2, 1);
%             srf = gobjects(2, 1);
%             cbr = gobjects(2, 1);
%             
%             % view mask of stationary form
%             fig_camera(1) = figure(21);
%             fig_camera(1).Color = [1 1 1];
%             fig_camera(1).Visible = 'off';
%             [~, M_s, ~] = camera.capture.getMask();
%             axs(1) = axes('Parent', fig_camera(1));
%             srf(1) = pcolor(axs(1), M_s);
%             srf(1).EdgeColor = 'none';
%             axs(1).FontName = 'Verdana';
%             axs(1).PlotBoxAspectRatio = [1, 1, 1];
%             axs(1).Title.String = {title; 'Stationary Form'};
%             axs(1).Title.FontName = 'Verdana';
%             axs(1).Title.Interpreter = 'none';
%             axs(1).XLabel.String = 'p [pixel]';
%             axs(1).YLabel.String = 'q [pixel]';
%             cbr(1) = colorbar(axs(1));
%             cbr(1).FontName = 'Verdana';
%             
%             % view mask gradient of stationary form
%             fig_camera(2) = figure(22);
%             fig_camera(2).Color = [1 1 1];
%             fig_camera(2).Visible = 'off';
%             [~, ~, gradM_s] = camera.capture.getMask();
%             axs(2) = axes('Parent', fig_camera(2));
%             srf(2) = pcolor(axs(2), abs(gradM_s.p) + abs(gradM_s.q));
%             srf(2).EdgeColor = 'none';
%             axs(2).FontName = 'Verdana';
%             axs(2).PlotBoxAspectRatio = [1, 1, 1];
%             axs(2).Title.String = {title; ...
%                 'Gradient Norm of Stationary Form'};
%             axs(2).Title.FontName = 'Verdana';
%             axs(2).Title.Interpreter = 'none';
%             axs(2).XLabel.String = 'p [pixel]';
%             axs(2).YLabel.String = 'q [pixel]';
%             cbr(2) = colorbar(axs(2));
%             cbr(2).FontName = 'Verdana';
%             
%             fig_camera(1).Visible = 'on';
%             fig_camera(2).Visible = 'on';
%         end%function
    end%methods
end%classdef


% end of module wdM.wdtM.WDCamera

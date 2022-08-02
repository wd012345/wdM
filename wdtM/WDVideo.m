classdef WDVideo < WDObject
    % WDVideo    represent sequence of WDImage objects
    %
    % The class WDVideo represents a sequence or stream of WDImage digital 
    % picture objects, and provides methods to process them.
    %
    %
    % Properties (private)
    % --------------------
    % N_b                       number of values per pixel
    % N_f                       number of frames per stream
    % a_f                       frame sampling index (frameListArray)
    % dt_f__s                   frame sampling time error (absolute value)
    % dtau_f                    frame sampling time error (relative value)
    % f__Hz                     frames per second
    % f_s                       frame index (frameListIllu)
    % f_red                     frame time index (idxReducedArray)
    % isDark                    dark (1) or illuminated (0) frame
    % name                      name of video
    % t_f__s                    frame sampling time
    %
    %
    % Methods
    % -------
    % WDVideo                   construct instance of this WDVideo class
    %
    % getChannel                get sequence of channel tensors
    % getFrame                  get image frame at given index
    % getFrameRate              get frequency of frame capturing
    % getIndexOfFrames          get illuminated or dark frames index
    % getIndexOfFrameTimes      return index of frame sampling times
    % getNumberOfFrames         get number of image frames
    % getNumberOfValuesPerPixel get number of values per pixel
    % getHeight                 get height of image frames in pixel
    % getSize                   get number of rows, colums, channels, frames
    % getStream                 get sequence of image frames
    % getTimeErrorOfFrame       get absolute and relative frame time error
    % getWidth                  get width of image frames in pixel
    % setIndexOfFrames          set illuminated or dark frames index
    % setIndexOfFrameTimes      set index of frame sampling times
    % setFrameRate              set frequency of frame capturing
    % setHeight                 set number of pixel rows
    % setNumberOfValuesPerPixel set number of values per pixel
    % setStream                 set video stream as sequence of image objects
    % setTimeErrorOfFrame       set absolute and relatve frame time error
    % setWidth                  set number of pixel colums
    %
    % computeDistribution       compute probabilities of pixel values
    % computeGradient           compute gradients of image sequence
    % computeRotatedIndex       compute rotated index of image sequence
    % computeThreshold          compute threshold of image sequence
    % generateBinary            generate binary image sequence
    % normalizeToUnity          compute unit norm of image sequence
    % rotateStream              rotate image sequence with specified angle
    % scaleBy                   scale image sequence channel tensors
    % scaleTo                   scale image sequence to specified value range
    % squeezeToChannel          average to single channel image sequence
    % trimPixelValues           trim pixel values of sequence to range
    % ________________________________________________________________________
    %
    % to do:
    % 1. ...
    %
    % author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    
%% properties (private) %%
    properties (Access = private)
        N_b             % number of values per pixel
        a_f             % frame sampling index (frameListArray)
        dt_f__s         % frame sampling time error (absolute value)
        dtau_f          % frame sampling time error (relative value)
        f__Hz           % frames per second
        f_s             % frame index (frameListIllu) of this video stream
        f_red           % frame time index (idxReducedArray)
        wd_image        % series of WDImage frames
        isDark          % dark (1) or illuminated (0)
%         t_f__s          % frame sampling time
%         url             % name of video file
%         frameListArray          % frame list
%         frameListIllu           % list of illuminated frames
%         idxFrameList            % frame index
%         idxReducedArray         % valid frame index
    end%properties
    
    
%% methods (public) %%
    methods
    %% constructor %%
        function wd_video = WDVideo(name, uid)
            % WDVideo    construct instance of this WDVideo class
            %
            % Initialize the image frame sequence of this WDVideo object. 
            % Each image frame of the sequence is a WDImage object.
            %
            % Arguments:
            % name              name of video file
            % uid               uid or url of video
            %
            % Returns:
            % wd_video          this WDVideo object
            wd_video = wd_video@WDObject(name, uid);
        end%function
        
    %% get and set methods %%
        function C = getChannel(wd_video, varargin)
             % getChannel    get sequence of channel tensors
             N_f = wd_video.getNumberOfFrames();
             C = cell(N_f, 1);
             for f = 1:N_f
                 C{f, 1} = wd_video.wd_image{f, 1}.getChannel(varargin);
             end%for
        end%function
        function wd_image = getFrame(wd_video, f)
            % getFrame    get frame image with at given index
            wd_image = wd_video.wd_image{f, 1};
        end%if
        function f__Hz = getFrameRate(wd_video)
            % getFrameRate    get frequency of frame capturing
            f__Hz = wd_video.f__Hz;
        end%function
        function [f_s, isDark] = getIndexOfFrames(wd_video)
            % getIndexOfFrames    get illuminated or dark frames index
            f_s = wd_video.f_s;
            isDark = wd_video.isDark;
        end%function
%         function a_f = getIndexOfFrameSampling(wd_video)
%             % getIndexOfFrameSampling   get frame sampling index
%             a_f = wd_video.a_f;
%         end%function
        function f_red = getIndexOfFrameTimes(wd_video)
            % getIndexOfFrameTimes    return index of frame sampling times
            f_red = wd_video.f_red;
        end%function
        function N_f = getNumberOfFrames(wd_video)
            % getNumberOfFrames     get number of image frames
            N_f = numel(wd_video.wd_image);
        end%function
        function N_b = getNumberOfValuesPerPixel(wd_video)
            % getNumberOfValuesPerPixel    get number of values per pixel
            N_b = wd_video.N_b;
        end%function
        function N_p = getHeight(wd_video)
            % getHeight    get height of image frames in pixel
            N_f = wd_video.getNumberOfFrames();
            N_p = zeros(N_f, 1);
            for f = 1:N_f
                N_p(f, 1) = wd_video.wd_image{f, 1}.getHeight();
            end%for
        end%function
        function [N_p, N_q, N_c, N_f] = getSize(wd_video)
            % getSize    get number of rows, colums, channels, frames
            N_f = wd_video.getNumberOfFrames();
            N_p = zeros(N_f, 1);
            N_q = zeros(N_f, 1);
            N_c = zeros(N_f, 1);
            for f = 1:N_f
                [N_p(f, 1), N_q(f, 1), N_c(f, 1)] = ...
                    wd_video.wd_image{f, 1}.getSize();
            end%for
        end%function
        function wd_image = getStream(wd_video)
            % getStream    get sequence of image frames
            wd_image = wd_video.wd_image;
        end%function
        function [dt_f__s, dtau_f] = getTimeErrorOfFrame(wd_video)
            % getTimeErrorOfFrame    get absolute and relative frame time error
            dt_f__s = wd_video.dt_f__s;
            dtau_f = wd_video.dtau_f;
        end%function
        function N_q = getWidth(wd_video)
            % getWidth    get width of image frames in pixel
            N_f = wd_video.getNumberOfFrames();
            N_q = zeros(N_f, 1);
            for f = 1:N_f
                N_q(f, 1) = wd_video.wd_image{f, 1}.getWidth();
            end%for
        end%function
%         function t_f__s = getTimeOfFrames(wd_video)
%             % getTimeOfFrames   return sampling times of frames
%             t_f__s = wd_video.t_f__s;
%         end%function
%         function setFlagOfIlluminationOff(wd_video, isDark)
%             % setFlagOfIlluminationOff   set flags indicating non-illumination
%             wd_video.isDark = isDark;
%         end%function
        function setFrameRate(wd_video, f__Hz)
            % setFrameRate    set frequency of frame capturing
            wd_video.f__Hz = f__Hz;
        end%function
        function setIndexOfFrames(wd_video, f_s, isDark)
            % setIndexOfFrames   set illuminated or dark frames index
            wd_video.f_s = f_s;
            wd_video.isDark = isDark;
        end%function
        function setIndexOfFrameSampling(wd_video, a_f)
            % setIndexOfFrameSampling   set frame sampling index
            wd_video.a_f = a_f;
        end%function
        function setIndexOfFrameTimes(wd_video, f_red)
            % setIndexOfFrameTimes    set index of frame sampling times
            wd_video.f_red = f_red;
        end%function
        function setNumberOfValuesPerPixel(wd_video, N_b)
            % setNumberOfValuesPerPixel    set number of values per pixel
            wd_video.N_b = N_b;
        end%function
%         function setHeight(wd_video, h)
%             % setHeight    set number of pixel rows
%             wd_video.N_p = h;
%         end%function
        function setStream(wd_video, wd_image)
            % setStream     set video stream as sequence of image objects
            wd_video.wd_image = wd_image;
        end%function
        function setTimeErrorOfFrame(wd_video, dt_f__s, dtau_f)
            % setTimeErrorOfFrame    set absolute and relatve frame time error
            wd_video.dt_f__s = dt_f__s;
            wd_video.dtau_f = dtau_f;
        end%function
%         function setTimeOfFrames(wd_video, t_f__s)
%             % setTimeOfFrameSampling   set sampling times of frames
%             wd_video.t_f__s = t_f__s;
%         end%function
%         function setWidth(wd_video, w)
%             % setWidth    set number of pixel colums
%             wd_video.N_q = w;
%         end%function
        
    %% computational methods %%
        function [P_h, N_b] = computeDistribution(wd_video)
            % computeDistribution    compute probabilities of pixel values
            %
            % Compute probability distribution sequence of the pixel values of
            % each channel.
            %
            % Returns
            % P_h       probability distribution of pixel values
            % N_b       number of values per pixel
            N_b = wd_video.N_b;
            N_f = wd_video.getNumberOfFrames();
            P_h = zeros(N_f, N_b);
            for f = 1:N_f
                [N_h, ~, ~] = wd_video ...
                    .wd_image{f, 1} ...
                    .computeDistribution(N_b);
                P_h(f, :) = N_h / sum(N_h, 2);
            end%for
        end%function
        function [grad_p, grad_q, grad_t] = computeGradient(wd_video)
            % computeGradient    compute gradients of image sequence
            %
            % Compute gradient along rows, colums and along the time axis of
            % the pixel value fields in each channel.
            %
            % Returns:
            % grad_p            pixel value gradient along rows
            % grad_q            pixel value gradient along colums
            % grad_t            gradient of image frames over time axis
            N_f = wd_video.getNumberOfFrames();
            grad_p = cell(N_f, 1);
            grad_q = cell(N_f, 1);
            grad_t = cell(N_f, 1);
            C = cell(N_f, 1);
            for f = 1:N_f
                [grad_p{f, 1}, grad_q{f, 1}, ~] = wd_video ...
                    .wd_image{f, 1} ...
                    .computeGradient([]);
                C{f, 1} = wd_video.wd_image{f, 1}.getChannel();
            end%for
            for f = 1:N_f
                if f == 1 || f == N_f
                   grad_t{f, 1} = 0.5 * C{f, 1} / wd_video.t_f__s(f);
                else
                    grad_t{f, 1} = 0.5 * (C{f + 1, 1} - C{f - 1, 1}) / ...
                        (wd_video.t_f__s(f + 1) - wd_video.t_f__s(f - 1));
                end%if
            end%for
        end%function
        function [p_r, q_r, C_r] = computeRotatedIndex(wd_video, ...
                phi__rad, p_c, q_c)
            % computeRotatedIndex    compute rotated index of image sequence
            %
            % Rotate the index of this image by the angle argument. The 
            % rotation of the index is computed about the given center 
            % coordinates. The rotation is applied to each channel of the 
            % image frames of the sequence.
            %
            % Arguments:
            % phi__rad          rotation angle in [rad]
            % p_c               row coordinate of center
            % q_c               colum coordinate of center
            % 
            % Returns:
            % p_r               rotated row index
            % q_r               rotated colum index
            % C_r               rotated channel tensor
            N_f = wd_video.getNumberOfFrames();
            p_r = cell(N_f, 1);
            q_r = cell(N_f, 1);
            C_r = cell(N_f, 1);
            for f = 1:N_f
                [p_r{f, 1}, q_r{f, 1}, C_r{f, 1}] = wd_video ...
                    .wd_image{f, 1} ...
                    .computeRotatedIndex(phi__rad, p_c(f, 1), q_c(f, 1));
            end%for
        end%function
        function [C_t, P_h] = computeThreshold(wd_video, method)
            % computeThreshold   compute threshold of image sequence
            %
            % Compute image threshold of each image frame in the image
            % sequence. The given method is applied. Also computes
            % the probability distribution of the pixel values if not done
            % yet.
            %
            % Arguments:
            % method            method applied to compute image thresholds
            %
            % Returns:
            % C_t               image thresholds
            N_f = wd_video.getNumberOfFrames();
            [P_h, ~] = wd_video.computeDistribution();
            C_t = zeros(N_f, 1);            
            for f = 1:N_f
                C_t(f, 1) = wd_video ...
                    .wd_image{f, 1} ...
                    .computeThreshold(method, P_h(f, :), wd_video.N_b);
            end%for
        end%function
        function binary = generateBinary(wd_video, C_t, isaPositive)
            % generateBinary    generate binary image sequence
            %
            % Generate binary image sequence using the given image threshold: 
            % pixel values below or equal threshold value become 0, all other
            % pixel values are 1.
            %
            % Arguments:
            % C_t               image threshold
            % isaPositive       if true generate positive binary images
            %
            % Returns:
            % binary            binary image sequence
            N_f = wd_video.getNumberOfFrames();
            binary = cell(N_f, 1);
            for f = 1:N_f
                binary{f, 1} = wd_video ...
                    .wd_image{f, 1} ...
                    .generateBinary(C_t(f, 1), isaPositive);
            end%for
        end%function
        function u = normalizeToUnity(wd_video)
            % normalizeToUnity    compute unit norm of image sequence
            N_f = wd_video.getNumberOfFrames();
            u = cell(N_f, 1);
            for f = 1:N_f
                u{f, 1} = wd_video.getFrame(f).normalizeToUnity();
            end%for
        end%function
        function C_r = rotateStream(wd_video, phi__rad, p_0, q_0)
            % rotateStream    rotate image sequence with specified angle
            %
            % Rotate all images of the sequence with the specified angle. The
            % rotation is done about the center specified by its pixel 
            % coordinates. If the center coordinates are empty, the image 
            % sequence is rotated about its center.
            %
            % Arguments:
            % phi__rad          rotation angle
            % p_0               row coordinate of rotation center
            % q_0               row coordinate of rotation center
            %
            % Returns:
            % C_r               sequence of rotated channel tensors
            N_f = wd_video.getNumberOfFrames();
            C_r = cell(N_f, 1);
            for f = 1:N_f
                [~, ~, C_r{f, 1}] = wd_video ...
                    .wd_image{f, 1} ...
                    .computeRotatedIndex(phi__rad, p_0, q_0);
                wd_video ...
                    .wd_image{f, 1} ...
                    .setChannel(C_r{f, 1});
            end%for
        end%function
        function C = scaleBy(wd_video, S)
            % scaleBy    scale image sequence channel tensors
            %
            % Multiply channel tensors of the image sequence with the 
            % specified scaling factor.
            %
            % Arguments:
            % S                 scaling factor
            %
            % Returns:
            % C                 sequence of scaled channel tensors
            N_f = wd_video.getNumberOfFrames();
            C = cell(N_f, 1);
            for f = 1:N_f
                C{f, 1} = wd_video.wd_image{f, 1}.scaleBy(S);
            end%for
        end%function
        function C = scaleTo(wd_video, C_min, C_max)
            % scaleTo    scale image sequence to specified value range
            %
            % Scale pixel values of the channel tensors of the image sequence
            % to the specified value range.
            %
            % Arguments:
            % C_min             lower limits of value range
            % C_max             upper limits of value range
            %
            % Returns:
            % C                 channel tensors scaled to range
            N_f = wd_video.getNumberOfFrames();
            C = cell(N_f, 1);
            for f = 1:N_f
                C{f, 1} = wd_video ...
                    .wd_image{f, 1} ...
                    .scaleTo(C_min(f, 1), C_max(f, 1));
            end%for
        end%function
        function G = squeezeToChannel(wd_video)
            % squeezeToChannel    average to single channel image sequence
            %
            % Compute the average value of the channels of each image in this
            % sequence.
            %
            % Returns:
            % G                 single channel image sequence
            N_f = wd_video.getNumberOfFrames();
            G = cell(N_f, 1);
            for f = 1:N_f
                G{f, 1} = wd_video ...
                    .wd_image{f, 1} ...
                    .squeezeToChannel();
            end%for
        end%function
        function C = trimByPadding(wd_video, dp_0, dq_0)
            % trimByPadding    trim image sequence border by given padding
            %
            % Trim the channel tensors of the image sequence by the top and 
            % bottom padding and by the left and right padding.
            %
            % Arguments:
            % dp_0              top and bottom padding
            % dq_0              left and right padding
            %
            % Returns:
            % C                 trimmed channel tensor
            N_f = wd_video.getNumberOfFrames();
            C = cell(N_f, 1);
            for f = 1:N_f
                C{f, 1} = wd_video.wd_image{f, 1}.trimByPadding(dp_0, dq_0);
            end%for
        end%function
        function C = trimPixelValues(wd_video, C_min, C_max)
            % trimPixelValues    trim pixel values of sequence to range
            %
            % Restrict pixel values in each channel of each image frame in
            % the sequence to the minimum and maximum value given in
            % the argument list. The image sequence is updated with the trimmed
            % image frames.
            %
            % Arguments:
            % C_min     minimum pixel value
            % C_max     maximum pixel value
            %
            % Returns:
            % C         channels trimmed to interval
            N_f = wd_video.getNumberOfFrames();
            C = cell(N_f, 1);
            for f = 1:N_f
                C{f, 1} = wd_video ...
                    .wd_image{f, 1} ...
                    .trimPixelValues(C_min, C_max(f));
             end%for
        end%function
        
    %% inactive methods %%
%         function preprocess(wd_video, f_cam, t__s, c_0)
%             % preprocess    preprocess image sequence
%             N_f = length(wd_video.wd_image_1d);
%             wd_video.f_s = nan(N_f, 1);
%             wd_video.isDark = false(N_f, 1);
%             for f = 1:length(wd_video.wd_image_1d)
%                 C = wd_video.wd_image_1d{f, 1}.squeezeToChannel([]);
%                 % get frame index of illuminated and dark images
%                 i_f = 0;
%                 codeW = 3;
%                 h = double(C(1, 1 + codeW * 66));  % colum 199
%                 l = double(C(1, 1 + codeW * 67));  % colum 202
%                 wd_video.isDark(f, 1) = (C(1, 1 + codeW * 63) > (h + l) / 2);
%                 for itmp = 1:63
%                     if double(C(1, 1 + codeW * (itmp - 1))) > (h + l) / 2
%                         i_f = i_f + 1; % found set "bit"
%                     end%if
%                     i_f = 2 * i_f; % convert binary to decimal
%                 end%for
%                 wd_video.f_s(f, 1) = i_f + 1 - wd_video.isDark(f, 1);
%                 % trim pixel values to 2^10 - 1 (10bit)
%                 wd_video.N_b = 2^10;
%                 C = (wd_video.N_b - 1) * double(C) / (2^16 - 1);
%                 % trim image padding
%                 C = C( ...
%                     c_0.dp_0{1, 1}(1):(end - c_0.dp_0{1, 1}(2) - 1), ...
%                     c_0.dq_0{1, 1}(1):(end - c_0.dq_0{1, 1}(2) - 1), :);
%                 % store image, rotated to p axis along X drive
%                 [~, ~, C] = ...
%                     WDImage(C, '').computeRotatedIndex(c_0.Phi_0__rad, [], []);
%                 wd_video.wd_image_1d{f, 1} = WDImage(C, wd_video.url);
%             end%for
%             % check for non-matching frame number
%             [~, f] = ismember(wd_video.f_s, f_cam);
%             f = nonzeros(f);
%             wd_video.f_s = f;
%             wd_video.a_f = f_cam(f);
%             wd_video.t_f__s = t__s.vbh(f); % - t__s.vbh(1);
%             % compute approximate times and time increments of frame capture
%             wd_video.dt_f__s = zeros(numel(wd_video.t_f__s), 1);
%             wd_video.dtau_f = zeros(numel(wd_video.t_f__s), 1);
%             wd_video.f_red = nan(numel(wd_video.t_f__s), 1);
%             for f = 1 : numel(wd_video.t_f__s)
%                 dt = wd_video.t_f__s(f) - t__s.ebh;
%                 [~, f_min] = min(abs(dt));
%                 wd_video.f_red(f) = f_min(1);
%                 wd_video.dt_f__s(f) = dt(f_min(1));
%                 if wd_video.dt_f__s(f) < 0.0
%                     wd_video.dtau_f(f) = wd_video.dt_f__s(f) / ( ...
%                         t__s.ebh(wd_video.f_red(f)) - ...
%                         t__s.ebh(wd_video.f_red(f) - 1));
%                 else % 0.0 <= wd_video.dt_f__s(f)
%                     wd_video.dtau_f(f) = wd_video.dt_f__s(f) / ( ...
%                         t__s.ebh(wd_video.f_red(f) + 1) - ...
%                         t__s.ebh(wd_video.f_red(f)));
%                 end%if
%             end%for
%             % reference frame times to common origin
%             wd_video.t_f__s = wd_video.t_f__s - t__s.vbh(1);
%         end%function
%         function c_b = computeHistogramMoment(wd_image, N_b, k)
%             % computeHistogramMoment   compute k moment of histogram
%             %   Compute the k moment of the histogram probability
%             %   distribution.
%             c_b = zeros(wd_image.N_c);
%             [N_h, e, N_b] = wd_image.computeHistogram(N_b);
%             for c = 1 : wd_image.N_c
%                 p = N_h(c, :) / (wd_image.N_p * wd_image.N_q);
%                 c_b(c) = round(sum(p .* (1 : N_b(c)) .^ k));
%                 c_b(c) = e(c_b(c));
%             end%for
%         end%function                    
%         function threshold = computeOtsu(wd_image, N_b)
%             % computeOtsu   compute image threshold with Otsu's method
%             %   Compute the image threshold to separate pixel value into
%             %   two classes. The threshold value is computed by inter-class
%             %   variance.
%             [N_h, ~, ~] = wd_image.computeHistogram(N_b);
%             total = sum(N_h); 
% %             top = N_b;  % 256
%             sumB = 0;
%             wB = 0;  % wB = pixels in lower bins
%             maximum = 0.0;
%             sum1 = dot(0:N_b - 1, N_h);
%             for b = 1:N_b
%                 wF = total - wB;  % wF = pixels in higher bins
%                 if wB > 0 && wF > 0
%                     mF = (sum1 - sumB) / wF;
%                     val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
%                     if ( val >= maximum )
%                         threshold = b;
%                         maximum = val;
%                     end
%                 end
%                 wB = wB + N_h(b);
%                 sumB = sumB + (b - 1) * N_h(b);
%             end
% %             threshold = threshold / top;
%         end
%         function M_p = computePaddingMask(image, delta, positive)
%             % computePaddingMask   compute mask of image padding
%             % Compute a mask of the partial image with an offset (padding)
%             % from its original boundaries. The value of the padding is
%             % given by a scalar argument (same offset off all sides), by a
%             % two-element argument (offset off p and q sides), or by a
%             % four-element argument (offset off each side). The padding
%             % mask is considered to be positive if the padding or offset
%             % regions of the image are masked. Otherwise it is the negative
%             % of the padding mask.
%             
%             % compute index
%             [p, q] = image.computeIndex();
%             
%             % generate padding mask
%             if numel(delta) == 1
%                 % generate mask index
%                 M_p = (delta < p & p <= image.N_p - delta) & ...
%                     (q <= delta & image.N_q - delta < q);
%             elseif numel(delta) == 2
%                 % generate mask index
%                 M_p = (delta(1) < p & p <= image.N_p - delta(1)) & ...
%                     (delta(2) < q & q <= image.N_q - delta(2));
%             elseif numel(delta) == 4
%                 % generate mask index
%                 M_p = (delta(1) < p & p <= image.N_p - delta(2)) & ...
%                     (delta(3) < q & q <= image.N_q - delta(4));
%             else
%                 throw(MException('wdtM:SizeError', 'Wrong argument size.'));
%             end%if
%             
%             % set positive or negative of padding mask
%             M_p = (M_p == positive);
%         end%function
%         function c_P = computeQuantile(image, C, P)
%             % computeQuantile   compute P-quantile of image channels
%             %   Compute the P-quantile of the pixel values in each channel
%             %   of the given pixel fields or of this image object.
%             epsilon = 1.0e-16;
%             
%             % handle empty channel argument
%             if isempty(C)
%                 C = image.C;
%             end%if
%             P = min(1, max(0, P));
%             
%             % get size
%             if size(C, 3) == 1
%                 if P == 0.0
%                     c_P = min(C, [], 'all');
%                     return;
%                 elseif P == 1.0
%                     c_P = max(C, [], 'all') + epsilon;
%                     return;
%                 else
%                     [n_p, n_q] = size(C);
%                     C = reshape(C, [n_p * n_q, 1]);
%                     C = sort(C);
%                     c_P = 0.5 * ...
%                         (C(floor(n_p * n_q * P)) + C(ceil(n_p * n_q * P)));
%                     return;
%                 end%if
%             end%if
%             c_P = zeros(size(C, 3), 1);
%             if P == 0.0
%                 c_P = min(C, [], 1 : 2);
%             elseif P == 1.0
%                 c_P = max(C, [], 1 : 2) + epsilon;
%             else
%                 for c = 1 : size(C, 3)
%                     [n_p, n_q] = size(C(:, :, c));
%                     x = reshape(C(:, :, c), [n_p * n_q, 1]);
%                     x = sort(x);
%                     c_P = 0.5 * ( ...
%                         x(floor(n_p * n_q * P)) + x(ceil(n_p * n_q * P)));
%                 end%for
%             end%if
%         end%function
%         function [u, U] = filterByColums(image, C, delta)
%             % filterByColums   filter channels by colums
%             %   Filter pixel values in each channel by a sum filter over
%             %   the given numer of colums. Return the filtered channels and
%             %   the sum of each colum of the filtered channels.
%             
%             % init return variable (if empty, take it from channel property
%             if isempty(C)
%                 C = image.C;
%             end%if
%             u = zeros(size(C));
%             U = zeros(size(C, 3), image.N_p);
%             
%             % handle length of filter in colums
%             delta = round(delta);
%             delta = max(1, delta);
%             delta = min(round(image.N_q / 3), delta);
%             
%             % construct filter mask and filter each channel by colums
%             H = ones(1, delta);
%             if size(C, 3) == 1
%                 u = filter2(H, C, 'same');
%                 U = sum(u, 1);
%                 return;
%             end%if
%             for c = 1 : size(C, 3)
%                 u(:, :, c) = filter2(H, C(:, :, c), 'same');
%                 U(c, :) = sum(u(:, :, c), 1);
%             end%for
%         end%function
%         
% %         function filterByRows(image)
% %             
% %         end%function
%         
%         
%     % get methods %
%         function [C_ls, C_lo, C_hi] = getLevelSet(image, level)
%             % getLevelSet   get pixel values at given level
%             % Get all pixels with a value equal to the level argument. If
%             % level is out of range return NaN results. Also returns the
%             % channels below and above the given level set.
%             C_lo = nan(size(image.C));
%             C_hi = nan(size(image.C));
%             C_ls = nan(size(image.C));
%             for c = 1 : image.N_c
%                 c_max = max(image.C(:, :, c), [], 'all');
%                 c_min = min(image.C(:, :, c), [], 'all');
%                 if c_min <= level && level <= c_max
%                     C_ls(:, :, c) = (image.C(:, :, c) == level) * level;
%                     C_lo(:, :, c) = (image.C(:, :, c) < level) .* ...
%                         image.C(:, :, c) + (level <= image.C(:, :, c)) * level;
%                     C_hi(:, :, c) = (image.C(:, :, c) > level) .* ...
%                         image.C(:, :, c);
%                 end%if
%             end%for
%         end%function
%         
%         
%         
%         function [p_m, q_m, p_n, q_n] = maskCircle(image, d)
%             % maskCircle   mask pixel index by circle of given diameter            
%             % Mask the pixel index of this image with a circle of specified
%             % diameter and return the remaining image index inside the
%             % circle. Then use the mask negative and return the remaining
%             % pixel index.
%             
%             % compute index
%             [p, q] = image.computeIndex();
%             
%             % compute positive of circle mask
%             M_c = image.computeCircleMask(d, true);
%             
%             % compute masked and non-masked index
%             p_m = M_c .* p;
%             q_m = M_c .* q;
%             p_n = ~M_c .* p;
%             q_n = ~M_c .* q;            
%         end%function
%         
%         function [p_m, q_m, p_0, q_0] = maskPadding(image, delta)
%             % maskPadding   mask pixel index by padding of given width
%             % Mask the image index with a frame of specified width and
%             % return the remaining image index inside the frame. The
%             % padding width can be a scalar (applies same value to all four
%             % sides), a two-element array (applies to p and q directions),
%             % or a four-element array (applies values to all four sides).
%             
%             % compute index
%             [p, q] = image.computeIndex();
%             
%             % compute padding mask positive
%             try
%                 M_p = image.computePaddingMask(delta, true);
%             catch ME
%                 disp(['ID: ', ME.identifier]);
%                 rethrow(ME);
%             end%try
%             
%             % compute masked and non-masked index
%             p_m = M_p .* p;
%             q_m = M_p .* q;
%             p_0 = ~M_p .* p;
%             q_0 = ~M_p .* q;
%         end%function
%         
%         
%         function C_l = recurConnectedComponent(wd_image)
%             %
%             %
%             
%             C_l = double(~wd_image.C);
%             label = 0;
%             findComponent(C_l, label);
%             
%             function findComponent(C_l, label)
%                 for p = 1 : wd_image.N_p
%                     for q = 1 : wd_image.N_q
%                         if C_l(p, q) == -1
%                             label = label + 1;
%                             searchLabel(C_l, label, p, q);
%                         end%if
%                     end%for
%                 end%for                            
%             end%function
%             
%             function searchLabel(C_l, label, p, q)
%                 C_l(p, q) = label;
%                 
%             end%function
%         end%function
%         
%         function c = computeGammaCorrection()
%         
%         function c = computeSigmoidCorrection(image, pixel_field, sigma)
%                         
%             % define function of sigmoid correction
%             h = @(x, s) 1.0 ./ (1.0 + exp(-s * (x - 0.5)));
%             
%             % init return variable (if empty, take it from channel property
%             if nargin == 0 || isempty(pixel_field)
%                 pixel_field = image.C;
%             end%if
%             c = zeros(size(pixel_field));
%             
%             % get size and compute for single channel
%             s = size(pixel_field);
%             if length(s) == 2
%                 c = h(pixel_field, sigma);
%                 
%             % compute for multi-channel
%             elseif 2 < length(s)
%                 for ch = 1 : s(3)
%                     c = h(pixel_field(:, :, ch), sigma);
%                 end%for
%                 
%             % handle unsupported cases
%             else
%                 c = [];
%             end%if                
%         end%function
%         
%         
%         function [grad, phi] = transformGradientToPolar(image, grad_p, grad_q)
%             % TRANSFORM_GRADIENT_TO_POLAR Gradient in polar coordinates
%             %   Transforms the components of the gradient from given
%             %   cartesian to polar coordinates. However, the Euclidian norm
%             %   is replaced by the mean squared values.
%             % Arguments:
%             %   * grad_p       cartesian rows gradient
%             %   * grad_q       cartesian colums gradient
%             % Returns:
%             %   * grad         norm approximated by mean squared values
%             %   * phi          polar angle [rad]
%             
%             
%             % define return variables
%             grad_p = squeeze(zeros(image.N_p, image.N_q, image.N_C));
%             grad_q = squeeze(zeros(image.N_p, image.N_q, image.N_C));
%            
%             % compute for single channel
%             if image.N_C == 1
%                 grad = 0.5 * (grad_p .^ 2 + grad_q .^ 2);
%                 phi = atan2(grad_q, grad_p);
%                 phi = unwrap(phi);
%                
%             % compute for multi-channel
%             elseif 1 < image.N_C
%                 for c = 1 : image.N_C
%                     grad(:, :, c) = ...
%                         0.5 * (grad_p(:, :, c) .^ 2 + grad_q(:, :, c) .^ 2);
%                     phi(:, :, c) = atan2(grad_q(:, :, c), grad_p(:, :, c));
%                     phi(:, :, c) = unwrap(phi(:, :, c));
%                 end%for
%            
%             % handle the unsupported case
%             else
%                 grad = [];
%                 phi = [];
%             end%if
%         end%function   

%function threshold = otsu_2D(hists, total)
%maximum = 0.0;
%threshold = 0;
%helperVec = 0:255;
%mu_t0 = sum(sum(repmat(helperVec',1,256).*hists));
%mu_t1 = sum(sum(repmat(helperVec,256,1).*hists));
%p_0 = zeros(256);
%mu_i = p_0;
%mu_j = p_0;
%for ii = 1:256
%    for jj = 1:256
%        if jj == 1
%            if ii == 1
%                p_0(1,1) = hists(1,1);
%            else
%                p_0(ii,1) = p_0(ii-1,1) + hists(ii,1);
%                mu_i(ii,1) = mu_i(ii-1,1)+(ii-1)*hists(ii,1);
%                mu_j(ii,1) = mu_j(ii-1,1);
%            end
%        else
%            p_0(ii,jj) = p_0(ii,jj-1)+p_0(ii-1,jj)-p_0(ii-1,jj-1)+hists(ii,jj);
%            mu_i(ii,jj) = mu_i(ii,jj-1)+mu_i(ii-1,jj)-mu_i(ii-1,jj-1)+
%(ii-1)*hists(ii,jj);
%            mu_j(ii,jj) = mu_j(ii,jj-1)+mu_j(ii-1,jj)-mu_j(ii-1,jj-1)+
%(jj-1)*hists(ii,jj);
%        end
%
%        if (p_0(ii,jj) == 0)
%            continue;
%        end
%        if (p_0(ii,jj) == total)
%            break;
%        end
%        tr = ((mu_i(ii,jj)-p_0(ii,jj)*mu_t0)^2 + (mu_j(ii,jj)-%p_0(ii,jj)*mu_t1)^2)/(p_0(ii,jj)*(1-p_0(ii,jj)));
%
%        if ( tr >= maximum )
%            threshold = ii;
%            maximum = tr;
%        end
%    end
%end
%end     
    end%methods
    
end%classdef


% end of module wdM.wdtM.WDVideo

classdef WDImage < handle
    % WDImage    class of digital image objects
    %
    % The class WDImage represents common digital pictures and provides the
    % methods to process the channels of pixel values.
    %
    % 
    % Properties (private)
    % --------------------
    % C                         2D (gray value) or 3D (color) channel tensor
    % name                      name of image or image file
    %
    %
    % Methods (static)
    % ----------------
    % createCircle              create binary WDImage object with circle
    % loadImage                 load image file and construct WDImage object
    %
    %
    % Methods (public)
    % ----------------
    % WDImage                   construct this WDImage object
    %
    % getChannel                get channel tensor
    % getHeight                 get number of pixel rows
    % getSize                   get number of pixel rows, colums and channels
    % getWidth                  get number of pixel colums
    % setChannel                set channel tensor
    %
    % computeDistribution       compute counts of pixel values in bins
    % computeGradient           compute gradient along rows and colums
    % computeIndex              compute index grid according to image size
    % computeRotatedIndex       compute index rotated about given angle
    % computeThreshold          compute image threshold
    % generateBinary            generate binary image from threshold
    % normalizeToUnity          normalize pixel values into unit interval
    % scaleBy                   scale pixel values by given factor
    % scaleTo                   scale pixel values to specified value range
    % squeezeToChannel          average all channels to single channel
    % trimByPadding             trim image border by given padding
    % trimPixelValues           trim pixel values to range
    
    
    % mask                      apply given mask to image channels
    
    
    %   computeCircleMask       compute circular mask with given diameter
    %   computeGradient         compute gradient along rows and colums
    %   computeHistogramMoment  compute k moment of histogram
    %   computeOtsu             compute image threshold with Otsu's method
    %   computePaddingMask      compute mask of image padding
    %   computeQuantile         compute P-quantile of image channels
    %   filterByColums          filter channels by colums
    %   getLevelSet             get pixel values at given level
    %   maskPadding             mask pixel index by padding of given width           
    % ________________________________________________________________________
    %
    % to do:
    % 1. implement filterByRows() method
    % 2. ...
    %
    % author:   Stefan Wittwer, info@wittwer-datatools.ch
    % _________________________________________________________________________
    
    
%% properties (private) %%
    properties (Access = private)
        C               % 2D or 3D tensor of pixel values
        name    string  % name of image
    end%properties
    
    
%% methods (static) %%
    methods (Static = true)
        function wd_image = createCircle(M, N, isFilled)
            % createCircle    create binary WDImage object with circle
            %
            % Create a WDImage object with a binary pixel field showing a
            % centered, filled or unfilled circle of diameter min(M, N).
            %
            % Arguments:
            % M         number of pixel rows >= 3
            % N         number of pixel colums >= 3
            % isFilled  fills circle if true, draws circle border otherwise
            % 
            % Returns:
            % wd_image  binary WDImage object with incircle
            arguments
                M        (1, 1) int
                N        (1, 1) int
                isFilled (1, 1) bool = true
            end%arguments
            M = max(3, M);
            N = max(3, N);
            C = zeros(M, N);
            % compute pixel coordinates of center and radius
            p_0 = M / 2;
            q_0 = N / 2;
            r = min(M, N) / 2;
            % create binary image
            p = 1 : M;
            q = 1 : N;
            if isFilled
                C(((p - p_0) .^ 2 + (q - q_0) .^ 2) <= r ^ 2) = 1;
                wd_image = WDImage(C, "CircleFilled");
            else
                C(((p - p_0) .^ 2 + (q - q_0) .^ 2) == r ^ 2) = 1;
                wd_image = WDImage(C, "CircleEmpty");
            end%if            
        end
        function wd_image = loadImage(wdio_image)
            % loadImage    load image file and construct WDImage object
            %
            % Load an image file via a WDIOImage interface object. A new 
            % WDImage object is constructed from the loaded channels with 
            % its pixel values.
            %
            % Arguments:
            % wdio_image        WDIOImage data I/O object
            %
            % Returns:
            % wd_image          WDImage object with loaded image data
            C = wdio_image.clone();
            name = wdio_image.getName();
            wd_image = WDImage(C, name);
        end%function        
    end%methods
    
    
%% methods (public) %%
    methods
    %% constructor %%
        function wd_image = WDImage(C, name)
            % WDImage    construct this WDImage object
            %
            % Initialize the channels for colours and opacity, and set the
            % name property of this image object.
            %
            % Arguments:
            % C     pixel value field as a tensor array object
            % name  string with image name
            %
            % Returns:
            % wd_image  this WDImage object
            wd_image.C = double(C);
            wd_image.name = name;
        end%function
        
    %% get and set methods %%
        function C = getChannel(wd_image, varargin) 
            % getChannel    get channel tensor
            if nargin == 1
                C = wd_image.C;
            else
                c = min(1, max(size(wd_image.C, 3), varargin{1, 1}));
                C = wd_image.C(:, :, c);
            end%if
            C = squeeze(C);
        end%function
        function N_p = getHeight(wd_image)
            % getHeight    get number of pixel rows
            N_p = size(wd_image.C, 1);
        end%function
        function [N_p, N_q, N_c] = getSize(wd_image)
            % getSize    get number of rows, colums and channels
            N_p = size(wd_image.C, 1); % number of rows = image height
            N_q = size(wd_image.C, 2); % number of cols = image width
            N_c = size(wd_image.C, 3);
        end%function
        function N_q = getWidth(wd_image)
            % getWidth    get number of pixel colums
            N_q = size(wd_image.C, 2);
        end%function
        function C = setChannel(wd_image, C)
            % setChannel    set channel tensor
            C = double(C);
            wd_image.C = C;
        end%function
        
    %% computational methods %%
        function [N_h, b, N_b] = computeDistribution(wd_image, N_b)
            % computeDistribution    compute counts of pixel values in bins
            %
            % Compute the counts or frequency of pixel values in a series
            % of value ranges (bins). All pixel values of this image are
            % considered.
            %
            % Arguments:
            % N_b               number of bins >= 3
            %
            % Returns:
            % N_h               counts per channel
            % b                 bins edges
            % N_b               number of bins
            N_b = max(3, N_b);
            N_h = zeros(size(wd_image.C, 3), N_b); % counts per channel
            b = zeros(size(wd_image.C, 3), N_b + 1); % bins edges
            for c = 1:size(wd_image.C, 3)
                [N_h(c, :), b(c, :)] = histcounts(wd_image.C(:, :, c), N_b);
            end%for
        end%function
        function [grad_p, grad_q, M] = computeGradient(wd_image, C, method)
            % computeGradient    compute gradient along rows and colums
            %
            % Compute gradient of pixel value fields in each channel of the
            % given channel tensor. If the channel tensor argument is empty
            % the gradients are computed on the channel tensor of this image
            % object. The gradients are computed along pixel rows and
            % colums.
            %
            % Arguments:
            % C                 user defined channels or empty
            % method            name of selected gradient method or empty
            %
            % Returns:
            % grad_p            pixel value gradient along rows
            % grad_q            pixel value gradient along colums
            % M                 gradient operator of selected method
            if isempty(C)
                C = wd_image.C;
            end%if
            grad_p = zeros(size(C));
            grad_q = zeros(size(C));
            % set scheme according selected method
            switch (method)
                case 'sobel'
                    M.p = [ ...
                        [-1.0, -2.0, -1.0]; ...
                        [ 0.0,  0.0,  0.0]; ...
                        [ 1.0,  2.0,  1.0]] / 6.0;
                    M.q = [ ...
                        [-1.0,  0.0,  1.0]; ...
                        [-2.0,  0.0,  2.0]; ...
                        [-1.0,  0.0,  1.0]] / 6.0;
                otherwise % default method is prewitt
                    M.p = [ ...
                        [-1.0, -1.0, -1.0]; ...
                        [ 0.0,  0.0,  0.0]; ...
                        [ 1.0,  1.0,  1.0]] / 6.0;
                    M.q = [ ...
                        [-1.0,  0.0,  1.0]; ...
                        [-1.0,  0.0,  1.0]; ...
                        [-1.0,  0.0,  1.0]];
            end%switch
            % compute gradients
            for c = 1 : size(C, 3)
                grad_p(:, :, c) = filter2(M.p, C(:, :, c)); % along row
                grad_q(:, :, c) = filter2(M.q, C(:, :, c)); % along col
            end%for
        end%function
        function [p, q] = computeIndex(wd_image)
            % computeIndex    compute index grid according to image size
            %
            % Compute the pixel index of this image along the pixel rows and
            % colums. Then create an index meshgrid.
            %
            % Returns:
            % p         index along p axis
            % q         index along q axis
            p = 1:size(wd_image.C, 1);
            q = 1:size(wd_image.C, 2);
            [q, p] = meshgrid(q, p);            
        end%function
        function [p_r, q_r, C_r] = computeRotatedIndex(wd_image, ...
                phi__rad, p_0, q_0)
            % computeRotatedIndex    compute index rotated about given angle
            %
            % Rotate the index of the image pixels by the negative of the 
            % angle argument. The rotated index is then applied to the 
            % channels and returns a pixel field rotated about the negative 
            % angle argument. If the center index arguments are empty, the 
            % center of this image is taken.
            %
            % Arguments:
            % phi__rad  rotation angle in [rad]
            % p_0       center pixel row
            % q_0       center pixel colum
            %
            % Returns:
            % p_r       rotated pixel row coordinates
            % q_r       rotated pixel colum coordinates
            % C_r       rotated channels
            C_r = zeros(size(wd_image.C(:, :, 1)));
            % handle angle and center index argument
            if isnan(phi__rad) || isempty(phi__rad) || ~isscalar(phi__rad)
                phi__rad = 0.0;
            end%if
            if isempty(p_0)
                p_0 = size(wd_image.C, 1) / 2.0;
            end%if
            if isempty(q_0)
                q_0 = size(wd_image.C, 2) / 2.0;
            end%if
            % compute image index and return it for zero rotation
            [p, q] = wd_image.computeIndex();
            if mod(phi__rad, 2.0 * pi) == 0.0
                p_r = p;
                q_r = q;
                C_r = wd_image.C;
                return;
            end%if
            % compute rotated index
            p_r = round( ...
                (p - p_0) * cos(phi__rad) - (q - q_0) * sin(phi__rad) + p_0);
            q_r = round( ...
                (p - p_0) * sin(phi__rad) + (q - q_0) * cos(phi__rad) + q_0);
            % make sure that rotated index is valid
            if ~isempty(p_r(p_r == 0))
                p_r = p_r + 1;
            end%if
            if ~isempty(q_r(q_r == 0))
                q_r = q_r + 1;
            end%if            
            p_r(isnan(p_r)) = 0;
            q_r(isnan(q_r)) = 0;
            p_r = (0 < p_r & p_r <= size(wd_image.C, 1)) .* p_r;
            q_r = (0 < q_r & q_r <= size(wd_image.C, 2)) .* q_r;
            % compute rotated pixel values
            for i = 1 : size(wd_image.C, 1)
                for j = 1 : size(wd_image.C, 2)
                    p_r(i, j) = round(p_r(i, j));
                    q_r(i, j) = round(q_r(i, j));
                    if p_r(i, j) && q_r(i, j)
                        C_r(i, j) = wd_image.C(p_r(i, j), q_r(i, j));
                    end%if
                end%for
            end%for
        end%function
        function C_t = computeThreshold(wd_image, method, P_h, N_b)
            % computeThreshold    compute image threshold
            %
            % Compute the image threshold of each channel of the channel 
            % tensor. The given method is applied.
            %
            % Arguments:
            % method            method applied to compute image threshold
            % P_h               probability distribution of pixel values
            % N_b               number of values per pixel
            %
            % Returns:
            % C_t               image threshold
            switch (method)
                case '0<C'
                    C_t = zeros(1, size(wd_image.C, 3));
                case '1<C'
                    C_t = ones(1, size(wd_image.C, 3));
                case '2^8-1<C'
                    C_t = ones(1, size(wd_image.C, 3)) * (2^8 - 1);
                case 'N_b-2<C'
                    C_t = ones(1, size(wd_image.C, 3)) * (N_b - 2);
                case '1.moment'
                    C_t = sum((0:N_b - 1) .* P_h, 2);
                case '2.moment'
                    C_t = sqrt(sum((0:N_b - 1) .^ 2.0 .* P_h, 2));
                case 'otsu'
                    C_t = wd_image.computeOtsu(N_b);
                case 'C(P=5%)'
                    C_t = find(0.05 < cumsum(P_h, 2), 1, 'first');
                case 'C(P=10%)'
                    C_t = find(0.10 < cumsum(P_h, 2), 1, 'first');
                case 'C(P=25%)'
                    C_t = find(0.25 < cumsum(P_h, 2), 1, 'first');
                case 'C(P=50%)'
                    C_t = find(0.50 < cumsum(P_h, 2), 1, 'first');
                case 'C(P=75%)'
                    C_t = find(0.75 < cumsum(P_h, 2), 1, 'first');
                case 'F<=1/C'
                    c_0 = find( ...
                        (1.0 ./ double(2:size(P_h, 2))) < P_h(:, 2:end), ...
                        1, 'first');
                    C_t = ~isempty(c_0) .* c_0;
                otherwise
                    msg = [ ...
                        'Image threshold method ' ...
                        method ...
                        ' not supported.'];
                    throw( ...
                        MException('WDImage:UnsupportedImageThreshold', msg));
            end%switch
        end%function
        function binary = generateBinary(wd_image, C_0, isaPositive)
            % generateBinary    generate binary image from threshold
            %
            % Generate binary image using the given image threshold: pixel
            % values below or equal threshold value become 0, all other
            % pixel values are 1.
            %
            % Arguments:
            % C_0               image threshold
            % isaPositive       if true generate positive binary image
            B = ((C_0 < mean(wd_image.C, 3)) == isaPositive);
            binary = WDImage(B, "binary");
        end%function
        function u = normalizeToUnity(wd_image, varargin)
            % normalizeToUnity    normalize pixel values into unit interval
            %
            % Normalize the pixel values into the unit interval [0; 1].
            % Normalization is done over all channels of the channel tensor. 
            % If the channel tensor argument is not empty the method computes
            % the unit norm of the given channel tensor. Otherwise it
            % normalizes the pixel field of this image.
            %
            % Arguments:
            % varargin          input channel tensor or empty
            %
            % Returns:
            % u                 unit norm channel
            if isempty(varargin)
                u = wd_image.getChannel();
            else
                u = varargin;
            end%if
            % compute norm of channels
            C_max = max(u, [], 'all');
            C_min = min(u, [], 'all');
            u = (u - C_min) / (C_max - C_min);
            if isempty(varargin)
                wd_image.C = u;
            end%if
        end%function
        function C = scaleBy(wd_image, S, varargin)
            % scaleBy    scale pixel values by given factor
            %
            % Multiply channel tensor with the specified scaling factor. The 
            % variable argument is empty or contains an external channel 
            % tensor. If it is empty scale the channel tensor of this image
            % object, otherwise scale the provided channel tensor.
            %
            % Arguments:
            % S                 scaling factor
            % varargin          external channel tensor or empty
            %
            % Returns:
            % C                 scaled channel tensor
            if isempty(varargin)
                C = wd_image.C * S;
                wd_image.C = C;
            else
                C = varargin{1, 1} * S;
            end%if
        end%function
        function C = scaleTo(wd_image, C_min, C_max, varargin)
            % scaleTo   scale pixel values to specified value range
            %
            % Scale pixel values of the channel tensor into the specified 
            % range. The variable argument is either empty or contains an
            % external channel tensor to be scaled.
            %
            % Arguments:
            % C_min             lower limit of value range
            % C_max             upper limit of value range
            % varargin          external channel tensor or empty
            %
            % Returns:
            % C                 channel tensor scaled to range
            if C_max <= C_min
                c = c_min;
                C_max = C_min;
                C_min = c;
            end%if
            % scale
            u = wd_image.normalizeToUnity(varargin);
            C = C_min + (C_max - C_min) .* u;
            if isempty(varargin)
                wd_image.setChannel(C);
            end%if
        end%function
        function G = squeezeToChannel(wd_image, varargin)
            % squeezeToChannel    average all channels to single channel
            %
            % Compute the mean value along the channels of this image
            % (default), or along the channels of the given non-empty
            % argument.
            %
            % Arguments:
            % varargin          input channels from another image or empty
            %
            % Returns:
            % G                 pixel field with averaged gray values
            if isempty(varargin)
                N_c = size(wd_image.C, 3);
                if 1 == N_c
                    G = squeeze(wd_image.C);
                elseif 1 < N_c
                    G = squeeze(mean(wd_image.C, 3));
                end%if
                wd_image.setChannel(G);
            else
                G = squeeze(mean(varargin{1, 1}, 3));
            end%if
        end%function
        function C = trimByPadding(wd_image, dp_0, dq_0, varargin) 
            % trimByPadding    trim image border by given padding
            %
            % Trim the channel tensor by the top and bottom padding and by 
            % the left and right padding. The channel tensor of this image
            % object is trimmed if the variable argument is empty, otherwise
            % the external channel tensor is trimmed.
            %
            % Arguments:
            % dp_0              top and bottom padding
            % dq_0              left and right padding
            % varargin          external channel tensor or empty
            %
            % Returns:
            % C                 trimmed channel tensor
            if isempty(varargin)
                C = wd_image.C( ...
                    dp_0(1):(end - dp_0(2) - 1), ...
                    dq_0(1):(end - dq_0(2) - 1), :);
                wd_image.C = C;
            else
                C = varargin{1, 1}( ...
                    dp_0(1):(end - dp_0(2) - 1), ...
                    dq_0(1):(end - dq_0(2) - 1), :);
            end%if
        end%function

        function C = trimPixelValues(wd_image, C_min, C_max)
            % trimPixelValues   trim pixel values to range
            % 
            % Trim pixel values in each channel to the minimum and
            % maximum limits given in the argument list. After trimming, pixel
            % values below the minimum are equal to zero, and pixel values
            % above the maximum are equal to the maximum.
            %
            % Arguments:
            % C_min     minimum pixel value
            % C_max     maximum pixel value
            %
            % Returns:
            % C         channels trimmed to range
            C_max = max(C_min, C_max);
            C_min = min(C_min, C_max);
            % trim pixel values to interval [C_min; C_max]
            wd_image.C(wd_image.C < C_min) = C_min;
            wd_image.C(C_max < wd_image.C) = C_max;
            C = wd_image.C;
        end%function
        
    %% masking methods %%
       function  mask(wd_image, M)
            % mask    apply given mask to image channels
            [N_p, N_q, N_c] = wd_image.getSize();
            N_M = size(M);
            if all([N_p N_q] == N_M)
                for c = 1:N_c
                    wd_image.C(:, :, c) = wd_image.C(:, :, c) .* M;
                end%for
            else
                throw(MException("WDImage:DimensionError", ...
                    "Size of mask must match image size."));
            end%if
        end%function
        function [M_c, p_0, q_0, d] = computeCircleMask(image, ...
                d, p_0, q_0, positive)
            % computeCircleMask   compute circular mask with given diameter
            % Compute a pixel mask with the shape of a circle. Its diameter
            % is the minimum of the width and the height by default, but it
            % can be specified by the argument d. The circle mask is
            % considered to be positive if the pixels outside the circle
            % are masked (false). Otherwise it is the negative of the mask.
            if isempty(d)
                d = min(image.N_p, image.N_q);
            end%if
            if isempty(p_0) || isempty(q_0)
                p_0 = 0.5 * image.N_p;
                q_0 = 0.5 * image.N_q;
            end%if
            [p, q] = image.computeIndex();            
            M_c = (p - p_0) .^ 2.0 + (q - q_0) .^ 2.0 < 0.25 * d ^ 2.0;
            M_c = (M_c == positive);
        end%function
        
        
        function c_b = computeHistogramMoment(wd_image, N_b, k)
            % computeHistogramMoment   compute k moment of histogram
            %   Compute the k moment of the histogram probability
            %   distribution.
            c_b = zeros(wd_image.N_c);
            [N_h, e, N_b] = wd_image.computeHistogram(N_b);
            for c = 1 : wd_image.N_c
                p = N_h(c, :) / (wd_image.N_p * wd_image.N_q);
                c_b(c) = round(sum(p .* (1 : N_b(c)) .^ k));
                c_b(c) = e(c_b(c));
            end%for
        end%function            
        
        function threshold = computeOtsu(wd_image, N_b)
            % computeOtsu   compute image threshold with Otsu's method
            %   Compute the image threshold to separate pixel value into
            %   two classes. The threshold value is computed by inter-class
            %   variance.
            [N_h, ~, ~] = wd_image.computeHistogram(N_b);
            total = sum(N_h); 
%             top = N_b;  % 256
            sumB = 0;
            wB = 0;  % wB = pixels in lower bins
            maximum = 0.0;
            sum1 = dot(0:N_b - 1, N_h);
            for b = 1:N_b
                wF = total - wB;  % wF = pixels in higher bins
                if wB > 0 && wF > 0
                    mF = (sum1 - sumB) / wF;
                    val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
                    if ( val >= maximum )
                        threshold = b;
                        maximum = val;
                    end
                end
                wB = wB + N_h(b);
                sumB = sumB + (b - 1) * N_h(b);
            end
%             threshold = threshold / top;
        end

        function M_p = computePaddingMask(wd_image, delta, positive)
            % computePaddingMask   compute mask of image padding
            % Compute a mask of the partial image with an offset (padding)
            % from its original boundaries. The value of the padding is
            % given by a scalar argument (same offset off all sides), by a
            % two-element argument (offset off p and q sides), or by a
            % four-element argument (offset off each side). The padding
            % mask is considered to be positive if the padding or offset
            % regions of the image are masked. Otherwise it is the negative
            % of the padding mask.
            
            % compute size and index
            [N_p, N_q, ~] = wd_image.getSize();
            [p, q] = wd_image.computeIndex();
            
            % generate padding mask
            if numel(delta) == 1
                % generate mask index
                M_p = (delta < p & p <= N_p - delta) & ...
                    (q <= delta & N_q - delta < q);
            elseif numel(delta) == 2
                % generate mask index
                M_p = (delta(1) < p & p <= N_p - delta(1)) & ...
                    (delta(2) < q & q <= N_q - delta(2));
            elseif numel(delta) == 4
                % generate mask index
                M_p = (delta(1) < p & p <= N_p - delta(2)) & ...
                    (delta(3) < q & q <= N_q - delta(4));
            else
                throw(MException('wdtM:SizeError', 'Wrong argument size.'));
            end%if
            
            % set positive or negative of padding mask
            M_p = (M_p == positive);
        end%function
        
        function c_P = computeQuantile(image, C, P)
            % computeQuantile   compute P-quantile of image channels
            %   Compute the P-quantile of the pixel values in each channel
            %   of the given pixel fields or of this image object.
            epsilon = 1.0e-16;
            
            % handle empty channel argument
            if isempty(C)
                C = image.C;
            end%if
            P = min(1, max(0, P));
            
            % get size
            if size(C, 3) == 1
                if P == 0.0
                    c_P = min(C, [], 'all');
                    return;
                elseif P == 1.0
                    c_P = max(C, [], 'all') + epsilon;
                    return;
                else
                    [n_p, n_q] = size(C);
                    C = reshape(C, [n_p * n_q, 1]);
                    C = sort(C);
                    c_P = 0.5 * ...
                        (C(floor(n_p * n_q * P)) + C(ceil(n_p * n_q * P)));
                    return;
                end%if
            end%if
            c_P = zeros(size(C, 3), 1);
            if P == 0.0
                c_P = min(C, [], 1 : 2);
            elseif P == 1.0
                c_P = max(C, [], 1 : 2) + epsilon;
            else
                for c = 1 : size(C, 3)
                    [n_p, n_q] = size(C(:, :, c));
                    x = reshape(C(:, :, c), [n_p * n_q, 1]);
                    x = sort(x);
                    c_P = 0.5 * ( ...
                        x(floor(n_p * n_q * P)) + x(ceil(n_p * n_q * P)));
                end%for
            end%if
        end%function
        
        
        function [u, U] = filterByColums(image, C, delta)
            % filterByColums   filter channels by colums
            %   Filter pixel values in each channel by a sum filter over
            %   the given numer of colums. Return the filtered channels and
            %   the sum of each colum of the filtered channels.
            
            % init return variable (if empty, take it from channel property
            if isempty(C)
                C = image.C;
            end%if
            u = zeros(size(C));
            U = zeros(size(C, 3), image.N_p);
            
            % handle length of filter in colums
            delta = round(delta);
            delta = max(1, delta);
            delta = min(round(image.N_q / 3), delta);
            
            % construct filter mask and filter each channel by colums
            H = ones(1, delta);
            if size(C, 3) == 1
                u = filter2(H, C, 'same');
                U = sum(u, 1);
                return;
            end%if
            for c = 1 : size(C, 3)
                u(:, :, c) = filter2(H, C(:, :, c), 'same');
                U(c, :) = sum(u(:, :, c), 1);
            end%for
        end%function
        
%         function filterByRows(image)
%             
%         end%function
        
        function gray = generateGrayScale(wd_image)
            %
            if 1 < wd_image.N_c
                gray = WDImage(mean(wd_image.C, 3), wd_image.name);
            else
                gray = wd_image;
            end%if
        end%function
        
        function [C_ls, C_lo, C_hi] = getLevelSet(image, level)
            % getLevelSet   get pixel values at given level
            % Get all pixels with a value equal to the level argument. If
            % level is out of range return NaN results. Also returns the
            % channels below and above the given level set.
            C_lo = nan(size(image.C));
            C_hi = nan(size(image.C));
            C_ls = nan(size(image.C));
            for c = 1 : image.N_c
                c_max = max(image.C(:, :, c), [], 'all');
                c_min = min(image.C(:, :, c), [], 'all');
                if c_min <= level && level <= c_max
                    C_ls(:, :, c) = (image.C(:, :, c) == level) * level;
                    C_lo(:, :, c) = (image.C(:, :, c) < level) .* ...
                        image.C(:, :, c) + (level <= image.C(:, :, c)) * level;
                    C_hi(:, :, c) = (image.C(:, :, c) > level) .* ...
                        image.C(:, :, c);
                end%if
            end%for
        end%function
        
        function [p_m, q_m, p_n, q_n] = maskCircle(image, d)
            % maskCircle   mask pixel index by circle of given diameter            
            % Mask the pixel index of this image with a circle of specified
            % diameter and return the remaining image index inside the
            % circle. Then use the mask negative and return the remaining
            % pixel index.
            
            % compute index
            [p, q] = image.computeIndex();
            
            % compute positive of circle mask
            M_c = image.computeCircleMask(d, true);
            
            % compute masked and non-masked index
            p_m = M_c .* p;
            q_m = M_c .* q;
            p_n = ~M_c .* p;
            q_n = ~M_c .* q;            
        end%function
        
        function [p_m, q_m, p_0, q_0] = maskPadding(image, delta)
            % maskPadding   mask pixel index by padding of given width
            % Mask the image index with a frame of specified width and
            % return the remaining image index inside the frame. The
            % padding width can be a scalar (applies same value to all four
            % sides), a two-element array (applies to p and q directions),
            % or a four-element array (applies values to all four sides).
            
            % compute index
            [p, q] = image.computeIndex();
            
            % compute padding mask positive
            try
                M_p = image.computePaddingMask(delta, true);
            catch ME
                disp(['ID: ', ME.identifier]);
                rethrow(ME);
            end%try
            
            % compute masked and non-masked index
            p_m = M_p .* p;
            q_m = M_p .* q;
            p_0 = ~M_p .* p;
            q_0 = ~M_p .* q;
        end%function
        
        function C_l = recurConnectedComponent(wd_image)
            %
            %
            
            C_l = double(~wd_image.C);
            label = 0;
            findComponent(C_l, label);
            
            function findComponent(C_l, label)
                for p = 1 : wd_image.N_p
                    for q = 1 : wd_image.N_q
                        if C_l(p, q) == -1
                            label = label + 1;
                            searchLabel(C_l, label, p, q);
                        end%if
                    end%for
                end%for                            
            end%function
            
            function searchLabel(C_l, label, p, q)
                C_l(p, q) = label;
                
            end%function
        end%function
            
        
    %% inactive methods %%
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


% end of module wdM.wdtM.WDImage

classdef WDIOImage < WDataIO
    % WDIOImage    Data I/O Interface for Image Data
    %
    % The WDIOImage class is subclassed from WDataIO and implements a data
    % I/O interface to image files. The list of the supported image 
    % formats is maintained inside the constructor.
    %
    % The current supported formats are:
    %   Image:
    %   *.bmp
    %   *.jpg
    %   *.png
    %    
    % The pixel fields of the loaded image is used as input of the 
    % wdM.wdtM.WDImage constructor.
    %
    %
    % Properties (constant)
    % ---------------------
    % IMAGE_EXTENSION   supported image file formats
    %
    %
    % Methods
    % -------
    % WDIOImage         construct instance of this WDIOImage class
    % clone             load image data from image file
    % commit            save image data to image file
    % connect           connect by checking existence of image file
    % disconnect        disconnect by resetting info property
    % getInfo           get info about image data
    %
    %
    % ________________________________________________________________________
    %
    % todo:
    % 1. ...
    % 
    % author:   Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    properties (Constant = true)
        IMAGE_FORMAT = {'.bmp', '.jpg', '.png'};
    end%properties
        
    methods
    % constructor %
        function wdio_image = WDIOImage(url)
            % WDIOImage    construct instance of this WDIOImage class
            %
            % Construct an object of this WDIOImage class to implement a
            % data I/O to an image file. The name of the image file is 
            % specified in the url argument. The url argument must contain a 
            % supported image file extension.
            %
            % Arguments:
            %   url             file path to image file
            %
            % Returns:
            % wdio_image        instance of this WDIOImage interface
            %
            % Throws:
            % ArgumentError     empty or unsupported url argument
            if isempty(url)
                throw(MException( ...
                    'WDIOImage:ArgumentError', 'Empty filename.'));
            end%if
            % check url match with supported format
            c = strfind(url, '.');
            is_image = any(strcmp(url(c(end):end), WDIOImage.IMAGE_FORMAT));
            if ~is_image
                throw(MException('WDIOImage:ArgumentError', ...
                    'Unsupported file type.'));
            end%if
            wdio_image = wdio_image@WDataIO(0, url);
        end%function
        
    % data I/O methods %
        function C = clone(wdio_image)
            % clone    load image or video data
            %
            % Load and return raw (unprocessed) image data from from the file
            % specified by the url property.
            %
            % Returns:
            % C     image color channels as 2D or 3D pixel field
            C = imread(wdio_image.url);
        end
        function commit(~, data, target)
            % commit    write image data to image file
            % 
            % Write image data in the data argument to the image file 
            % specified by the target argument. The value of the target
            % argument may differ from the url, but should be a valid image 
            % file name of supported type.
            %
            % Arguments:
            % data      2D or 3d pixel field of image
            % target    url of target image file
            %
            % Throws:
            % ArgumentError     empty or unsupported target image filename
            if isempty(target) 
                throw(MException('WDIOImage:ArgumentError', ...
                    'Empty filename.'));
            end%if
            % check target url match with supported format
            c = strfind(target, '.');
            is_image = any( ...
                strcmp(target(c(end):end), WDIOImage.IMAGE_FORMAT));
            % eventually commit
            if is_image
                imwrite(data, target, image_extension);
            else
                throw(MException('WDIOImage:ArgumentError', ...
                    'Unsupported target image file type.'));
            end%if
        end%function
        
    % connection methods %
        function is_connected = connect(wdio_image)
            % connect    connect by checking existence of image file
            %
            % Connect to the current image file specified in the url
            % property by checking its existence. If the file could not be
            % opened the method returns false.
            is_connected = isfile(wdio_image.url);
        end%function
        function is_connected = disconnect(wdio_image)
            % disconnect    disconnect by resetting info property
            % Disconnect by clearing the info property. Set the return
            % value to false.
            wdio_image.info = WDFieldSet('', []);
            is_connected = false;
        end%function
        
    % get methods %
        function info = getInfo(wdio_image)
            % getInfo    get info about image data            
            %
            % Get info about the image file and image properties.
            %
            % Returns: 
            % info     structure with info data
            info = imfinfo(wdio_image.url);
        end%function
    end%methods
end%classdef


% end of module wdM.wdbM.WDIOImage

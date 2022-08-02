classdef WDIOVideo < WDataIO
    % WDIOVideo  Data I/O Interface for Video File Data
    %
    % The WDIOVideo class is a subclass of WDataIO and implements a data I/O 
    % interface to video files. The list of the video formats currently 
    % supported is maintained in a static property belos.
    %    
    % The pixel fields of the loaded image series are used to construct an 
    % array of WDImage objects.
    %
    %
    % Constant Properties
    % -------------------
    % VIDEO_EXTENSION   list of supported video formats
    %
    %
    % Methods
    % -------
    % WDIOVideo         construct instance of this WDIOVideo class
    % clone             load video data from file
    % commit            save video data to file
    % connect           connect by checking existence of image file
    % disconnect        disconnect by resetting info property
    % getInfo           get info about video data
    % ________________________________________________________________________
    %
    % ToDo      1. fix bug in commit method
    %           2. ...
    %
    % Author    Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    properties (Constant = true)
        VIDEO_EXTENSION = [ ...
            ".avi", ...
            ".h265", ...
            ".mov"];
    end%properties
    
    methods
    % constructor %
        function wdio_video = WDIOVideo(url)
            % WDIOVideo  construct instance of this WDIOVideo class
            %
            % Construct an object of this WDIOVideo class to implement a 
            % data I/O interface to a video file. The full path of the video 
            % file is specified in the url argument. The url argument must 
            % have a supported extension.
            %
            % Arguments:
            %   url             full file path to video file
            %
            % Return:
            % wdio_video        instance of this WDIOVideo interface
            %
            % Throws:
            % ArgumentError     empty url argument
            % ArgumentError     unsupported video file format
            if isempty(url)
                throw(MException( ...
                    'WDIOVideo:ArgumentError', 'Empty filename.'));
            end%if
            % check url match with supported format
            is_supported = false;
            for video_extension = WDIOVideo.VIDEO_EXTENSION
                is_supported = (is_supported | endsWith(url, video_extension));
            end%for
            if ~is_supported
                throw(MException('WDIOVideo:ArgumentError', ...
                    'Unsupported video file format.'));
            end%if
            wdio_video = wdio_video@WDataIO(0, url);
        end%function
        
    % data I/O methods %
    function data = clone(wdio_video)
            % clone  load video data from file
            %
            % Load and return raw (unprocessed) video data from the file 
            % specified by the url property. The ffmpeg branch of the method
            % first extracts the video to a series of image files in a 
            % subdirectory of the url. Then reads the image files and collect
            % to an array of WDImage objects.
            %
            % Return
            % data      aray of WDImage object objects
            try % extract and read frames using ffmpeg
                % eventually create extraction target directory
                target = [wdio_video.url(1 : end - 5) filesep];
                mkdir(target);
                % extract image series from stream and count frames
                cmd = [ ...
                    'ffmpeg -hide_banner -loglevel quiet -i ' ...
                    wdio_video.url ' ' target 'f%06d.png'];
                dos(cmd);
                % read to image series
                frames = dir([target '*.png']);
                data = cell(numel(frames), 1);
                for f = 1:numel(frames)
                    data{f, 1} = WDImage(imread([target frames(f).name]), ...
                        frames(f).name);
                end%for
                % free disk space
                delete([target '*.png']);
                rmdir(target);
            catch % alternative: try video reader
                vr = VideoReader(wdio_video.url);
                data = cell(vr.NumFrames, 1);
                f = 0;
                while hasFrame(vr)
                    f = f + 1;
                    data{f, 1} = WDImage(readFrame(vr), sprintf('f%06d', f));
                end%while
            end%try
        end%function
        function commit(~, data, target)
            % commit  write data to target video file
            %
            % Write WDImage objects array provided in the data argument to the 
            % target video file. The value of the target argument may differ
            % from the url, but should be a valid video file name of supported
            % format.
            %
            % Arguments:
            % data              array of WDImage objects
            % target            url of target image file
            %
            % Throws:
            % ArgumentError     empty target filename argument
            % ArgumentError     unsupported target filename format
            if isempty(target) 
                throw(MException('WDIOHDF5:ArgumentError', 'Empty filename.'));
            end%if
            % check target url with supported format
            is_supported = false;
            for video_extension = WDIOVideo.VIDEO_EXTENSION
                is_supported = ...
                    (is_supported | endsWith(target, video_extension));
            end%for
            if ~is_supported
                throw(MException('WDIOVideo:ArgumentError', ...
                    'Unsupported target file format.'));
            end%if
            try % commit using ffmpeg
                tmp = ['.' filesep 'tmp' filesep];
                if ~isfolder(tmp)
                    mkdir(tmp);
                end%if
                % write image files
                for f = 1:data.getLength()
                    imwrite([tmp sprintf("f%06d", f)], data(f).getChannels());
                end%for
                % write image series to video file
                cmd = [ ...
                    'ffmpeg -hide_banner -loglevel quiet -i ' ...
                    tmp 'f%06d ' target];
                dos(cmd);
                delete([tmp '*.png']);                
            catch % alternative: commit using video writer
                vw = VideoWriter(target);
                open(vw);
                for f = data.C
                    writeVideo(vw, f);
                end%for
                close(vw);
            end%try
        end%function
        
    % connection methods %
        function is_connected = connect(wdio_video)
            % connect  connect by checking existence of video file
            %
            % Connect to the current video file specified in the url argument
            % by checking its existence. If the file could not be opened the
            % method returns false.
            %
            % Returns
            % is_connected      url specifies a valid video file
            is_connected = isfile(wdio_video.url);
        end%function
        function is_connected = disconnect(wdio_video)
            % disconnect disconnect by resetting info property
            %
            % Disconnect by clearing the info property. Set the return
            % value to false.
            %
            % Returns
            % is_connected      always equal to false
            wdio_video.info = [];
            is_connected = false;
        end%function
        
    % get methods %
        function info = getInfo(wdio_video)
            % getInfo  get info about video data
            % 
            % Get info about the entire video file and initialize the
            % info property with the info structure. An info structure is 
            % returned.
            %
            % Return
            % info      structure with info about video data
            try % to get video info using ffprobe
                logfile = ['.' filesep 'tmp' filesep 'ffprobe.log'];
                cmd = [ ...
                    'ffprobe -hide_banner -loglevel quiet ' ...
                    '-show_streams -select_streams v:0 -of flat ' ...
                    wdio_video.url ' > ' logfile];
                dos(cmd);
                % parse main video specs
                info = struct();
                lines = readlines(logfile, 'EmptyLineRule', 'skip');
                for line = lines'
                    [C, ~, e] = sscanf(line, 'streams.stream.%d.width=%d');
                    if isempty(e)
                        info.N_q = max(0, C(end));
                    end%if
                    [C, ~, e] = sscanf(line, 'streams.stream.%d.height=%d');
                    if isempty(e)
                        info.N_p = max(0, C(end));
                    end%if
                    [C, ~, e] = sscanf(line, ...
                        'streams.stream.%d.avg_frame_rate="%f/%d"');
                    if isempty(e)
                        info.f__Hz = C(2) / C(3);
                    end%if
                    [C, ~, e] = sscanf(line, 'streams.stream.%d.pix_fmt="%s"');
                    if isempty(e)
                        info.format = char(C(2:end - 1)');
                    end%if                        
                end%for                
            catch % alternative: get video info using video reader
                vr = VideoReader(wdio_video.url);
                info = struct( ...
                    'N_b', vr.BitsPerPixel, ...
                    'N_f', vr.NumFrames, ...
                    'N_p', vr.Height, ...
                    'N_q', vr.Width, ...
                    'f__Hz', vr.FrameRate);
            end%try
            % save info to property
            wdio_video.info = info;
        end%function
    end%methods
end%classdef


% end of module WDIOVideo

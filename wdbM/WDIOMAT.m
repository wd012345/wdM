classdef (Sealed = true) WDIOMAT < WDataIO
    % WDIOMAT  Data I/O Interface for MATLAB Data    
    % The WDIOMAT class is a final subclass of WDataIO and implement the data 
    % I/O interface to MATLAB *.mat file database. MATLAB is a commercial
    % software product of Mathworks (www.mathworks.com).
    %
    %
    % WDIOMAT Public Methods
    %   WDIOMAT                 construct an instance of this WDIOMAT object
    %   connect                 connect by checking existence of MAT file
    %   disconnect              disconnect by resetting info property
    %   getInfo                 get info about data organisation
    %   clone                   load data from source group or dataset
    %   commit                  write data to target group or dataset
    %
    %
    % author:           Stefan Wittwer
    % bug reports:      info@wittwer-datatools.ch
    % backlog items: 
    % 
    
    methods
        function wdio_mat = WDIOMAT(url)
            % WDIOMAT  construct an instance of this WDIOMAT class
            % Construct an object of this WDIOMAT class to implement a
            % data I/O to the MAT file specified in the url argument. The
            % url argument must have the correct extension ('*.mat').            
            % Arguments:
            %   url             path to MAT file
            % Return:
            % wdio_mat          instance of this WDIOMAT interface
            if ~endsWith(url, '.mat')
                throw(MException('WDIOMAT:ArgumentError', ...
                    'Empty filename or incorrect file extension.'));
            end%if
            wdio_mat = wdio_mat@WDataIO(0, url);
        end%function
        
        function raw = clone(wdio_mat, source)
            % clone  load data from source variable
            % Load and return raw data from the group or dataset specified 
            % by the source argument.
            if isempty(source)
                clones = load(wdio_mat.url);
                raw = WDFieldSet('MAT', clones);
            else
                clones = load(wdio_mat.url, source);
                raw = WDFieldSet(source, clones);
            end%if
        end%function
        
        function commit(wdio_mat, data, graph, video)
            % commit  write data to MAT file
            % Write data provided in the data argument to the MAT file
            % specified in the url property
            save(wdio_mat.url, 'data', '-v7.3');
            if ~isempty(graph)
                url = strrep(wdio_mat.url, '.mat', '.fig');
                savefig(graph(ishandle(graph)), url);
            end%if
            % do this via WDIOVideo interface in future...
            if ~isempty(video)
                url = strrep(wdio_mat.url, '.mat', '.avi');
                v = VideoWriter(url);
                open(v);
                writeVideo(v, video);
                close(v);
            end%if
        end%function
        
        function is_connected = connect(wdio_mat)
            % connect  connect by checking existence of MAT file 
            % Connect to the current MAT file specified in the url
            % property by checking its existence. If the file could not be
            % opened the method returns false.
            is_connected = isfile(wdio_mat.url);
        end%function
        
        function is_connected = disconnect(wdio_mat)
            % disconnect disconnect by resetting info property
            % Disconnect by clearing the info property. Set the return
            % value to false.
            wdio_mat.info = WDFieldSet('', []);
            is_connected = false;
        end%function
        
        function info = getInfo(wdio_mat)
            % getInfo  get MAT file object
            % Get a MAT file object which allows to access and change the
            % variables in the MAT file database directly.
            info = matfile(wdio_mat.url);
        end%function        
    end%method
    
end%classdef

classdef WDataIO < handle
    % WDataIO   generic data I/O interface
    %
    % The WDataIO is a generic data I/O interface object. It basically
    % provides
    % * the loading (cloning) of data from the database into the data tools
    %   objects, 
    % * the saving (committing) from the data tool objects back
    %   into the database.
    %
    % The database can be
    % * a non-SQL database with collections and documents,
    % * a SQL database with tables and records,
    % * a directory containing formatted text data files or binary data files,
    % * a compound data file (e.g. HDF5, video files) with subdata entries.
    % The implementation of a data I/O interface for a specific database is
    % done by subclassing WDataIO. There are abstract methods in this
    % WDataIO interface which have to be defined in each subclass.
    %
    %
    % Properties
    % ----------
    % info              structure with database info
    % port              port of database
    % url               url of database
    %
    %
    % Methods
    % -------
    % WDataIO           construct an instance of this WDataIO interface
    % getName           get name of this data I/O interface
    % getPath           get the directory path to the database
    % getURL            get url of the database
    %
    %
    % Abstract Methods
    % ----------------
    % connect           connect to database
    % disconnect        disconnect from database
    % clone             read from data source
    % commit            write to data target
    % getInfo           get info about database
    %
    % ________________________________________________________________________
    %
    % ToDo
    % 1. ...
    %
    % Authors:
    % Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    properties (Access = protected)
        info                % data fields with info on data structure
        port    uint16      % port of database (0-65535)
        url                 % name or file path of database
    end
    
    methods
    % constructor %
        function wd_io = WDataIO(port, url)
            % WDataIO  construct an instance of this WDataIO interface
            %
            % Construct an interface instance to the database or
            % directory provided in the arguments.
            %
            % Arguments:
            % port              port number of non-SQL or SQL database
            % url               url or path to database or directory
            %
            % Returns:
            % an instance of this WDataIO interface
            wd_io.port = port;
            wd_io.url = url;
        end
        
    % get methods %
        function name = getName(wd_io)
            % getName  get name of database
            %
            % By definition the name of the database (source or target) is 
            % the tail (last part) of the url.
            %
            % Returns:
            % name of this data I/O interface
            name = strsplit(wd_io.url, filesep);
            name = name{end};
        end%function
        function path = getPath(wd_io)
            % getPath  get directory path to the database
            %
            % By definition the base path to the database is the url minus the
            % last part. This is the directory of the database.
            %
            % Returns:
            % the directory path to the database
            path = strsplit(wd_io.url, filesep);
            path = strjoin(path(1 : end - 1), filesep);
            path = [path filesep];
        end%function            
        function url = getURL(wd_io)
            % getURL  get url of the database
            %
            % Return the value of the url property containing the database
            % url.
            %
            % Returns:
            % url of the database
            url = wd_io.url;
        end%function
    end%methods
    
    methods (Abstract)
        data = clone(wd_io, source)             % read from data source
        commit(wd_io, data, target)             % write to data target
        is_connected = connect(wd_io)           % connect to database
        is_connected = disconnect(wd_io)        % disconnect from database
        info = getInfo(wd_io)                   % return info about database
    end%methods
end%classdef


% end of module wdM.wdbM.WDataIO

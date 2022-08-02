classdef WDIOHDF5 < WDataIO
    % WDIOHDF5  Data I/O Interface for HDF5 Data
    %
    % The WDIOHDF5 class is a subclass of WDataIO and implement the data 
    % I/O interface to a HDF5-formatted file database. HDF is the 
    % hierarchical data format supported by the HDF Group (www.hdfgroup.org).
    %
    %     /                     root group
    %     /group                group as member of root
    %     /group/dataset        dataset as member of group (field set)
    %     /dataset              dataset (field set) with fields (features)
    %
    % The data of a HDF5 file database model is organized in groups and
    % datasets. The top group is the root group (\) followed by
    % zero or more group levels. A group is either followed by other groups
    % or by datasets.
    %
    %
    % Methods
    % -------
    % WDIOHDF5          construct an instance of this WDIOHDF5 object
    % connect           connect by checking existence of HDF5 file
    % disconnect        disconnect by resetting info property
    % getInfo           get info about data organisation
    % clone             load data from source group or dataset
    % commit            write data to target group or dataset
    % ________________________________________________________________________
    %
    % ToDo      1. ...
    %
    % Author    Stefan Wittwer, info@wittwer-datatools.ch
    % ________________________________________________________________________
    
    methods
    % constructor %
        function wdio_hdf5 = WDIOHDF5(url)
            % WDIOHDF5  construct an instance of this WDIOHDF5 class
            %
            % Construct an object of this WDIOHDF5 class to implement a
            % data I/O to the HDF5 file specified in the url argument. The
            % url argument must have the correct extension ('*.h5').
            %
            % Arguments:
            % url               path to HDF5 file
            %
            % Return:
            % wdio_hdf5         instance of this WDHDF5IO interface
            %
            % Throws:
            % exception on empty filename of incorrect file extension
            if ~endsWith(url, '.h5')
                throw(MException('WDIOHDF5:ArgumentError', ...
                    'Empty filename or incorrect file extension.'));
            end%if
            wdio_hdf5 = wdio_hdf5@WDataIO(0, url);
        end%function
        
    % data I/O methods %
    function data = clone(wdio_hdf5, source)
            % clone  load data from source group or dataset
            %
            % Load and return raw data from the group or dataset specified 
            % by the source argument.
            %
            % Arguments
            % source    name of data source
            %
            % Returns
            % data      raw data
            data = h5read(wdio_hdf5.url, source);
        end%function
        function commit(wdio_hdf5, data, target)
            % commit  write data to target group or dataset
            %
            % Write data provided in the data argument to the group or
            % dataset specified by the target argument.
            %
            % Arguments
            % data      target data
            % target    name of target dataset
            h5write(wdio_hdf5.url, target, data);
        end%function
        
    % connection methods %
        function is_connected = connect(wdio_hdf5)
            % connect  connect by checking existence of HDF5 file 
            %
            % Connect to the current HDF5 file specified in the url
            % argument by checking its existence. If the file could not be
            % opened the method returns false.
            %
            % Returns
            % is_connected      url specifies a valid HDF5 database file
            is_connected = isfile(wdio_hdf5.url);
        end%function
        function is_connected = disconnect(wdio_hdf5)
            % disconnect disconnect by resetting info property
            %
            % Disconnect by clearing the info property. Set the return
            % value to false.
            %
            % Returns
            % is_connected      always equal to false
            wdio_hdf5.info = [];
            is_connected = false;
        end%function
        
    % get methods %
        function info = getInfo(wdio_hdf5, group)
            % getInfo  get info about database organisation
            %
            % Get info about the entire HDF5 file if the group argument is
            % empty. Otherwise get info about the group specified by the
            % argument. In both cases an info structure is returned.
            %
            % Arguments
            % group     dataset or data group
            %
            % Returns
            % info      structure with database info
            if isempty(group)
                group = '/';
            end%if
            info = h5info(wdio_hdf5.url, group);
            wdio_hdf5.info = info;
        end%function
    end%methods
end%classdef


% end of module wdM.wdbM.WDIOHDF5

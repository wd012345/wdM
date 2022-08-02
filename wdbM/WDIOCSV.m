classdef WDIOCSV < WDataIO
    % WDIOCSV   data I/O interface for CSV formatted text files
    % ...
    %
    %
    % ToDo:     ...
    %
    % Author:   Stefan Wittwer, info@wittwer-datatools.ch
    
    properties (Access = private, Constant = true)
        separator = ',';
        filetype = ".csv";
    end%properties
    
    methods
    % constructor %
        function wdio_csv = WDIOCSV(port, url)
            wdio_csv = wdio_csv@WDataIO(port, url);
            if ~isfolder(url)
                message = "URL argument must be valid path to CSV files.";
                throw(MException("WDIOCSV:InvalidURL", message));
            end%if
        end%function
        
    % I/O methods %
        function wd_field = clone(wd_io, source)
            if endsWith(source, WDIOCSV.filetype)
                data_string = wd_io.info.get;
                
            end%if
%                         try
%                 line = split(line(2:end), WDIOCSV.separator);
%             catch ME
%                 rethrow(ME);
%             end%try
        end%function
        
    % get methods %
        function info = getInfo(wdio_csv)
            line = readlines(wdio_csv.url, "EmptyLineRule", "skip");
            header = split(line(1), WDIOCSV.separator);
            N_f = numel(header);
            N_r = size(line, 1) - 1;
            info = WDField(header, "HeaderString");
            info = [info WDField(line(2:end), "DataString")];
            info = [info WDField(N_f, "NumberOfHeaderFields")];
            info = [info WDField(N_r, "NumberOfRecords")];
            wdio_csv.info = info;
        end%function
    end%methods
end%classdef


% end of module wdbM.WDIOCSV

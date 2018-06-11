function [data] = getFile(InputFile)
%getFile Converts a text file to a cell array where the line of each file
%        is stored in a row of the cell array.
%   InputFile: Full path of the file to be parsed (or just name assuming
%              function is executed in working directory)
%   data:      file stored as cell array, or - if file doesn't exist -
%              returns and empty matrix
fid = fopen(InputFile);
try
    output = textscan(fid,'%s','delimiter','\n','whitespace','');
    fclose(fid);
    %Extract contents of cell
    data = output{1};
catch
    %If file does not exist, return empty matrix
    data = [];    
end
end


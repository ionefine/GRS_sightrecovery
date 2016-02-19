function [vertex,x,y,z]=readCortexLabels(filename, startRow, endRow)
% read in a freesurfer label file
% written by IF 2015/12/27 

%% Initialize variables.
delimiter = {',',' '};
if nargin<=2
    startRow = 3;
    endRow = inf;
end

formatSpec = '%f%s%f%s%f%s%f%f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
vertex = dataArray{:, 1};
x = dataArray{:, 3};
y = dataArray{:, 5};
z = dataArray{:, 7};




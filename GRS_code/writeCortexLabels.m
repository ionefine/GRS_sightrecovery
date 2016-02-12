function writeCortexLabels(filename, data)

%% Initialize variables.
delimiter = {' '};

%% Open the text file.

fp = fopen(filename,'wt');
fprintf(fp, '%s','#!ascii label  , from subject mmSR20151113 vox2ras=TkReg');
fprintf(fp, '\n%d\n', length(data.vertices));
for n=1:length(data.vertices)
    fprintf(fp, '%d%s%6.3f%s%6.3f%s%6.3f%s%1.10f\n', data.vertices(n), delimiter{1}, data.vertXYZ(n, 1), delimiter{1}, ...
        data.vertXYZ(n, 2), delimiter{1}, data.vertXYZ(n, 3), delimiter{1}, 0);
end

%% Close the text file.
fclose(fp);



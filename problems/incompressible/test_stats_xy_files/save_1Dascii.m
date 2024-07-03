function [] = save_1Dascii(data,filename)
    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Check if the file was successfully opened
    if fileID == -1
        error('File could not be opened.');
    end
    
    % Write each element of the data array to the file
    for i = 1:length(data)
        fprintf(fileID, '%.16f\n', data(i));
    end
    
    % Close the file
    fclose(fileID);
    
    % Confirm writing success
    disp(['Data successfully written to ', filename]);

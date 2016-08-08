function mat2bin
% function mat2bin
% A function that stores the data from kdd_cup_1999.mat
% in a binary file column-wise.
% 
% author: Nikos Sismanis
% date: Jan 2014

%load('kdd_cup_1999.mat'); % load the data from a .mat file
%data = x.corpus;



fid = fopen('data524288.bin', 'w'); % open a binary file
fwrite(fid, d, 'single'); % write data to file
fclose(fid); % close the file

end
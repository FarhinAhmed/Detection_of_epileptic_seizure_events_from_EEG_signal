%% read the edf files and convert them to matfiles

clc
clear 
close all

files=dir('*.edf');
for i=1:length(files)
    filename=files(i).name;
    [hdr,record]=edfread(filename);
    filename=erase(filename,'.edf')
%     f = sprintf('%s.mat',filename);
    f=strcat(filename,'.mat')
    save(f,'record')
end
function [mriStack, info] = readDicom(dir)
%READDICOM Returns both DICOM images and meta data for each image
%   The user has to choose a directory which contains the DICOM files. It
%   does not matter if the files are sorted or not as the method checks the
%   instance number of each slice and sorts them, as well as the
%   metadata, into the sequences.

if ~exist('dir', 'var')
    disp("No path given as argument. Please choose directory with DICOM data.")
    sourcetable = dicomCollection(uigetdir);
else
    sourcetable = dicomCollection(dir);
end
rows = sourcetable.Rows;
cols = sourcetable.Columns;
frames = sourcetable.Frames;
files = sourcetable.Filenames;

info = cell(1, frames);
mriStack = zeros(double([rows cols frames]));
for i = 1:frames
    mrInfo = dicominfo(files{1, 1}{i, 1});
    info{1, mrInfo.InstanceNumber} = mrInfo;
    mriStack(:,:, mrInfo.InstanceNumber) = dicomread(files{1, 1}{i, 1});
end
% Normalize gray values
mriStack = mriStack - min(mriStack(:));
mriStack = mriStack ./ max(mriStack(:));
end

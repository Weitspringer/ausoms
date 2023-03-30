%% Parameters
executeControlPointReg = 0;
executeAutoReg = 0;

%% Read DICOM images
folder = strcat(uigetdir(), "\"); % Choose folder containing MRI sequence
mrSlices = dir(folder);
mrSlices([1 2]) = [];
mr = zeros(256, 256, length(mrSlices), 'double');
sliceLocations = double.empty(length(mrSlices), 0);
for i = 1:length(mrSlices)
    mrInfo = dicominfo(strcat(folder, mrSlices(i).name));
    mr(:,:,mrInfo.AcquisitionNumber) = dicomread(strcat(folder, mrSlices(i).name));
    sliceLocations(i) = mrInfo.SliceLocation;
    if i == 3
        pixelSpacing = mrInfo.PixelSpacing;
        rows = mrInfo.Rows;
        columns = mrInfo.Columns;
    end
end
[file, path] = uigetfile('*.*'); % Choose the corresponding CT file
ct_file = strcat(path, file);
ct = squeeze(dicomread(ct_file));

%% Interpolation of the MR slices
% Caution: Assumption that the CT has a isometric resolution of 0.125/0.125/0.125 (mm) 
[x0, y0, z0] = ndgrid(0:pixelSpacing(1):double(columns-1)*pixelSpacing(1),0:pixelSpacing(2):double(rows-1)*pixelSpacing(2), sliceLocations);
[x1, y1, z1] = ndgrid(0:0.125:double(columns-1)*pixelSpacing(1),0:0.125:double(rows-1)*pixelSpacing(2),sliceLocations(1):0.125:sliceLocations(end));
mrtInterp = interpn(x0,y0,z0,mr,x1,y1,z1, 'makima');

%% Registration

% Registration with control points
% Note that it could be necessary to adjust the unregistered slices as the
% heads in the ct images are all positioned differently
if executeControlPointReg == 1
    sizeMrtInterp = size(mrtInterp);
    sizeCt = size(ct);
    
    % Axial Registration
    fixedAxial = mat2gray(mrtInterp(:,:,round(sizeMrtInterp(3)/2)));
    unregisteredAxial = ct(:,:,278);
    [mp,fp] = cpselect(unregisteredAxial,fixedAxial,'Wait',true);
    transformationAxial = fitgeotrans(mp,fp,'similarity'); % 2D
    RfixedAxial = imref2d(size(fixedAxial));
    registeredAxial = imwarp(unregisteredAxial,transformationAxial,'OutputView',RfixedAxial);
    
    % Use transformationAxial to warp all axial slices
    ctAxialRegistered = zeros([size(registeredAxial) sizeCt(3)], 'int16');
    for i = 1:sizeCt(3)
        % Slices are placed back exactly as they were taken out of the
        % volume
        axialSlice = ct(:,:,i);
        ctAxialRegistered(:,:,i) = imwarp(axialSlice,transformationAxial,'OutputView',RfixedAxial);
    end
    
    % Coronal Registration
    sizeCtAxialRegistered = size(ctAxialRegistered);
    fixedCoronal = squeeze(mat2gray(mrtInterp(round(sizeMrtInterp(1)/2),:,:)));
    unregisteredCoronal = squeeze(ctAxialRegistered(round(sizeCtAxialRegistered(1)/2)+10,:,:));
    [mp2,fp2] = cpselect(unregisteredCoronal,fixedCoronal,'Wait',true);
    transformationCoronal = fitgeotrans(mp2,fp2,'similarity');
    RfixedCoronal = imref2d(size(fixedCoronal));
    registeredCoronal = imwarp(unregisteredCoronal,transformationCoronal,'OutputView',RfixedCoronal);
    
    % Use transformationCoronal to warp all coronal slices
    ctCoronalRegistered = zeros([sizeCtAxialRegistered(1) size(registeredCoronal)], 'int16');
    for i = 1:sizeCtAxialRegistered(1)
        coronalSlice = ctAxialRegistered(i,:,:);
        ctCoronalRegistered(i,:,:) = imwarp(squeeze(coronalSlice),transformationCoronal,'OutputView',RfixedCoronal);
    end
    
    % Sagittal Registration
    sizeCtCoronalRegistered = size(ctCoronalRegistered);
    fixedSagittal = squeeze(mat2gray(mrtInterp(:,round(sizeCtCoronalRegistered(2)/2),:)));
    unregisteredSagittal = squeeze(ctCoronalRegistered(:,round(sizeCtCoronalRegistered(2)/2),:));
    [mp3,fp3] = cpselect(unregisteredSagittal,fixedSagittal,'Wait',true);
    transformationSagittal = fitgeotrans(mp3,fp3,'similarity');
    RfixedSagittal = imref2d(size(fixedSagittal));
    registeredSagittal = imwarp(unregisteredSagittal,transformationSagittal,'OutputView',RfixedSagittal);
    
    % Use transformationSaggital to warp all sagittal slices
    sizeRegisteredSagittalReference = size(registeredSagittal);
    sizeRegisteredCoronalReference = size(registeredCoronal);
    ctSagittalRegistered = zeros([sizeRegisteredSagittalReference(1) sizeRegisteredCoronalReference(2) sizeRegisteredSagittalReference(2)], 'int16');
    for i = 1:sizeCtCoronalRegistered(2)
        sagittalSlice = ctCoronalRegistered(:,i,:);
        ctSagittalRegistered(:,i,:) = imwarp(squeeze(sagittalSlice),transformationSagittal,'OutputView',RfixedSagittal);
    end
    
    ctControlPointRegistered = ctSagittalRegistered;
    figure; imshowpair(fixedAxial, registeredAxial);
    title('Control Point Guided Registration - Axial');
    figure; imshowpair(fixedCoronal, registeredCoronal);
    title('Control Point Guided Registration - Coronal');
    figure; imshowpair(fixedSagittal, registeredSagittal);
    title('Control Point Guided Registration - Sagittal');
end

%% Automatic Registration
if executeAutoReg == 1
    fixedVolume = mrtInterp;
    movingVolume = ctControlPointRegistered;
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 0.005;
    optimizer.GrowthFactor = 1.01;
    Rfixed  = imref3d(size(fixedVolume),0.125,0.125,0.125);
    Rmoving = imref3d(size(movingVolume),0.125,0.125,0.125);
    movingRegisteredVolume = imregister(movingVolume,Rmoving, fixedVolume,Rfixed, 'rigid', optimizer, metric);
    figure; imshowpair(squeeze(fixedVolume(:,80,:)), squeeze(movingRegisteredVolume(:,80,:)));
    figure; imshowpair(squeeze(fixedVolume(80,:,:)), squeeze(movingRegisteredVolume(80,:,:)));
    figure; imshowpair(fixedVolume(:,:,60), movingRegisteredVolume(:,:,60));
end

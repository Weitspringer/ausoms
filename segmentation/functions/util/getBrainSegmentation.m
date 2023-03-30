function [brain] = getBrainSegmentation(mriStack, pixelSpacing, sliceThickness, SEsize, repetitions)
%GETBRAINSEGMENTATION Retrieve rodent brain tissue from the MR volume
%   Uses the 3D PCNN method introduced by N. Chou, J. Wu, J. Bai Bingren, 
%   A. Qiu and K. Chuang, "Robust Automatic Rodent Brain Extraction Using 
%   3-D Pulse-Coupled Neural Networks (PCNN)," in IEEE Transactions on 
%   Image Processing, vol. 20, no. 9, pp. 2554-2564, Sept. 2011, 
%   doi: 10.1109/TIP.2011.2126587
%   Code retrieved by usage of web.archive.org:
%   (1) https://web.archive.org/web/20151031152204/http://www.a-star.edu.sg/sbic/RESOURCES/Software.aspx
%   Last checked: 30.03.2021

%   The method uses bias field correction beforehand to improve the
%   results. In our use case, the voxel dimensions were distorted. The
%   result of the brain segmentation is postprocessed to compensate errors.
%   The 3D-PCNN method also introduces a random component. To compensate
%   the variance, the optimal solution is determined automatically.

%   === Parameters ===
%   volume - Normalized grayscale volume containing the slices of the MRI sequence
%   pixelSpacing - 2x2 array containing the size of a pixel (in mm)
%                  Attribute "pixel spacing" of DICOM meta data
%   sliceThickness - Attribute "slice thickness" of DICOM meta data to
%                    calculate the volume of a voxel
%   SEsize - Size of the structuring element for the 3D-PCNN segmentation
    
    mriStackCorr = biasFieldCorrection(mriStack); % Bias field correction
    mriStackCast = mriStackCorr .* 255;
    % 3D PCNN
    maxIterations = 200;
    optimalBrainVolume = [100 550]; % Default brain size range for mice. Try [1200 4400] for adult rat brain (1)
    voxelVolume = [pixelSpacing(1), pixelSpacing(2), sliceThickness];
    pcnn_segmentations = cell(1, repetitions);
    % Repeat the 3D PCNN segmentation to compensate variance
    for idx = 1:repetitions
        [border, G_I] = PCNN3D(mriStackCast, SEsize, voxelVolume, optimalBrainVolume, maxIterations);
        close;
        its = findit(G_I, optimalBrainVolume);
        bestSolution = border{1, its}; % Find best iteration of the algorithm.

        brain = zeros(size(mriStackCorr, 1), size(mriStackCorr, 2), size(mriStackCorr, 3));
        % Postprocess the segmentation and only choose the brain segment in
        % each slice. Utilization of the mouse head symmetry. 
        for i = 1:size(brain, 3)
            brainMask = full(bestSolution{i});
            try
                [y, x] = ndgrid(1:size(brainMask, 1), 1:size(brainMask, 2));
                centroid = mean([x(logical(brainMask)), y(logical(brainMask))]);
                row = centroid(2);
                col = centroid(1);
                brain(:,:,i) = bwselect(brainMask, col, row);
            catch
                brain(:,:,i) = brainMask;
            end
        end
        pcnn_segmentations{1, idx} = brain;
    end
    
    % Choose the best slices of the different 3D-PCNN brain segmentations
    % automatically
    for i = 1:size(mriStack, 3)
        avg_ones = 0;
        for j = 1:repetitions
            seg = pcnn_segmentations{1, j};
            slice = seg(:,:,i);
            avg_ones = avg_ones + nnz(slice);
        end
        avg_ones = avg_ones / repetitions;
        diffs = zeros(1, repetitions);
        for j = 1:repetitions
            seg = pcnn_segmentations{1, j};
            slice = seg(:,:,i);
            ones = nnz(slice);
            diffs(j) = abs(avg_ones - ones);
        end
        [~,best] = min(diffs);
        vol = pcnn_segmentations{1, best};
        brain(:,:,i) = vol(:,:,i);
    end
    
    brain = cast(brain, 'logical');
end

function [combination] = combineSegmentations(foregroundSegmentation, brainSegmentation, pixelSpacing, maxSkullThickness)
%COMBINESEGMENTATIONS Combine foreground segmentation and brain segmentation

% Create a region of interest (ROI) around the brain
dilated = zeros(size(brainSegmentation));
ps = pixelSpacing(1);
parfor i = 1:size(brainSegmentation, 3)
    brain = brainSegmentation(:,:,i);
    SEsize = ceil(maxSkullThickness / ps);
    SE = strel('disk', SEsize);
    dilatedBrain = imdilate(brain, SE);
    dilated(:,:,i) = dilatedBrain;
end
roi = dilated & ~brainSegmentation;
% Cut out dark structures of the brain in the foreground segmentation,
% then select only pixels in the ROI of the foreground segmentation.
combination = (foregroundSegmentation & ~brainSegmentation) & roi;
end


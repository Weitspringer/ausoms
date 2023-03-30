function [foregroundSegmentation] = getForegroundSegmentation(images, SEBg, SEClose, nSegments)
%GETFOREGROUNDSEGMENTATION Remove background from image and segment the
%dark structures of the mouse head in each MRI slice. The foreground is
%defined as the mouse head in the MR image.
%   At first, a morphological background removal is performed on the input
%   image. Closing the binarized foreground provides a rough binary image
%   which fills potentially underestimated areas. To counter the errors in
%   this step, the next step is to binarize the input image adaptively and
%   combine these two results.
%   Before the two masks are combined, the second mask gets inverted. The
%   overlap of these to images is the segmented foreground.

foregroundSegmentation = cast(zeros(size(images)), 'logical');
parfor i = 1:size(images, 3)
    %% Morphological background removal
    I = im2uint8(images(:,:,i)); % imbinarize uses 256 bins
    Ibg = imclose(I, SEBg);
    ISub = Ibg - I; % Retrieve dark structures of the slice which were supplanted
    BWSub = imclose(imbinarize(ISub), SEClose);

    %% Adaptive binarizing of image, foreground specified as dark
    BW = ~imbinarize(I, 'adaptive', 'ForegroundPolarity', 'dark');

    %% Combine BWSub and BW to further improve the mask
    foregroundSegmentation(:,:,i) = bwareafilt(BW & BWSub, nSegments);
end
end


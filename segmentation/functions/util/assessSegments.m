function [assessment] = assessSegments(segmentation, pixelSpacing)
%ASSESSSEGMENTS Exclude segments by means of shape
%   The method uses the maximum Feret diameter to distinguish between skull
%   segments and small fragments.

assessment = segmentation;
ps = pixelSpacing(1);
parfor i = 1:size(assessment, 3)
    slice = assessment(:,:,i);
    cc = bwconncomp(slice); % Get connected components (segments)
    [out,~] = bwferet(slice,'MaxFeretProperties'); % Max Feret properties
    diameters = out.MaxDiameter(:);
    minDiameter = 4; % in mm
    g = diameters >= round(minDiameter / ps); % Filter out fragments
    segments = cc.PixelIdxList(g);
    corr = zeros(size(slice)); % New empty slice
    % "Paste" the filtered segments onto the empty slice
    for j = 1:size(segments, 2)
        points = segments{1,j};
        corrIter = zeros(size(slice));
        corrIter(points) = 1;
        corr = corr | corrIter;
    end
    assessment(:,:,i) = corr;
end
end


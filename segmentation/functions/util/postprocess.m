function [postprocessed] = postprocess(segmentation, pixelSpacing)
%POSTPROCESS Postprocessing of the segmentation result
%   The postprocessing attempts to close small gaps between 2 segments.
%   They can be due to partial volume artifacts (PVEs) or other types of
%   artifacts (e.g. movement artifacts).
%   Each slice is checked if it has 2 segments. If there are 2 segments,
%   the gap between the segments is assessed. The postprocessing
%   approximates the skull in small gaps.
%
%   pixelSpacing is an 1x2 array containing the pixel dimensions in the x
%   and y dimension [x,y].

postprocessed = segmentation;
ps = pixelSpacing(1); % Assumption: The slices are isometric
for i = 1:size(postprocessed, 3)
    % Connect segments
    Ioriginal = double(postprocessed(:,:,i));
    cc = bwconncomp(Ioriginal);
    if cc.NumObjects == 2
        points = cc.PixelIdxList;
        seg1 = zeros(size(Ioriginal));
        seg2 = zeros(size(Ioriginal));
        seg1(points{1, 1}) = 1; % Slice showing the first segment
        seg2(points{1, 2}) = 1; % Slice showing the second segment

        % Get distances to the corresponding segments in each image
        d1 = bwdist(seg1);
        d2 = bwdist(seg2);
        
        d = d1 + d2; % Combine distances to paths via all pixels in the slice
        minval = min(d(:)); % Shortest path length
        mmthresh = 1.5; % gap maximum (mm)
        goodness = 0.5; % deviation of the distance from the direct path (px)
        
        % Classify additional pixels as 1 if gap is small enough
        if minval <= ceil(mmthresh / ps)
            [minrows, mincols] = find(d <= minval + goodness);
            indizes = sub2ind(size(Ioriginal), minrows, mincols);
            Ioriginal(indizes) = 1;
            postprocessed(:,:,i) = Ioriginal;
        end
    end
end
end


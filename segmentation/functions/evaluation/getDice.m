function [dice] = getDice(segmentation, groundTruth, visualize)
%GETDICE Quantitative evaluation of the segmentation with the Dice metric.
%   Both segmentation and groundTruth have to be provided as logical 2D
%   image sequences.
%   Also, the result is visualized if needed.
if ~exist('visualize', 'var')
    visualize = false;
end

if visualize
    [~, ~, frames] = size(segmentation);
    frameScores = zeros(1, frames);
    for i = 1:frames
        seg = segmentation(:,:,i);
        gT = groundTruth(:,:,i);
        TP = nnz(seg & gT);
        FP = nnz(seg & ~gT);
        FN = nnz(~seg & gT);
        frameScores(i) = 2*TP / (2*TP + FP + FN);
    end
    figure;
    plot(frameScores);
    legend("dice", 'Location', 'south');
    title("Dice score (separate slices)");
    xlabel("Slice");
end
TP = nnz(segmentation & groundTruth);
FP = nnz(segmentation & ~groundTruth);
FN = nnz(~segmentation & groundTruth);
dice = 2*TP / (2*TP + FP + FN); % Formula
end


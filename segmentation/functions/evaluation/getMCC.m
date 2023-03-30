function [mcc] = getMCC(segmentation, groundTruth, visualize)
%GETMCC  Quantitative evaluation of the segmentation with the MCC metric.
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
        TN = nnz(~seg & ~gT);
        FN = nnz(~seg & gT);
        % Catch edge cases if only one class is non-zero
        if (TP == 0 && TN ~= 0 && FP == 0 && FN == 0) || (TP ~= 0 && TN == 0 && FP == 0 && FN == 0)
            frameScores(i) = 1.0; % Perfect classification
        elseif (TP == 0 && TN == 0 && FP ~= 0 && FN == 0) || (TP == 0 && TN == 0 && FP == 0 && FN ~= 0)
            frameScores(i) = -1.0; % Perfect misclassification
            % Otherwise, catch edge cases if exactly 2 classes are non-zero
        elseif (FN == 0 && TN == 0) || (FP == 0 && TN == 0) || (TP == 0 && FN == 0) || (TP == 0 && FP == 0)
            frameScores(i) = 0.0; % Not better than a toin coss
        else
            % Only if nothing of the above applies, use the formula
            frameScores(i) = (TP * TN - FP * FN) / (sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)));
        end
    end
    figure;
    plot(frameScores);
    legend("mcc", 'Location', 'south');
    title(strcat("Matthews Correlation Coefficient (separate slices)"));
    xlabel("Slice");
end
TP = nnz(segmentation & groundTruth);
FP = nnz(segmentation & ~groundTruth);
TN = nnz(~segmentation & ~groundTruth);
FN = nnz(~segmentation & groundTruth);
% Catch edge cases if only one class is non-zero
if (TP == 0 && TN ~= 0 && FP == 0 && FN == 0) || (TP ~= 0 && TN == 0 && FP == 0 && FN == 0)
    mcc = 1.0; % Perfect classification
elseif (TP == 0 && TN == 0 && FP ~= 0 && FN == 0) || (TP == 0 && TN == 0 && FP == 0 && FN ~= 0)
    mcc = -1.0; % Perfect misclassification
    % Otherwise, catch edge cases if exactly 2 classes are non-zero
elseif (FN == 0 && TN == 0) || (FP == 0 && TN == 0) || (TP == 0 && FN == 0) || (TP == 0 && FP == 0)
    mcc = 0.0; % Not better than a toin coss
else
    % Only if nothing of the above applies, use the formula
    mcc = (TP * TN - FP * FN) / (sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)));
end
end


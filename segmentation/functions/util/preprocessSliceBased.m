function [stack] = preprocessSliceBased(inputStack, pixelSpacing)
%PREPROCESSSLICEBASED The preprocessing improves the contrast of each slice
%in the MRI sequence.
%   inputStack - Normalized stack of MRI slices
%   pixelSpacing - 1x2 array containing the pixel dimensions

stack = inputStack;
ps = pixelSpacing(1); % Assumption: The separate slices are isometric
parfor i = 1:size(stack, 3)
    I = stack(:,:,i);
    tileSize = 1; % in mm
    tilesRow = round(size(inputStack, 1) / (tileSize / 0.0781));
    tilesCol = round(size(inputStack, 2) / (tileSize / 0.0781));
    histEq = adapthisteq(I, 'Distribution', 'rayleigh', 'ClipLimit', 0, 'NumTiles', [tilesRow tilesCol]);
    I = imadjust(histEq);
    % Heuristic: 0.8591 mm for searchWinSize and 0.2342 mm for compWinSize
    searchWinSize = 2*round(0.8591 / ps / 2) + 1; % Round to next odd value (expected by imnlmfilt)
    compWinSize = 2*round(0.2343 / ps / 2) + 1;
    stack(:,:,i) = imnlmfilt(I, 'SearchWindowSize', searchWinSize,  'ComparisonWindowSize', compWinSize);
end
end

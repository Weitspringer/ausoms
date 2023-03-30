function [corrected] = removeOutliers(segmentation)
%REMOVEOUTLIERS Make use of neighboring information between the slices to
%remove outliers
%   Check each segment if it is overlapping with at least one other segment
%   in at least one neighboring slice. If this is not the case, remove the
%   segement. Write corrected slice directly back into the volume. The
%   algorithm makes use of the geodesic distance of binary images.

corrected = segmentation;
for i = 1:size(corrected, 3)
    if size(corrected, 3) > 1
        if i == 1 % First slice has only one neighbor
            current = corrected(:,:,i);
            next = corrected(:,:,i+1);
            distNext = bwdistgeodesic(current, next);
            D = distNext == inf;
            corrected(:,:,i) = current .* ~D;
        else
            if i < size(corrected)
                last = corrected(:,:,i-1); % Previous neighbor
                current = corrected(:,:,i); % Slice whose segments are assessed
                next = corrected(:,:,i+1); % Next neighbor
                distLast = bwdistgeodesic(current, last);
                distNext = bwdistgeodesic(current, next);
                D = (distLast == inf) & (distNext == inf); % Pixels of outliers: inf values (value of background pixels: NaN)
                corrected(:,:,i) = current .* ~D;
            else
                if i == size(corrected, 3) % Last slice has only one neighbor
                    last = corrected(:,:,i-1);
                    current = corrected(:,:,i);
                    distLast = bwdistgeodesic(current, last);
                    D = distLast == inf;
                    corrected(:,:,i) = current .* ~D;
                end
            end
        end
    end
end
end


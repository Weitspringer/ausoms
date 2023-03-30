function [eval] = evaluate(skullSeg, sGt, beta, brainSeg, bGt)
%EVALUATE Calculate scores for segmentation(s) and store them in struct
%   skullSeg: Segmentation volume of skull, logical
%   sGT: Skull ground truth volume, logical
%   beta: Beta parameter for F-Score
%   brainSeg: Segmentation volume of brain, logical
%   bGt: Brain ground truth volume, logical
%   The result is a struct which contains the evaluation for the skull and,
%   if the parameters were given, the evaluation for the brain.

assert(isequal(size(skullSeg), size(sGt)));
eval = struct();
if isa(skullSeg, 'logical') && isa(sGt, 'logical')
    skull = struct();
    skull.whole.fScore = getFScore(skullSeg, sGt, beta);
    skull.whole.mcc = getMCC(skullSeg, sGt);
    skull.whole.dice = getDice(skullSeg, sGt);
    halfSize = round(size(skullSeg, 1) / 2);
    skull.upperHalf.fScore = getFScore(skullSeg(1:halfSize,:,:), sGt(1:halfSize,:,:), beta);
    skull.upperHalf.mcc = getMCC(skullSeg(1:halfSize,:,:), sGt(1:halfSize,:,:));
    skull.upperHalf.dice = getDice(skullSeg(1:halfSize,:,:), sGt(1:halfSize,:,:));
    eval.skull = skull;
    if exist('brainSeg', 'var') && exist('bGt', 'var') 
        if isa(brainSeg, 'logical') && isa(bGt, 'logical')
            assert(isequal(size(brainSeg), size(bGt)));
            brain = struct();
            brain.fScore = getFScore(brainSeg, bGt, beta);
            brain.mcc = getMCC(brainSeg, bGt);
            brain.dice = getDice(brainSeg, bGt);
            eval.brain = brain;
        else
            error("One of the input brain volumes for evaluation is not a logical volume")
        end
    end
else
    error("One of the input skull volumes for evaluation is not a logical volume")
end
end


function [skull, brain, eval] = segmentVolume(dir, options)
%SEGMENTVOLUME Automatic segmentation of the mouse skull in MRI DICOM data
%   It is recommended to use slices in the transversal plane. The skullcap
%   should be visible in every slice. The segmentation is fully automatic.
%   However, feel free to adjust some options to adjust the method to your
%   use case.

%   === Mandatory ===
%   dir - Directory containing DICOM sequence / volume
%   options = struct

%   === Optional ===
%   SEBg - Structuring element used for background removal
%   SEBgClose - Structuring element used for closing of foreground
%   sizeSEBrain - Structuring element size for 3D PCNN brain segmentation
%   resizeFactor - Resizing factor of volume. Default: 1
%   cutoff - Number of slices ignored at the start as well as the end
%            of the sequence. Default: 0
%   nSegmentsForeground - Number of segments for the foreground.
%                         Default: 4
%   brainSegmentation - Logical volume of existing brain segmentation,
%                       same size as input volume! Alternative to using
%                       built in PCNN segmentation of the brain
%   maxSkullThickness - Define thickness of skull (in mm). Default:
%                           0.7
%   evalFile - .mat file which contains ground truth of skull
%              'skulltruth', scale factor 'factor' in relation to scale 
%               of input slices and - optionally - ground truth of
%               brain 'braintruth'
%   fScoreBeta - Indicates how 'many times' recall is more important
%                than precision. Default: 2
%   pcnnRepetitions - Repetitions for averaging of brain segmentation.
%                     Default: 1
%   outFile - Specify file for saving the results

%% ASCII art

disp("==================================================================")
disp(" Automatic Segmentation of the Mouse Skull (AuSoMS) in MR images  ")
disp("==================================================================")
disp("                        _____         _____                       ")
disp("                       /     \_______/     \                      ")
disp("                      |  (_             _)  |                     ")
disp("                       \__               __/                      ")
disp("                          \  (  ) (  )  /                         ")
disp("                           \           /                          ")
disp("                           ->(  (_)  )<-                          ")
disp("                               |_|_|                              ")
disp("                                                                  ")

%% Check if default options are needed
if ~isfield(options, 'resizeFactor')
    options.resizeFactor = 1;
end
if ~isfield(options, 'cutoff')
    options.cutoff = 0;
end
if ~isfield(options, 'nSegmentsForeground')
    options.nSegmentsForeground = 4;
end
if ~isfield(options, 'maxSkullThickness')
    options.maxSkullThickness = 0.7;
end
if ~isfield(options, 'fScoreBeta')
    options.fScoreBeta = 2;
end
if ~isfield(options, 'pcnnRepetitions')
    options.pcnnRepetitions = 1;
end

%% Read DICOM images
[mriStack, info] = readDicom(dir);

%% Trim sequence
cutoff = options.cutoff;
mriStack = mriStack(:,:,cutoff+1:size(mriStack, 3)-cutoff);

referenceInfo = info{1, 1}; % The metadata of the first image

%% Scale images
mriStackUpscale = imresize(mriStack, options.resizeFactor, 'lanczos3');
pixelSpacing = (1/options.resizeFactor) * referenceInfo.PixelSpacing; % Adjust pixel spacing accordingly
sliceThickness = referenceInfo.SliceThickness; % Slice thickness stays the same
if ~isfield(options, 'SEBg')
    options.SEBg = strel("disk", round(0.3905 / pixelSpacing(1)), 8);
end
if ~isfield(options, 'SEBgClose')
    options.SEBgClose = strel("disk", round(0.1562 / pixelSpacing(1)), 8);
end
if ~isfield(options, 'sizeSEBrain')
    options.sizeSEBrain = round(0.3905 / pixelSpacing(1));
end

%% Preprocessing
segmentationTimer = tic;
disp("Preprocessing slices...")
preprocessedStack = preprocessSliceBased(mriStackUpscale, pixelSpacing);
%% Segmentation
foregroundSegmentation = getForegroundSegmentation(preprocessedStack, ...
    options.SEBg, options.SEBgClose, options.nSegmentsForeground);
if ~isfield(options, 'brainSegmentation')
    % No brain segmentation is given
    disp("Segmenting brain using 3D-PCNN...")
    brain = getBrainSegmentation(preprocessedStack, pixelSpacing, sliceThickness, options.sizeSEBrain, options.pcnnRepetitions);
else
    % Brain segmentation is given
    disp("Using brain segmentation given as parameter...")
    if ~isequal(size(options.brainSegmentation, 3), size(foregroundSegmentation))
        disp("Applying cutoff to brain segmentation...")
        brain = options.brainSegmentation(:,:,cutoff+1:size(options.brainSegmentation, 3)-cutoff);
    end
    if ~isequal(size(brain), size(foregroundSegmentation)) % Size of the volumes must agree
        error("Please make sure the seperate slices of the brain segmentation are of the same size as the input MR slices.")
    end
    brain = imresize(brain, options.resizeFactor, 'lanczos3'); % Resize brain segmentation according to resize factor
end
roi = combineSegmentations(foregroundSegmentation, brain, pixelSpacing, options.maxSkullThickness);
roiWithoutFragments = assessSegments(roi, pixelSpacing);
segmentation = removeOutliers(roiWithoutFragments);

%% Postprocessing the segmentation
disp("Postprocessing...")
skull = postprocess(segmentation, pixelSpacing);
disp("Segmentation successful!")
elapsedTime = toc(segmentationTimer);
disp(strcat("Elapsed time (s): ", num2str(elapsedTime)));

%% Evaluation and saving
eval = struct();
if isfield(options, 'evalFile')
    disp("Evaluating results...")
    load(options.evalFile, 'skulltruth', 'factor', 'braintruth');
    if (exist('skulltruth', 'var') && exist('factor', 'var'))
        skulltruth = skulltruth(:,:,cutoff+1:size(skulltruth, 3)-cutoff);
        skulltruth = cast(imresize(skulltruth, options.resizeFactor / factor, 'lanczos3'), 'logical');
        if exist('braintruth', 'var')
            braintruth = braintruth(:,:,cutoff+1:size(braintruth, 3)-cutoff);
            braintruth = cast(imresize(braintruth, options.resizeFactor / factor, 'lanczos3'), 'logical');
            eval = evaluate(skull, skulltruth, options.fScoreBeta, brain, braintruth);
        else
            disp("Evaluating skull only...")
            eval = evaluate(skull, skulltruth, options.fScoreBeta);
        end
    else
        warning("Evaluation failed. Please make sure the .mat file contains at least the variables 'skulltruth' and 'factor'.")
    end
end

if isfield(options, 'outFile')
    save(outFile, 'skull', 'brain', 'eval');
    disp(strcat("Saved result in ", outFile));
else
    warning("Result was not saved");
end
end

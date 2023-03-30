%% Files
dicomDir = uigetdir();
[name, folder] = uigetfile();
evalFile = fullfile(folder, name);

%% Read data
[stack, inf] = readDicom(dicomDir);
referenceInfo = inf{1, 1};
ps = referenceInfo.PixelSpacing;
st = referenceInfo.SliceThickness;
load(evalFile);

%% Trim stacks
cutoff = 2;
stack = stack(:,:,cutoff+1:size(stack, 3)-cutoff);
braintruthCut = braintruth(:,:,cutoff+1:size(braintruth, 3)-cutoff);
skulltruth = skulltruth(:,:,cutoff+1:size(skulltruth, 3)-cutoff);

%% Define image size
assert(size(stack, 1) == size(stack, 2)) % Image has square base
s = 256;
rF = s / size(stack, 1);
trF = rF / factor;
v = imresize(stack, rF, 'lanczos3'); % Rescale stack
btruthCut = imresize(braintruthCut, trF, 'lanczos3'); % Rescale brain ground truth
struth = imresize(skulltruth, trF, 'lanczos3'); % Rescale skull ground truth
ps = (1/rF) * ps; % Rescale pixel spacing

%% Brain segmentation (non-deterministic)
seSizes = [4 5 6 7 10]; % For 256
beta = 2;
pre = preprocessSliceBased(v, ps);
reps = 1000;
brainresults = cell(1, reps * size(seSizes, 2));
eval = cell(4, 1);
idx = 1;
for i = 1:size(seSizes, 2)
    for iter = 0:reps-1
        bseg = getBrainSegmentation(pre, ps, st, seSizes(i), 1);
        dice = getDice(bseg, btruthCut);
        fscore = getFScore(bseg, btruthCut, beta);
        mcc = getMCC(bseg, btruthCut);
        eval{1, 1} = dice;
        eval{2, 1} = fscore;
        eval{3, 1} = mcc;
        eval{4, 1} = seSizes(i);
        brainresults{1, idx} = eval;
        idx = idx + 1;
    end
end
save('pcnn_results', 'brainresults');

%% Skull segmentation (deterministic)
seBg = (1:1:18);
seBgClose = (1:1:18);
resizes = [0.5 1];
maxSt = 0.7;
nF = 4;
beta = 2;
skullresults = cell(size(seBg, 2) * size(seBgClose, 2) * size(resizes, 2), 6);
options = struct();
options.evalFile = evalFile;
idx = 1;
btResize = imresize(braintruth, 1/factor, 'lanczos3');
options.brainSegmentation = btResize;
for k = 1:size(resizes, 2)
    options.resizeFactor = resizes(k);
    for i = 1:size(seBg, 2)
        for j = 1:size(seBgClose, 2)
            options.SEBg = strel('disk', round(seBg(i) .* resizes(k)));
            options.SEBgClose = strel('disk', round(seBgClose(j) .* resizes(k)));
            [~,~,eva] = segmentVolume(dicomDir, options);
            % Evaluate
            skullresults{idx, 1} = seBg(i);
            skullresults{idx, 2} = seBgClose(j);
            skullresults{idx, 3} = resizes(k);
            skullresults{idx, 4} = maxSt;
            skullresults{idx, 5} = nF;
            skullresults{idx, 6} = eva.skull;
            idx = idx + 1;
        end
    end
end
resultTable = cell2table(skullresults, 'VariableNames', ...
    {'SE Bg [px]', 'SE Bg Close [px]', 'Resize factor', ...
    'Max. Skull Thickness [mm]', 'Segments foreground', 'Result'});

save('skull_results', 'skullresults');

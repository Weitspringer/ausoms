%% Load data
[stack, ~] = readDicom;
[name, folder] = uigetfile;
evalFile = fullfile(folder, name);
load(evalFile);
clear name folder evalFile
slicesPerRow = 4;
slicesPerColumn = ceil(size(stack, 3) / slicesPerRow);

%% Prepare ground truth data
skulltruth = imresize(skulltruth, 1/factor, 'lanczos3');
braintruth = imresize(braintruth, 1/factor, 'lanczos3');
clear factor

%% Print histograms of ground truth sections and the sections itself
figure;
for i = 1:size(stack, 3)
    I = stack(:,:,i);
    Itruth = skulltruth(:,:,i);
    cutted = I .* Itruth;
    cutted = cutted ./ max(cutted(:)); % Normalize gray values to [0, 1]
    subplot(slicesPerColumn, slicesPerRow, i);
    [counts, bins] = imhist(cutted(Itruth));
    imhist(cutted(Itruth));
    xlabel("Intensity");
    ylabel("Frequency")
    if max(counts) ~= 0
        set(gca, 'YLim', [0 max(counts)]);
    end
    title(strcat("Slice ", num2str(i)));
end
figure;
for i = 1:size(stack, 3)
    I = stack(:,:,i);
    Itruth = skulltruth(:,:,i);
    cutted = I .* Itruth;
    cutted = cutted ./ max(cutted(:)); % Normalize gray values to [0, 1]
    subplot_tight(slicesPerColumn, slicesPerRow, i, [0.04, 0.0001]);
    imshow(cutted);
    title(strcat("Slice ", num2str(i)));
end

%% Bias Field (MICO)
iterNum = 10;
N_region=3;  q=1;
figure;
for i = 1:size(stack, 3)
    Img = stack(:,:,i);
    A=255;
    Img_original = Img;
    [nrow,ncol] = size(Img);n = nrow*ncol;
    %load ROI
    ROI = (Img>5); ROI = double(ROI);
    
    tic
    
    Bas=getBasisOrder3(nrow,ncol);
    N_bas=size(Bas,3);
    for ii=1:N_bas
        ImgG{ii} = Img.*Bas(:,:,ii).*ROI;
        for jj=ii:N_bas
            GGT{ii,jj} = Bas(:,:,ii).*Bas(:,:,jj).*ROI;
            GGT{jj,ii} = GGT{ii,jj} ;
        end
    end
    
    
    energy_MICO = zeros(3,iterNum);
    
    b=ones(size(Img));
    for ini_num = 1:1
        C=rand(3,1);
        C=C*A;
        M=rand(nrow,ncol,3);
        a=sum(M,3);
        for k = 1 : N_region
            M(:,:,k)=M(:,:,k)./a;
        end
        
        [e_max,N_max] = max(M,[], 3);
        for kk=1:size(M,3)
            M(:,:,kk) = (N_max == kk);
        end
        
        M_old = M; chg=10000;
        energy_MICO(ini_num,1) = get_energy(Img,b,C,M,ROI,q);
        
        
        for n = 2:iterNum
            
            [M, b, C]=  MICO(Img,q,ROI,M,C,b,Bas,GGT,ImgG,1, 1);
            energy_MICO(ini_num,n) = get_energy(Img,b,C,M,ROI,q);
            
            if(mod(n,1) == 0)
                PC=zeros(size(Img));
                for k = 1 : N_region
                    PC=PC+C(k)*M(:,:,k);
                end
                img_bc = Img./b;  % bias field corrected image
            end
        end
    end
    
    [M,C]=sortMemC(M,C);
    seg=zeros(size(Img));
    for k = 1 : N_region
        seg=seg+k*M(:,:,k);   % label the k-th region 
    end
    subplot_tight(slicesPerColumn, slicesPerRow, i, [0.04, 0.0001]);
    imshow(b);
    title(strcat("Slice ", num2str(i)));
end

%% Visualize brain segmentation results
spacing = 1000;

medians = zeros(2, 5);
fails = 0;
for i = 0:size(brainresultsthree, 2)/spacing - 1
    results3 = brainresultsthree(1, i*spacing + 1:(i+1)*spacing);
    results4 = brainresultsfour(1, i*spacing + 1:(i+1)*spacing);
    resultsZI = brainresultsZI(1, i*spacing + 1:(i+1)*spacing);
    resultsME = brainresultsME(1, i*spacing + 1:(i+1)*spacing);
    
    dice = zeros(1,4*spacing);
    mcc = zeros(1,4*spacing);
    
    seSize = results3{1, 1}{4, 1};
    mouses = [repmat(('Mouse 03'), 1000, 1); repmat(('Mouse 04'), 1000, 1); repmat(('Mouse ZI'), 1000, 1); repmat(('Mouse ME'), 1000, 1)];
    for j = 1:size(results3, 2)
        if ~isa(results3{1, j}, 'cell')
            dice(j) = NaN;
            mcc(j) = NaN;
            fails = fails + 1;
        else
            dice(j) = results3{1, j}{1, 1};
            mcc(j) = results3{1, j}{3, 1};
        end
        if ~isa(results4{1, j}, 'cell')
            dice(spacing+j) = NaN;
            mcc(spacing+j) = NaN;
            fails = fails + 1;
        else
            dice(spacing+j) = results4{1, j}{1, 1};
            mcc(spacing+j) = results4{1, j}{3, 1};
        end
        if ~isa(resultsZI{1, j}, 'cell')
            dice(2*spacing+j) = NaN;
            mcc(2*spacing+j) = NaN;
            fails = fails + 1;
        else
            dice(2*spacing+j) = resultsZI{1, j}{1, 1};
            mcc(2*spacing+j) = resultsZI{1, j}{3, 1};
        end
        if ~isa(resultsME{1, j}, 'cell')
            dice(3*spacing+j) = NaN;
            mcc(3*spacing+j) = NaN;
            fails = fails + 1;
        else
            dice(3*spacing+j) = resultsME{1, j}{1, 1};
            mcc(3*spacing+j) = resultsME{1, j}{3, 1};
        end
    end
    figure;
    subplot(1,2,1)
    boxplot(dice, mouses);
    title(strcat('Brain SE size = ', 32, num2str(seSize)));
    ylabel('Dice')
    subplot(1,2,2);
    boxplot(mcc, mouses);
    title(strcat('Brain SE size = ', 32, num2str(seSize)));
    ylabel('MCC')
    
    medians(1,i+1) = median(dice, 'omitnan');
    medians(2,i+1) = median(mcc, 'omitnan');
end
figure;
str={'Median Dice score'; 'Median MCC score'};
n = size(medians, 2);
cmap = colormap(summer(5));
b=bar(medians, 'FaceColor','flat');
for k = 1:n
    d = b(k);
    d.CData(1,:) = cmap(k,:);
    d.CData(2,:) = cmap(k,:);
end
hAx=gca;
hAx.XTickLabel=str;
grid on
l = cell(1,n);
l{1}='4'; l{2}='5'; l{3}='6'; l{4}='7'; l{5}='10';   
legend(b,l);
title("Median scores for the brain segmentation on all 4 sequences for different SE sizes")

%% Skull parameters
diceWhole = zeros(1, size(skullresults, 1));
diceUpper = zeros(1, size(skullresults, 1));
mccWhole = zeros(1, size(skullresults, 1));
mccUpper = zeros(1, size(skullresults, 1));
bg = zeros(1, size(skullresults, 1));
bgClose = zeros(1, size(skullresults, 1));

for i = 1:size(skullresults, 1)
    result = skullresults{i, 6};
    diceWhole(i) = result.whole.dice;
    diceUpper(i) = result.upperHalf.dice;
    mccWhole(i) = result.whole.mcc;
    mccUpper(i) = result.upperHalf.mcc;
    bg(i) = skullresults{i, 1};
    bgClose(i) = skullresults{i, 2};
end

figure;

% Dice
ZW = reshape(diceWhole(1:324), 18, 18);
[ZWmax, IW] = max(ZW(:));
[ZWmaxRow,ZWmaxCol] = ind2sub(size(ZW), IW);
ZU = reshape(diceUpper(1:324), 18, 18);
[ZUmax, IU] = max(ZU(:));
[ZUmaxRow,ZUmaxCol] = ind2sub(size(ZU), IU);
[X, Y] = meshgrid(1:18, 1:18);

subplot(2, 2, 1);
s1 = surf(X, Y, ZW, 'FaceAlpha',0.8);
s1.EdgeColor = 'none';
colorbar;
zlabel('Dice score')
xlabel('SE Size Bg')
ylabel('SE Size Close')

title('Dice score for whole sequence')
subplot(2, 2, 3);
s2 = surf(X, Y, ZU, 'FaceAlpha',0.8);
s2.EdgeColor = 'none';
colorbar;
xlabel('SE Size Bg')
ylabel('SE Size Close')
zlabel('Dice score')
title('Dice score for upper half of sequence')

% MCC
ZWM = reshape(mccWhole(1:324), 18, 18);
[ZWMmax, IWM] = max(ZWM(:));
[ZWMmaxRow,ZWMmaxCol] = ind2sub(size(ZWM), IWM);
ZUM = reshape(mccUpper(1:324), 18, 18);
[ZUMmax, IUM] = max(ZUM(:));
[ZUMmaxRow,ZUMmaxCol] = ind2sub(size(ZUM), IUM);

subplot(2, 2, 2);
s5 = surf(X, Y, ZWM, 'FaceAlpha',0.8);
s5.EdgeColor = 'none';
colorbar;
zlabel('MCC score')
xlabel('SE Size Bg')
ylabel('SE Size Close')

title('MCC score for whole sequence')
subplot(2, 2, 4);
s6 = surf(X, Y, ZUM, 'FaceAlpha',0.8);
s6.EdgeColor = 'none';
colorbar;
xlabel('SE Size Bg')
ylabel('SE Size Close')
zlabel('MCC score')
title('MCC score for upper half of sequence')

% Difference
ZD = abs(ZU - ZW);
ZDM = abs(ZUM - ZWM);

figure;
subplot(1,2,1);
s7 = surf(X, Y, ZD, 'FaceAlpha',0.8);
s7.EdgeColor = 'none';
colorbar;
zlabel('Dice difference')
xlabel('SE Size Bg')
ylabel('SE Size Close')

subplot(1,2,2);
s9 = surf(X, Y, ZDM, 'FaceAlpha',0.8);
s9.EdgeColor = 'none';
colorbar;
zlabel('MCC difference')
xlabel('SE Size Bg')
ylabel('SE Size Close')

%% Evaluation results
iterations = 1000;
diceUpper = zeros(1,4*iterations);
diceBrain = zeros(1,4*iterations);
diceWhole = zeros(1,4*iterations);
mccUpper = zeros(1,4*iterations);
mccBrain = zeros(1,4*iterations);
mccWhole = zeros(1,4*iterations);

mouses = [repmat(('Mouse 03'), 1000, 1); repmat(('Mouse 04'), 1000, 1); repmat(('Mouse ZI'), 1000, 1); repmat(('Mouse ME'), 1000, 1)];

for j = 1:iterations
    eval1 = m1results{1,j}; % 03
    eval2 = m2results{1,j}; % 04
    eval3 = m3results{1,j}; % ZI
    eval4 = m4results{1,j}; % ME
    
    diceUpper(j) = eval1.skull.upperHalf.dice;
    diceUpper(iterations+j) = eval2.skull.upperHalf.dice;
    diceUpper(2*iterations+j) = eval3.skull.upperHalf.dice;
    diceUpper(3*iterations+j) = eval4.skull.upperHalf.dice;
    
    diceBrain(j) = eval1.brain.dice;
    diceBrain(iterations+j) = eval2.brain.dice;
    diceBrain(2*iterations+j) = eval3.brain.dice;
    diceBrain(3*iterations+j) = eval4.brain.dice;
    
    diceWhole(j) = eval1.skull.whole.dice;
    diceWhole(iterations+j) = eval2.skull.whole.dice;
    diceWhole(2*iterations+j) = eval3.skull.whole.dice;
    diceWhole(3*iterations+j) = eval4.skull.whole.dice;
    
    
    mccUpper(j) = eval1.skull.upperHalf.mcc;
    mccUpper(iterations+j) = eval2.skull.upperHalf.mcc;
    mccUpper(2*iterations+j) = eval3.skull.upperHalf.mcc;
    mccUpper(3*iterations+j) = eval4.skull.upperHalf.mcc;
    
    mccBrain(j) = eval1.brain.mcc;
    mccBrain(iterations+j) = eval2.brain.mcc;
    mccBrain(2*iterations+j) = eval3.brain.mcc;
    mccBrain(3*iterations+j) = eval4.brain.mcc;
    
    mccWhole(j) = eval1.skull.whole.mcc;
    mccWhole(iterations+j) = eval2.skull.whole.mcc;
    mccWhole(2*iterations+j) = eval3.skull.whole.mcc;
    mccWhole(3*iterations+j) = eval4.skull.whole.mcc;
end

figure;

subplot(1,2,1);
boxplot(diceUpper, mouses);
title('Dice score calculated on upper half of the sequence');
ylabel('Dice')

subplot(1,2,2);
boxplot(mccUpper, mouses);
title('MCC score calculated on upper half of the sequence');
ylabel('MCC')


figure;

subplot(1,2,1);
boxplot(diceBrain, mouses);
title('Dice score for brain segmentation');
ylabel('Dice')

subplot(1,2,2);
boxplot(mccBrain, mouses);
title('MCC score for brain segmentation');
ylabel('MCC')


figure;

subplot(1,2,1);
boxplot(diceWhole, mouses);
title('Dice score calculated on whole sequence');
ylabel('Dice')

subplot(1,2,2);
boxplot(mccWhole, mouses);
title('MCC score calculated on whole sequence');
ylabel('MCC')

meanArray = [median(diceBrain,'omitnan'), median(diceWhole,'omitnan'), median(diceUpper,'omitnan');
    median(mccBrain,'omitnan'), median(mccWhole,'omitnan'), median(mccUpper,'omitnan')];
figure;
str={'Median Dice score'; 'Median MCC score'};
n = size(meanArray, 2);
cmap = colormap(summer(n));
b=bar(meanArray, 'FaceColor','flat');
for k = 1:n
    d = b(k);
    d.CData(1,:) = cmap(k,:);
    d.CData(2,:) = cmap(k,:);
end
hAx=gca;
hAx.XTickLabel=str;
grid on
l = cell(1,n);
l{1}='Brain'; l{2}='Whole skull'; l{3}='Skullcap';   
legend(b,l);
title("Median scores of the segmentation for sequences 1-4")
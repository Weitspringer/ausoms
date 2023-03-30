%% Files
[name, folder] = uigetfile();
dataFile = fullfile(folder, name);
[name, folder] = uigetfile();
evalFile = fullfile(folder, name);

%% Options
options = struct();
options.cutoff = 2;
options.pcnnRepetitions = 1;
iterations = 1000;

%% Segment and calculate scores
results = zeros(1, iterations);
options.evalFile = evalFile;
for i = 1:iterations
    [~,~,eval] = segmentVolume(dataFile, options);
    results(i) = eval;
end

save('evaluation', 'results');
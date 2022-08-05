% coherent averaging

subjectList = {'nl-1624','nl-jg','nl-1928','nl-0055','nl-1909','nl-1908','nl-0051','nl-0034','nl-myl','nl-rta','nl-0057','nl-2274','nl-2275',...
    'nl-2277','nl-2278','nl-2279','nl-2280','nl-2281','nl-2282', 'nl-2283','nl-2285', 'nl-2287', 'nl-2300','nl-1539','nl-2113','nl-2126','nl-2215',...
    'nl-2329','nl-2331','nl-2332','nl-2286'};%,'nl-2335','nl-2336','nl-2343','nl-2344','nl-2346','nl-2347','nl-2349','nl-2350','nl-2351','nl-2353', 'nl-2374',...
    %'nl-2378', 'nl-2379','nl-2381', 'nl-2383'};
fileDir = '~/data/orientationRCA/snr';
dataDir = '~/data/orientationtuning';
matfileDir = 'Exp_MATL_HCN_128_Avg';
taskList = [zeros(31,1)];%;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];

[orderedSNR, rankInd] = sort(meanSNR(1:nSubj,5), 'descend');
isOver2 = orderedSNR>2;

subjectList = subjectList(rankInd);
subjectList = subjectList(isOver2);
taskList = task(rankInd);
taskList = taskList(isOver2);

nSubj= length(subjectList);

target_frequencies = [8, 2];

output_cos = []; output_sin = []; output_amp=[]; output_side_cos=[]; output_side_sin=[];
for s = 1:length(subjectList)
    thisFolder = dir(fullfile(dataDir, [subjectList{s},'*']));
    thisDir = fullfile(dataDir, thisFolder.name, matfileDir);
        
    [output_cos{s}, output_sin{s}, output_amp{s}, output_side_cos{s}, output_side_sin{s}] =...
        getPhases_occipital(thisDir, target_frequencies, taskList(s));
    
end
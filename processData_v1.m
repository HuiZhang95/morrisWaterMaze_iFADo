%% check the data 22July2021
% checkBehData_v2.m
%%
clear all
saveDataPath = 'F:\project\navigationParadigm\MorrisMaze\analyzeBehData\';
load([saveDataPath,'behInfo_ifado_v2.mat']);
addpath('F:\project\navigationParadigm\MorrisMaze\matlabCode\privateFun');
nameIdx = 2; % 1 is age; 2 is genotype
matchIdealPath_flag = 0; % 1 match; 0 no match; to match the ideal path length between carrier and non-carriers
matlabIdealPath_trialNo = 10;
subjNoAll = length(behInfo);
[ageInfo,geneInfo,genderInfo,matchIdealPath_vectorAll,subjIdxPool,condiName] = morrisWaterMaze_extract_subInfo_04Nov2021(behInfo,nameIdx,matchIdealPath_flag,matlabIdealPath_trialNo);  

%% plot demograph data
% plot gene information
[carrier_idx,nonCarrier_idx,e23Carrier_idx] = plotDemographicData_morrisWaterMaze(ageInfo,geneInfo,genderInfo);

%% check the behavioral performance for lure and teleported conditions
diffBetweenLureAndTeleport(genderInfo,matchIdealPath_vectorAll,behInfo)
%% to plot the distribution of target locations for carriers and noncarriers
subjPool = {'carrier_idx','nonCarrier_idx'};%,'e23Carrier_idx'};
target_and_StartSearchLoc_distribution_morris(subjPool,subjIdxPool,behInfo,carrier_idx,nonCarrier_idx);

%% to plot the number of trials finished by carriers and nonCarriers
[data_carrier,data_non] = trialNumbers(subjPool,subjIdxPool,behInfo,carrier_idx,nonCarrier_idx);
% carrier_idx(data_carrier<131) = [];
% nonCarrier_idx(data_non<120) = [];

%% to investigate if there is any difference of behavioral performance 
% between carriers and non-carriers
fieldNames = {'targetLoc_to_center','dropLoc_to_center','startSearchLoc_to_center',...
    'ideal_distance_lenth','actural_distance_lenth','dropError','distanceError',...
    'distance_underEstimateRatio','angleError','actural_path_lenth','elapsedSearchTime',...
    'picMemory','RT_itemMemory'};
% fieldNames = {'targetLoc_to_landmark','dropLoc_to_landmark','startSearchLoc_to_landmark'};
subjData = behaviralDiff_carrierVSnonCarrier(fieldNames,matchIdealPath_vectorAll,behInfo,subjIdxPool,ageInfo,geneInfo);
behaviralDiff_carrierVSnonCarrierVSe23Carrier(fieldNames,matchIdealPath_vectorAll,behInfo,subjIdxPool,ageInfo,geneInfo);
behaviralDiff_carrierVSnonCarrier_withTrial(fieldNames,matchIdealPath_vectorAll,behInfo,subjIdxPool,ageInfo,geneInfo);
%% correlation between performance and age
age_corr_performance(behInfo,fieldNames,subjIdxPool,ageInfo,subjData);

%% correlate between performance and the distance between target and the nearest landmark/center
fieldNames = {'dropError','distanceError','distance_underEstimateRatio',...
    'angleError','picMemory','RT_itemMemory'};
fieldNames_design = {'targetLoc_to_center','targetLoc_to_landmark'};
performance_corr_targetLandmarkCenterDis(fieldNames,fieldNames_design,matchIdealPath_vectorAll,behInfo,subjIdxPool)


%% to investigate deviation of the navigation trajactory; the navigation pattern 
fieldNames = {'dis_to_ideal','dis_to_actural','inboundPath_to_center'};
% fieldNames = {'dropError','distanceError','distance_underEstimateRatio',...
%     'angleError','picMemory','RT_itemMemory'};
navigationPattern_carrierVSnonCarrier(fieldNames, behInfo,matchIdealPath_vectorAll,subjIdxPool)

%% correlate between performance and navigation pattern
fieldNames = {'dropError','distanceError','distance_underEstimateRatio',...
    'angleError','picMemory','RT_itemMemory'};
fieldNames_pattern = {'dis_to_ideal','dis_to_actural','inboundPath_to_center','ideal_distance_lenth'};
performance_corr_naviPattern(subjIdxPool,fieldNames,fieldNames_pattern,behInfo,matchIdealPath_vectorAll)


%% check the propotion of non-moving period and their difference between 
% carriers and noncarriers
nonMoving_carrierVSnonCarrier(matchIdealPath_vectorAll,behInfo,subjIdxPool);

%% plot the in-bound trajactory of carriers and non-carriers
clear all
saveDataPath = 'F:\project\navigationParadigm\MorrisMaze\analyzeBehData\';
load([saveDataPath,'behInfo_ifado_v2.mat']);
plotInbundPath(behInfo)

%% multilevel regression analysis
fieldNames = {'dropError','distanceError','distance_underEstimateRatio',...
    'angleError','picMemory','RT_itemMemory','targetLoc_to_center',...
    'targetLoc_to_landmark','ideal_distance_lenth','inboundPath_to_center_per',...
    'dis_to_actural','actural_path_lenth','elapsedSearchTime','nonMovingRatio'};
multilevelRegression_v3(behInfo,matchIdealPath_vectorAll,fieldNames);
% multilevelRegression(behInfo,matchIdealPath_vectorAll,fieldNames);



%% plot data from R






































































































%% to investigate the navigation pattern (near boundary or not)
fieldNames = {'inboundPath_to_center'};
subjData = NaN(subjNoAll,3,16);
subjData2 = NaN(subjNoAll,3,16);
for subjNo = 1:subjNoAll
    matchIdealPath_vector = matchIdealPath_vectorAll{subjNo};
    for fieldNo = 1:length(fieldNames)
        trialNoAll = length(behInfo(subjNo).trialInfo_furtherProcess);
        if ~isempty(behInfo(subjNo).trialInfo_furtherProcess)
            trialdata = zeros(trialNoAll,16);
            trialdata_idx = zeros(trialNoAll,1);
            objLoc = NaN(trialNoAll,1);
            for trialNo = 1:trialNoAll
                linearSpeed = behInfo(subjNo).trialInfo_furtherProcess(trialNo).linearSpeed;
                linearSpeed(1) = 0;
                if length(behInfo(subjNo).trialInfo_furtherProcess(trialNo).dis_to_ideal)>1
                    trialdata_idx(trialNo) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
                    tmp_objLoc = behInfo(subjNo).trialInfo_furtherProcess(trialNo).targetLoc_to_center;
                    varName = eval(strcat('behInfo(subjNo).trialInfo_furtherProcess(trialNo).',fieldNames{fieldNo}));
                    varName(linearSpeed==0) = [];
                    varName1 = round(varName); 
                    varName1(varName1>16) = [];
                    tmpRange = [0:1:16];
                    for i = 1:length(tmpRange)-1
                        trialdata(trialNo,i) = sum(varName1>=tmpRange(i) & varName1 < tmpRange(i+1))*100/length(varName1);
                    end
                    objLoc(trialNo) = find(tmpRange==round(tmp_objLoc));
                else
                    trialdata_idx(trialNo) = NaN;
                    trialdata(trialNo,:) = NaN;
                    objLoc(trialNo) = NaN;
                end
            end

            tmpidx = find(trialdata_idx==1 & matchIdealPath_vector(:,1)==1);
            subjData(subjNo,1,:) = nanmean(trialdata(tmpidx,:));
            tmpidx = find(trialdata_idx==0 & matchIdealPath_vector(:,2)==1);
            subjData(subjNo,2,:) = nanmean(trialdata(tmpidx,:));
            tmpidx = find((trialdata_idx==0 | trialdata_idx==1) & matchIdealPath_vector(:,3)==1);
            subjData(subjNo,3,:) = nanmean(trialdata(tmpidx,:));

            subjData2_tmp = NaN(3,trialNoAll,16);
            for locIdx = 1:16
                tmpidx = find(objLoc==locIdx & trialdata_idx==1 & matchIdealPath_vector(:,1)==1);
                if ~isempty(tmpidx)
                    tmpidx2 = abs([1:16]-locIdx)+1; tmpidx3 = unique(tmpidx2);
                    for tmpidx4 = 1:length(tmpidx3)
                        subjData2_tmp(1,tmpidx,tmpidx3(tmpidx4)) = nanmean(trialdata(tmpidx,find(tmpidx2==tmpidx3(tmpidx4))),2);
                    end
                end

                tmpidx = find(objLoc==locIdx & trialdata_idx==0 & matchIdealPath_vector(:,2)==1);
                if ~isempty(tmpidx)
                    tmpidx2 = abs([1:16]-locIdx)+1; tmpidx3 = unique(tmpidx2);
                    for tmpidx4 = 1:length(tmpidx3)
                        subjData2_tmp(2,tmpidx,tmpidx3(tmpidx4)) = nanmean(trialdata(tmpidx,find(tmpidx2==tmpidx3(tmpidx4))),2);
                    end
                end

                tmpidx = find(objLoc==locIdx & (trialdata_idx==0 | trialdata_idx==1) & matchIdealPath_vector(:,3)==1);
                if ~isempty(tmpidx)
                    tmpidx2 = abs([1:16]-locIdx)+1; tmpidx3 = unique(tmpidx2);
                    for tmpidx4 = 1:length(tmpidx3)
                        subjData2_tmp(3,tmpidx,tmpidx3(tmpidx4)) = nanmean(trialdata(tmpidx,find(tmpidx2==tmpidx3(tmpidx4))),2);
                    end
                end
            end
            subjData2(subjNo,:,:) = squeeze(nanmean(subjData2_tmp,2));
        end
    end
end

for fieldNo = 1:length(fieldNames)
    vLabel = fieldNames{fieldNo};
    meanD = zeros(3,2,size(subjData,3));
    stdD = zeros(3,2,size(subjData,3));
    pValues = zeros(3,size(subjData,3));
    tmpData_memory = cell(3,2);
    for condiNo = 1:3 % 1:lure; 2: teleport; 3: lure+teleport
        for i = 1:2
            subjIdx = subjIdxPool{i};
            tmpData = squeeze(subjData(subjIdx,condiNo,:));
            tmpData_memory{condiNo,i} = tmpData;
            meanD(condiNo,i,:) = nanmean(tmpData);
            stdD(condiNo,i,:) = nanstd(tmpData)/sqrt(length(subjIdx));
        end
        [~,p,~,~] = ttest2(tmpData_memory{condiNo,1},tmpData_memory{condiNo,2});
        pValues(condiNo,:) = squeeze(p);
    end
    
    figure;
    inputColor = {[0 0.5 0.5],[0.5 0 0.5]};
    inputXlabel = {'lure','teleport','combine'};
    inputLegend = condiName;
    for condiNo = 1:3
        subplot(3,1,condiNo);
        for i = 1:2
            tmpData_mean = squeeze(meanD(condiNo,i,:));
            tmpData_std = squeeze(stdD(condiNo,i,:));
            plot(1:length(tmpData_mean),tmpData_mean,'-','color',inputColor{i}); hold on
            errorbar(1:length(tmpData_mean),tmpData_mean,tmpData_std,'.','color',[0 0.5 0.5]); 
        end
        
        tmpP = pValues(condiNo,:) ; tmpP_idx = find(tmpP<0.05);
        if ~isempty(tmpP_idx)
            hold on;
            plot(tmpP_idx,tmpData_mean(tmpP_idx),'r*');
        end
        ylabel('% of inbound path')
        xlabel('Distance to the center of the island')
        title(inputXlabel{condiNo})
    end
    
    disRange = [0 5 10 16];
    meanD = zeros(3,2,length(disRange)-1);
    stdD = zeros(3,2,length(disRange)-1);
    pValues = cell(4,length(disRange)-1);
    tmpData_memory = cell(3,2);
    for condiNo = 1:3 % 1:lure; 2: teleport; 3: lure+teleport
        for i = 1:2
            subjIdx = subjIdxPool{i};
            tmpData1 = NaN(length(subjIdx),length(disRange)-1);
            for disRangeNo = 1:length(disRange)-1
                tmpData = squeeze(nansum(subjData(subjIdx,condiNo,disRange(disRangeNo)+1:disRange(disRangeNo+1)),3));
                tmpData1(:,disRangeNo) = tmpData;
                meanD(condiNo,i,disRangeNo) = nanmean(tmpData);
                stdD(condiNo,i,disRangeNo) = nanstd(tmpData,[],1)/sqrt(length(subjIdx));
            end
            tmpData_memory{condiNo,i} = tmpData1;
        end
        [~,p,~,~] = ttest2(tmpData_memory{condiNo,1},tmpData_memory{condiNo,2});
        p = squeeze(p);
        for disRangeNo = 1:length(disRange)-1
            p1 = num2str(p(disRangeNo));  pValues{condiNo,disRangeNo} = p1(1:4);
        end
    end
    [~,p,~,~] = ttest2(tmpData_memory{1,1}-tmpData_memory{2,1},tmpData_memory{1,2}-tmpData_memory{2,2});
    p = squeeze(p);
    for disRangeNo = 1:length(disRange)-1
        p1 = num2str(p(disRangeNo));  pValues{condiNo+1,disRangeNo} = p1(1:4);
    end
    
    addpath('H:\matlabCode\with&withoutpic\privateFunctions');
    figure;
    areaName = {'center','middle','outer'};
    for disRangeNo = 1:length(disRange)-1
        subplot(2,2,disRangeNo);
        meanD1 = squeeze(meanD(:,:,disRangeNo));
        stdD1 = squeeze(stdD(:,:,disRangeNo));
        rawDataOrNot = 2;
        inputXlabel = {'lure','teleport','combine'};
        inputColor = {[0 0.5 0.5],[0.5 0 0.5]};
        inputLegend = NaN;%condiName;
        [~,~,hb]=plotBar_color(meanD1,inputColor, inputLegend, inputXlabel, rawDataOrNot,stdD1);
        ylabel('% of inbound path')
        title({strcat(areaName{disRangeNo},' ',condiName{1},' vs.',condiName{2}) strcat(' Lure p=',pValues{1,disRangeNo},'teleport p=',pValues{2,disRangeNo}) ...
            strcat('interaction p=',pValues{4,disRangeNo},'mainEffect p=',pValues{3,disRangeNo})});
    end
end



for fieldNo = 1:length(fieldNames)
    vLabel = fieldNames{fieldNo};
    meanD = zeros(3,2,size(subjData2,3)-1);
    stdD = zeros(3,2,size(subjData2,3)-1);
    pValues = zeros(3,size(subjData2,3)-1);
    tmpData_memory = cell(3,2);
    for condiNo = 1:3 % 1:lure; 2: teleport; 3: lure+teleport
        for i = 1:2
            subjIdx = subjIdxPool{i};
            tmpData = squeeze(subjData2(subjIdx,condiNo,1:end-1));
            tmpData_memory{condiNo,i} = tmpData;
            meanD(condiNo,i,:) = nanmean(tmpData);
            stdD(condiNo,i,:) = nanstd(tmpData)/sqrt(length(subjIdx));
        end
        [~,p,~,~] = ttest2(tmpData_memory{condiNo,1},tmpData_memory{condiNo,2});
        pValues(condiNo,:) = squeeze(p);
    end
    
    figure;
    inputColor = {[0 0.5 0.5],[0.5 0 0.5]};
    inputXlabel = {'lure','teleport','combine'};
    inputLegend = condiName;
    for condiNo = 1:3
        subplot(3,1,condiNo);
        for i = 1:2
            tmpData_mean = squeeze(meanD(condiNo,i,:));
            tmpData_std = squeeze(stdD(condiNo,i,:));
            plot(1:length(tmpData_mean),tmpData_mean,'-','color',inputColor{i}); hold on
            errorbar(1:length(tmpData_mean),tmpData_mean,tmpData_std,'.','color',[0 0.5 0.5]); 
        end
        
        tmpP = pValues(condiNo,:) ; tmpP_idx = find(tmpP<0.05);
        if ~isempty(tmpP_idx)
            hold on;
            plot(tmpP_idx,tmpData_mean(tmpP_idx),'r*');
        end
        ylabel('% of inbound path')
        xlabel('Distance between target location and inbound path')
        title(inputXlabel{condiNo})
    end
end



%% to investigate the target-location to center vs. drop location to center
fieldNames1 = {'targetLoc_to_center','targetLoc_to_landmark','dropError','dropError','distanceError'};
fieldNames2 = {'dropLoc_to_center','dropLoc_to_landmark','distanceError','angleError','angleError'};
subjData = NaN(subjNoAll,3,length(fieldNames1));
subjData2 = NaN(subjNoAll,3,length(fieldNames1));
for subjNo = 1:subjNoAll
    matchIdealPath_vector = matchIdealPath_vectorAll{subjNo};
    for fieldNo = 1:length(fieldNames1)
        trialNoAll = length(behInfo(subjNo).trialInfo_furtherProcess);
        if ~isempty(behInfo(subjNo).trialInfo_furtherProcess)
            trialdata = zeros(trialNoAll,16);
            trialdata_idx = zeros(trialNoAll,1);
            targetLoc = NaN(trialNoAll,1);
            dropLoc = NaN(trialNoAll,1);
            for trialNo = 1:trialNoAll
                linearSpeed = behInfo(subjNo).trialInfo_furtherProcess(trialNo).linearSpeed;
                linearSpeed(1) = 0;
                if length(behInfo(subjNo).trialInfo_furtherProcess(trialNo).dis_to_ideal)>1
                    trialdata_idx(trialNo) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
                    tmp_targetLoc = eval(strcat('behInfo(subjNo).trialInfo_furtherProcess(trialNo).',fieldNames1{fieldNo}));
                    tmp_dropLoc = eval(strcat('behInfo(subjNo).trialInfo_furtherProcess(trialNo).',fieldNames2{fieldNo}));

                    if fieldNo == 1 | fieldNo>2
                        targetLoc(trialNo) = tmp_targetLoc;
                        dropLoc(trialNo) = tmp_dropLoc;
                    elseif fieldNo == 2
                        [tmp1, tmp2] = min(tmp_targetLoc);
                        targetLoc(trialNo) = tmp1;
                        dropLoc(trialNo) = tmp_dropLoc(tmp2);
                    end
                else
                    trialdata_idx(trialNo) = NaN;
                    trialdata(trialNo,:) = NaN;
                    targetLoc(trialNo) = NaN;
                    dropLoc(trialNo) = NaN;
                end
            end

            tmpidx = find(trialdata_idx==1 & matchIdealPath_vector(:,1)==1);
            var1 = targetLoc(tmpidx); var2 = dropLoc(tmpidx);
            subjData(subjNo,1,fieldNo) = atanh(corr(var1,var2,'type','spearman'));
            tmpidx = find(trialdata_idx==0 & matchIdealPath_vector(:,2)==1);
            var1 = targetLoc(tmpidx); var2 = dropLoc(tmpidx);
            subjData(subjNo,2,fieldNo) = atanh(corr(var1,var2,'type','spearman'));
            tmpidx = find((trialdata_idx==0 | trialdata_idx==1) & matchIdealPath_vector(:,3)==1);
            var1 = targetLoc(tmpidx); var2 = dropLoc(tmpidx);
            subjData(subjNo,3,fieldNo) = atanh(corr(var1,var2,'type','spearman'));
        end
    end
end

figure;
for fieldNo = 1:length(fieldNames1)
    meanD = zeros(3,2);
    stdD = zeros(3,2);
    tmpData_memory = cell(3,2);
    pAll = zeros(3,2);
    pValues = cell(1,4);
    for condiNo = 1:3 % 1:lure; 2: teleport; 3: lure+teleport
        for i = 1:2
            subjIdx = subjIdxPool{i};
            tmpData = squeeze(subjData(subjIdx,condiNo,fieldNo));
            [~,p,~,~] = ttest(tmpData);
            pAll(condiNo,i) = p;
            tmpData_memory{condiNo,i} = tmpData;
            meanD(condiNo,i) = nanmean(tmpData);
            stdD(condiNo,i) = nanstd(tmpData)/sqrt(length(subjIdx));
        end
        [~,p,~,~] = ttest2(tmpData_memory{condiNo,1},tmpData_memory{condiNo,2});
        p = num2str(p); pValues{condiNo} = p(1:4);
    end
    [~,p,~,~] = ttest2(tmpData_memory{1,1}-tmpData_memory{2,1},tmpData_memory{1,2}-tmpData_memory{2,2});
    p = num2str(p); pValues{condiNo+1} = p(1:4);
    addpath('H:\matlabCode\with&withoutpic\privateFunctions');
    subplot(2,3,fieldNo);%subplot(2,1,count);%
    rawDataOrNot = 2;
    inputXlabel = {'lure','teleport','combine'};
    inputColor = {[0 0.5 0.5],[0.5 0 0.5]};
    inputLegend = NaN;%condiName;
    [~,~,hb]=plotBar_color(meanD,inputColor, inputLegend, inputXlabel, rawDataOrNot,stdD);
    for ib = 1:numel(hb)
        xData = hb(ib).XData+hb(ib).XOffset;
        for condiNo = 1:size(tmpData_memory,1)
            [~,p,~,~] = ttest(tmpData_memory{condiNo,ib});
            if p<0.05
                text(xData(condiNo),meanD(condiNo,ib)+stdD(condiNo,ib),...
                    '*','color','r','fontSize',20);
            end
        end
    end
%         pbaspect([1.5 1 1]);
    title({strcat(fieldNames1{fieldNo},' corr. with ', fieldNames2{fieldNo},' ',condiName{1},' vs.',condiName{2}) ...
        strcat(' Lure p=',pValues{1},'teleport p=',pValues{2}) ...
        strcat('interaction p=',pValues{4},'mainEffect p=',pValues{3})});
end















































































































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% to investigate if the behavioral performance
% is correlated with the experimental setups and how they related with
% APOE4
fieldNames1 = {'elapsedSearchTime'};
fieldNames2 = {'actural_distance_lenth','dropError','distanceError',...
    'distance_underEstimateRatio','angleError','actural_path_lenth',...
    'targetLoc_to_center','dropLoc_to_center','startSearchLoc_to_center',...
    'ideal_distance_lenth','targetLoc_to_landmark','startSearchLoc_to_landmark'}; %'dropLoc_to_landmark'


% fieldNames1 = {'actural_distance_lenth','dropError','distanceError',...
%     'distance_underEstimateRatio','angleError','actural_path_lenth','elapsedSearchTime'};
% fieldNames2 = {'targetLoc_to_center','dropLoc_to_center','startSearchLoc_to_center',...
%     'ideal_distance_lenth','targetLoc_to_landmark','startSearchLoc_to_landmark','elapsedSearchTime'}; %'dropLoc_to_landmark'
subjData = NaN(subjNoAll,3,length(fieldNames1),length(fieldNames2));
for subjNo = 1:subjNoAll
    matchIdealPath_vector = matchIdealPath_vectorAll{subjNo};    
    for fieldNo1 = 1:length(fieldNames1)
        for fieldNo2 = 1:length(fieldNames2)
            trialNoAll = length(behInfo(subjNo).trialInfo_furtherProcess);
            if ~isempty(behInfo(subjNo).trialInfo_furtherProcess)
                trialdata = zeros(trialNoAll,2,3);
                for trialNo = 1:trialNoAll
                    if length(behInfo(subjNo).trialInfo_furtherProcess(trialNo).dis_to_ideal)>1
                        trialdata(trialNo,1) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
                        if strfind(fieldNames2{fieldNo2},'landmark')
                            varName = eval(strcat('behInfo(subjNo).trialInfo_furtherProcess(trialNo).',fieldNames1{fieldNo1}));
                            trialdata(trialNo,2) = min(varName);
                            varName = eval(strcat('behInfo(subjNo).trialInfo_furtherProcess(trialNo).',fieldNames2{fieldNo2}));
                            trialdata(trialNo,3) = min(varName);
                        else
                            varName = eval(strcat('behInfo(subjNo).trialInfo_furtherProcess(trialNo).',fieldNames1{fieldNo1}));
                            trialdata(trialNo,2) = varName;
                            varName = eval(strcat('behInfo(subjNo).trialInfo_furtherProcess(trialNo).',fieldNames2{fieldNo2}));
                            trialdata(trialNo,3) = varName;
                        end
                    else
                        trialdata(trialNo,1) = NaN;
                        trialdata(trialNo,2) = NaN;
                        trialdata(trialNo,3) = NaN;
                    end
                end

                tmpidx = find(trialdata(:,1)==1 & matchIdealPath_vector(:,1)==1);
                var1 = trialdata(tmpidx,2); var2 = trialdata(tmpidx,3);
                subjData(subjNo,1,fieldNo1,fieldNo2) = atanh(corr(var1,var2,'type','spearman'));
                tmpidx = find(trialdata(:,1)==0 & matchIdealPath_vector(:,2)==1);
                var1 = trialdata(tmpidx,2); var2 = trialdata(tmpidx,3);
                subjData(subjNo,2,fieldNo1,fieldNo2) = atanh(corr(var1,var2,'type','spearman'));
                tmpidx = find((trialdata(:,1)==0 | trialdata(:,1)==1) & matchIdealPath_vector(:,3)==1);
                var1 = trialdata(tmpidx,2); var2 = trialdata(tmpidx,3);
                subjData(subjNo,3,fieldNo1,fieldNo2) = atanh(corr(var1,var2,'type','spearman'));
            end
        end
    end
end

for fieldNo1 = 1:length(fieldNames1)
    figure;
    count = 0;
    for fieldNo2 = 1:length(fieldNames2)
        count = count+1;
        meanD = zeros(3,2);
        stdD = zeros(3,2);
        tmpData_memory = cell(3,2);
        pAll = zeros(3,2);
        pValues = cell(1,4);
        pCorrAge = zeros(3,2,2);
        for condiNo = 1:3 % 1:lure; 2: teleport; 3: lure+teleport
            for i = 1:2
                subjIdx = subjIdxPool{i};
                tmpData = squeeze(subjData(subjIdx,condiNo,fieldNo1,fieldNo2));
                [r,p] = corr(tmpData,ageInfo(subjIdx),'type','spearman');
                pCorrAge(condiNo,i,1) = r;pCorrAge(condiNo,i,2) = p;
                [~,p,~,~] = ttest(tmpData);
                pAll(condiNo,i) = p;
                tmpData_memory{condiNo,i} = tmpData;
                meanD(condiNo,i) = nanmean(tmpData);
                stdD(condiNo,i) = nanstd(tmpData)/sqrt(length(subjIdx));
            end
            [~,p,~,~] = ttest2(tmpData_memory{condiNo,1},tmpData_memory{condiNo,2});
            p = num2str(p); pValues{condiNo} = p(1:4);
        end
        [~,p,~,~] = ttest2(tmpData_memory{1,1}-tmpData_memory{2,1},tmpData_memory{1,2}-tmpData_memory{2,2});
        p = num2str(p); pValues{condiNo+1} = p(1:4);
        addpath('H:\matlabCode\with&withoutpic\privateFunctions');
        subplot(3,3,fieldNo2);%subplot(2,1,count);%
        rawDataOrNot = 2;
        inputXlabel = {'lure','teleport','combine'};
        inputColor = {[0 0.5 0.5],[0.5 0 0.5]};
        inputLegend = NaN;%condiName;
        [~,~,hb]=plotBar_color(meanD,inputColor, inputLegend, inputXlabel, rawDataOrNot,stdD);
        % differ from 0 or not
        for ib = 1:numel(hb)
            xData = hb(ib).XData+hb(ib).XOffset;
            for condiNo = 1:size(tmpData_memory,1)
                [~,p,~,~] = ttest(tmpData_memory{condiNo,ib});
                if p<0.05
                    text(xData(condiNo),meanD(condiNo,ib)+stdD(condiNo,ib)*1.5,...
                        '*','color','g','fontSize',20);
                end
            end
        end
        % correlate with age or not
        for ib = 1:numel(hb)
            xData = hb(ib).XData+hb(ib).XOffset;
            for condiNo = 1:size(tmpData_memory,1)
                if pCorrAge(condiNo,ib,2)<0.05
                    if pCorrAge(condiNo,ib,1)>0
                        text(xData(condiNo),meanD(condiNo,ib)+stdD(condiNo,ib),...
                            '*','color','r','fontSize',20);
                    elseif pCorrAge(condiNo,ib,1)<0
                        text(xData(condiNo),meanD(condiNo,ib)+stdD(condiNo,ib),...
                            '*','color','b','fontSize',20);
                    end
                end
            end
        end
        
%         pbaspect([1.5 1 1]);
        title({strcat(fieldNames2{fieldNo2},' ',condiName{1},' vs.',condiName{2}) ...
            strcat(' Lure p=',pValues{1},'teleport p=',pValues{2}) ...
            strcat('interaction p=',pValues{4},'mainEffect p=',pValues{3})});
    end
    suptitle(strcat(fieldNames1{fieldNo1},' corr. with '))
end
    
%% previous check before 22July2021
subjNoAll = length(behInfo);
ageInfo = NaN(subjNoAll,1);
geneInfo = NaN(subjNoAll,1);
genderInfo = {behInfo(:).sex};
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    trialdata = NaN(trialNoAll,5);
    
    ageInfo(subjNo) = behInfo(subjNo).age;
    if strcmp(behInfo(subjNo).geneType,'e3/e4')
        geneInfo(subjNo) = 1;
    elseif strcmp(behInfo(subjNo).geneType,'e3/e3')
        geneInfo(subjNo) = 2;
    elseif strcmp(behInfo(subjNo).geneType,'e4/e4')
        geneInfo(subjNo) = 3;
    elseif strcmp(behInfo(subjNo).geneType,'e2/e4')
        geneInfo(subjNo) = 4;
    elseif strcmp(behInfo(subjNo).geneType,'e2/e3')
        geneInfo(subjNo) = 5;
    end
    trialdata = zeros(trialNoAll,4);

    for trialNo = 1:trialNoAll
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;

        pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
        dropError = pdist([targetLoc;dropLoc],'euclidean');
        travelLenth = pdist([startSearchLoc;dropLoc],'euclidean');
        trialdata(trialNo,1) = travelLenth;
        trialdata(trialNo,2) = dropError;%dropError;
        trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
        trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
        trialdata(trialNo,5) = (behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ...
            behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime)*dropError/pathLenth;
    end
    
    tmpidx = find(trialdata(:,4)==1);
    tmpidx1 = find(trialdata(tmpidx,3)==1);
    subjData(subjNo,1) = nanmean(trialdata(tmpidx,2));
    subjData(subjNo,2) = length(tmpidx1)/length(tmpidx);
    tmpidx = find(trialdata(:,4)==0);
    tmpidx1 = find(trialdata(tmpidx,3)==1);
    subjData(subjNo,3) = nanmean(trialdata(tmpidx,2));
    subjData(subjNo,4) = length(tmpidx1)/length(tmpidx);
end

% check the age inormation between carrier and non-carriers
carrierAge = ageInfo(geneInfo==1);
nonCarrAge = ageInfo(geneInfo==2);
[~,p,~,stats]= ttest2(carrierAge,nonCarrAge);
figure;
bar(1:2,[nanmean(carrierAge) nanmean(nonCarrAge)],0.4,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5); hold on
errorbar(1:2,[nanmean(carrierAge) nanmean(nonCarrAge)],[nanstd(carrierAge)/sqrt(length(carrierAge)) nanstd(nonCarrAge)/sqrt(length(nonCarrAge))],'.','LineWidth',3, 'color',[0 .9 .9]); hold on
plot(ones(size(carrierAge)),carrierAge,'ro','lineWidth',3); hold on
plot(2*ones(size(nonCarrAge)),nonCarrAge,'ro','lineWidth',3);
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'carrier','non-carrier'});
ylim([20 75])
title(strcat('age p=',num2str(p),' carrierNo=',num2str(length(carrierAge)),...
    ' noncarrierNo=',num2str(length(nonCarrAge))));

% to investigate the difference of item- and spatial memory between young and old, carrier and
% non-carriers
for nameIdx = 1:2
    if nameIdx == 1
        ageInfo1 = ageInfo(~isnan(ageInfo));
        medianData = median(ageInfo1);
        subjIdx1 = find(ageInfo<medianData & geneInfo<3); 
        subjIdx2 = find(ageInfo>=medianData & geneInfo<3);
        condiName = {'young','old'};
    elseif nameIdx == 2
        subjIdx1 = find(geneInfo==1); % carrier
        subjIdx2 = find(geneInfo==2); % non-carrier
        condiName = {'Carrier','noncarrier'};
    end
    subjIdx3 = 1:subjNoAll;
    subjIdxPool = {'subjIdx1','subjIdx2','subjIdx3'};

    for memoryNo = 1:2
        figure;
        tmpData_memory = cell(2,2);
        for i = 1:2
            subjIdx = eval(subjIdxPool{i});
            if memoryNo == 1
                tmpData_memory{i,1} = subjData(subjIdx,2); % lure
                tmpData_memory{i,2} = subjData(subjIdx,4); % teleport
                vLabel = 'item memory accuracy';
                ylimRange = [0.5 1];
            elseif memoryNo == 2
                tmpData_memory{i,1} = subjData(subjIdx,1); % lure
                tmpData_memory{i,2} = subjData(subjIdx,3); % teleport
                vLabel = 'drop error';
                ylimRange = [0 16];
            end
            [~,p,~,stats] = ttest(tmpData_memory{i,1},tmpData_memory{i,2});
            meanD = [nanmean(tmpData_memory{i,1}) nanmean(tmpData_memory{i,2})];
            stdD = [nanstd(tmpData_memory{i,1}) nanstd(tmpData_memory{i,2})];
            subplot(1,2,i);
            bar(1:2,meanD,0.4,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5); hold on
            errorbar(1:2,meanD,stdD,'.','LineWidth',3, 'color',[0 .9 .9]); hold on
            for subjNo = 1:length(tmpData_memory{i,1})
                plot(1:2,[tmpData_memory{i,1}(subjNo) tmpData_memory{i,2}(subjNo)],'-','LineWidth',2); hold on
            end
            ylim(ylimRange);
            set(gca,'xtick',[1 2])
            set(gca,'xticklabel',{'lure','teleport'});
            ylabel(vLabel)
            title({condiName{i} strcat('p=',num2str(p))});
        end
        [~,p1,~,stats] = ttest2(tmpData_memory{1,1},tmpData_memory{2,1});
        [~,p2,~,stats] = ttest2(tmpData_memory{1,2},tmpData_memory{2,2});
        [~,p3,~,stats] = ttest2(tmpData_memory{1,1}-tmpData_memory{1,2},tmpData_memory{2,1}-tmpData_memory{2,2});
        suptitle({strcat(condiName{1},' vs.',condiName{2}) strcat(' Lure p=',num2str(p1),'teleport p=',num2str(p2)) ...
            strcat('interaction p=',num2str(p3))});
    end
end

% perform mixed ANOVA
addpath('H:\matlabCode\with&withoutpic\privateFunctions');
subjIdx = find(geneInfo<3);
medianData = median(ageInfo);
ageInfo1 = ageInfo(subjIdx);
genderInfo1 = genderInfo(subjIdx);
genderInfo2 = zeros(length(subjIdx),1);
genderInfo2(ismember(genderInfo1,'f')) = 1;
genderInfo2(ismember(genderInfo1,'m')) = 2;
ageInfo1(ageInfo1<medianData) = 1;
ageInfo1(ageInfo1>=medianData) = 2;
betweenSubVar = [];
betweenSubVar = cat(2,betweenSubVar,geneInfo(subjIdx),ageInfo1,genderInfo2);
dataANOVA = NaN(sum(geneInfo<3),2,2);
dataANOVA(:,1,1) = subjData(subjIdx,1); % spatial lure
dataANOVA(:,2,1) = subjData(subjIdx,2); % item lure
dataANOVA(:,1,2) = subjData(subjIdx,3); % spatial Tel
dataANOVA(:,2,2) = subjData(subjIdx,4); % item Tel

anovaResult = simple_mixed_anova(dataANOVA, betweenSubVar, ...
    {'spatial_item','lure_teleport'},{'carrierOrnot','youngOld','gender'});
    


% correlation between performance and age seperately for carrier and non
% carriers
geneCondiName = {'Carrier','noncarrier'};
ageCondiName = {'young','old'};
dataColor = {'r','g'};
testCondiName = {'lure_spatial','lure_item','teleport_spatial','teleport_item'};
figure;
for condiNo = 1:length(testCondiName)
    subplot(2,2,condiNo);
    for geneNo = 1:2
        subjIdx = find(geneInfo==geneNo);
        ageData = ageInfo(subjIdx);
        behVar = subjData(subjIdx,condiNo);
        [r,p] = corr(behVar,ageData,'type','Spearman');
        hold on
        plot(ageData,behVar,strcat(dataColor{geneNo},'.'),'markerSize',10,'HandleVisibility','off'); hold on;
        P = polyfit(ageData,behVar,1); 
        yfit = polyval(P,ageData); 
        hold on
        plot(ageData,yfit,strcat(dataColor{geneNo},'-'));
        text(ageData(1),behVar(1),{strcat('r=',num2str(r)) strcat('p=',num2str(p))},'Color',dataColor{geneNo},'FontSize',14)
    end
    legend(geneCondiName)
    title(testCondiName{condiNo});
end

% the relation between item memory and spatial memory
ageData_median = median(ageInfo(geneInfo<3));
figure;
for geneNo = 1:2
    subjIdx = find(geneInfo==geneNo);
    for ageNo = 1:2
        subplot(2,2,(geneNo-1)*2+ageNo);
        if ageNo == 1
            subjIdx = find(geneInfo==geneNo & ageInfo<=ageData_median);
        else
            subjIdx = find(geneInfo==geneNo & ageInfo>ageData_median);
        end
        for i = 1:2
            var1 = subjData(subjIdx,(i-1)*2+1);
            var2 = subjData(subjIdx,(i-1)*2+2);
            [r,p] = corr(var1,var2,'type','Spearman');
            hold on
            plot(var1,var2,strcat(dataColor{i},'.'),'markerSize',10,'HandleVisibility','off'); hold on;
            P = polyfit(var1,var2,1); 
            yfit = polyval(P,var1); 
            hold on
            plot(var1,yfit,strcat(dataColor{i},'-'));
            xlabel('drop error');
            ylabel('item memory');
            text(var1(1),var2(1),{strcat('r=',num2str(r)) strcat('p=',num2str(p))},'Color',dataColor{i},'FontSize',14)
        end
        legend({'lure','teleport'});
        title(strcat(geneCondiName{geneNo},'_ ',ageCondiName{ageNo}));
    end
end


figure;
for geneNo = 1:2
    subjIdx = find(geneInfo==geneNo);
    subplot(1,2,geneNo);
    for i = 1:2
        var1 = subjData(subjIdx,(i-1)*2+1);
        var2 = subjData(subjIdx,(i-1)*2+2);
        [r,p] = corr(var1,var2,'type','Spearman');
        hold on
        plot(var1,var2,strcat(dataColor{i},'.'),'markerSize',10,'HandleVisibility','off'); hold on;
        P = polyfit(var1,var2,1); 
        yfit = polyval(P,var1); 
        hold on
        plot(var1,yfit,strcat(dataColor{i},'-'));
        xlabel('drop error');
        ylabel('item memory');
        text(var1(1),var2(1),{strcat('r=',num2str(r)) strcat('p=',num2str(p))},'Color',dataColor{i},'FontSize',14)
    end
    legend({'lure','teleport'});
    title(geneCondiName{geneNo});
end


%% to calculate teh relation between the incoming distance and the navigation time
subjNoAll = size(behInfo,2);
dataAll = struct;
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    trialdata = zeros(trialNoAll,8);
    for trialNo = 1:trialNoAll
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;
        incomingXY = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
        incomingTime = behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ... 
            behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime;

        pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
        dropError = pdist([targetLoc;dropLoc],'euclidean');
        travelLenth = pdist([startSearchLoc;dropLoc],'euclidean');
        incomingPath = sum(arrayfun(@(x) pdist([incomingXY(x,:);incomingXY(x+1,:)],'euclidean'),1:size(incomingXY,1)-1));
        trialdata(trialNo,1) = travelLenth;
        trialdata(trialNo,2) = dropError;%dropError;
        trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
        trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
        trialdata(trialNo,5) = pathLenth;
        trialdata(trialNo,6) = incomingTime;
        trialdata(trialNo,7) = incomingPath;
    end
    dataAll(subjNo).trialdata = trialdata;
end

conditionLabel = {'lure','teleport'};
conditionNoLabel = [1 0];


for nameIdx = 1:2
    if nameIdx == 1
        ageInfo1 = ageInfo(~isnan(ageInfo));
        medianData = median(ageInfo1);
        subjIdx1 = find(ageInfo<medianData & geneInfo<3); 
        subjIdx2 = find(ageInfo>=medianData & geneInfo<3);
        condiName = {'young','old'};
    elseif nameIdx == 2
        subjIdx1 = find(geneInfo==1); % carrier
        subjIdx2 = find(geneInfo==2); % non-carrier
        condiName = {'Carrier','noncarrier'};
    end
    subjIdx3 = 1:subjNoAll;
    subjIdxPool = {'subjIdx1','subjIdx2','subjIdx3'};

    figure;
    for subjIdxNo = 1:2
        subjIdx = eval(subjIdxPool{subjIdxNo});
        for conditionNo = 1:length(conditionLabel)
            tmp_label = conditionLabel{conditionNo};
            tmpConNo = conditionNoLabel(conditionNo);
            
            tmpr = NaN(1,length(subjIdx)); 
            for subjNo = 1:length(subjIdx)
                tmpidx1 = find(dataAll(subjIdx(subjNo)).trialdata(:,4)==tmpConNo);
                dropError = dataAll(subjIdx(subjNo)).trialdata(tmpidx1,2);
                incomingDist = dataAll(subjIdx(subjNo)).trialdata(tmpidx1,5);
                incomingTime = dataAll(subjIdx(subjNo)).trialdata(tmpidx1,6);
                tmpr(subjNo) = partialcorr(dropError,incomingDist,incomingTime);
            end
            subplot(2,2,(subjIdxNo-1)*2+conditionNo);
            bar(1,nanmean(tmpr)); hold on
            plot(ones(1,length(tmpr)),tmpr,'.','markerSize',30);
            ylim([-0.2 0.8])
            [~,p,~,stat] = ttest(tmpr);
            title(strcat(condiName{subjIdxNo},tmp_label,'-partial correlation-',num2str(p)));
        end
    end
    suptitle('within sub drop error corr incoming Dis partial out time')
    
    
    
    for subjIdxNo = 1:2
        subjIdx = eval(subjIdxPool{subjIdxNo});
        dropErrorDis_itemMem = NaN(length(subjIdx),4,2);
        for conditionNo = 1:length(conditionLabel)
            tmp_label = conditionLabel{conditionNo};
            tmpConNo = conditionNoLabel(conditionNo);
            for correctNo = 1:2
                for subjNo = 1:length(subjIdx)
                    tmpidx1 = find(dataAll(subjIdx(subjNo)).trialdata(:,4)==tmpConNo & dataAll(subjIdx(subjNo)).trialdata(:,3)==(2-correctNo));
                    dropErrorDis_itemMem(subjNo,1,conditionNo,correctNo) = nanmean(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,6));
                    dropErrorDis_itemMem(subjNo,2,conditionNo,correctNo) = nanmean(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,5));
                    dropErrorDis_itemMem(subjNo,3,conditionNo,correctNo) = nanmean(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,7));
                    dropErrorDis_itemMem(subjNo,4,conditionNo,correctNo) = nanmean(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,2));
                end
            end
        end
        % incoming time
        indexPool = {'incoming time','incoming dist_ideal','incoming dist_prac','drop error'};
        figure;
        for condiNo = 1:length(conditionLabel)
            tmp_label = conditionLabel{condiNo};
            for idxNo = 1:length(indexPool)
                tmpData = squeeze(dropErrorDis_itemMem(:,idxNo,condiNo,:));
                [~,p,~,stat] = ttest(tmpData(:,1),tmpData(:,2));
                meanData = nanmean(tmpData,1);
                stdData = nanstd(tmpData)/sqrt(subjNoAll);
                subplot(2,4,(condiNo-1)*4+idxNo);
                bar(1:2,meanData); hold on
                errorbar(1:2,meanData,stdData,'k.');
                set(gca,'xtick',[1 2]);
                set(gca,'xticklabel',{'Remember','Forgotten'});
                xlabel('item memory');
                ylabel(indexPool{idxNo});
                title(strcat(tmp_label,num2str(p)));
            end
        end
        suptitle(strcat(condiName{subjIdxNo}))
    end
end

figure;
for nameIdx = 1:2
    if nameIdx == 1
        ageInfo1 = ageInfo(~isnan(ageInfo));
        medianData = median(ageInfo1);
        subjIdx1 = find(ageInfo<medianData & geneInfo<3); 
        subjIdx2 = find(ageInfo>=medianData & geneInfo<3);
        condiName = {'young','old'};
    elseif nameIdx == 2
        subjIdx1 = find(geneInfo==1); % carrier
        subjIdx2 = find(geneInfo==2); % non-carrier
        condiName = {'Carrier','noncarrier'};
    end
    subjIdx3 = 1:subjNoAll;
    subjIdxPool = {'subjIdx1','subjIdx2','subjIdx3'};

    for subjIdxNo = 1:2
        subjIdx = eval(subjIdxPool{subjIdxNo});
        dropErrorDis_itemMem = NaN(length(subjIdx),4,2);
        tmpr = NaN(length(subjIdxNo),2,2);
        for conditionNo = 1:length(conditionLabel)
            tmp_label = conditionLabel{conditionNo};
            tmpConNo = conditionNoLabel(conditionNo);
            for subjNo = 1:length(subjIdx)
                tmpidx1 = find(dataAll(subjIdx(subjNo)).trialdata(:,4)==tmpConNo);
                tmpr(subjNo,conditionNo,1) = corr(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,2),dataAll(subjIdx(subjNo)).trialdata(tmpidx1,1),'type','Pearson');
                tmpr(subjNo,conditionNo,2) = partialcorr(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,2),...
                    dataAll(subjIdx(subjNo)).trialdata(tmpidx1,1),...
                    dataAll(subjIdx(subjNo)).trialdata(tmpidx1,6),'type','Pearson');
            end
        end

        for i = 1:2
            subplot(2,4,(nameIdx-1)*2+subjIdxNo+(i-1)*4);
            [~,p,~,stat] = ttest(tmpr(:,1,i),tmpr(:,2,i));
            meanData = squeeze(nanmean(tmpr(:,:,i)));
            stdData=squeeze(nanstd(tmpr(:,:,i)))/sqrt(length(subjIdx));
            bar(1:2,meanData);
            hold on
            errorbar(1:2,meanData,stdData,'k.');
            set(gca,'xtick',[1 2]);
            set(gca,'xticklabel',conditionLabel);
            if i==1
                ylabel('travel path corr drop error')
            else
                ylabel('travel path corr drop error out time')
            end
            title(strcat(condiName{subjIdxNo},'p=',num2str(p)))
        end
    end
end
    
%% distance to the nearest landmark
mountainLoc = [-16 16]; %[-20 20];
largeShip = [-16 -16];%[-20 -20];
smallShip = [16 -16];%[20 -20];

dataAll = struct;
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);

    trialdata = zeros(trialNoAll,8);
    for trialNo = 1:trialNoAll
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;
        incomingXY = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
        incomingTime = behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ... 
            behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime;

        dist_landmark = min([pdist([targetLoc;mountainLoc],'euclidean'),...
            pdist([targetLoc;largeShip],'euclidean') ...
            pdist([targetLoc;smallShip],'euclidean')]);
        pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
        dropError = pdist([targetLoc;dropLoc],'euclidean');
        travelLenth = pdist([startSearchLoc;dropLoc],'euclidean');
        incomingPath = sum(arrayfun(@(x) pdist([incomingXY(x,:);incomingXY(x+1,:)],'euclidean'),1:size(incomingXY,1)-1));
        trialdata(trialNo,1) = travelLenth;
        trialdata(trialNo,2) = dropError;%dropError;
        trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
        trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
        trialdata(trialNo,5) = dist_landmark;
        trialdata(trialNo,6) = incomingTime;
        trialdata(trialNo,7) = incomingPath;
    end
    dataAll(subjNo).trialdata = trialdata;
end



figure;
for nameIdx = 1:2
    if nameIdx == 1
        ageInfo1 = ageInfo(~isnan(ageInfo));
        medianData = median(ageInfo1);
        subjIdx1 = find(ageInfo<medianData & geneInfo<3); 
        subjIdx2 = find(ageInfo>=medianData & geneInfo<3);
        condiName = {'young','old'};
    elseif nameIdx == 2
        subjIdx1 = find(geneInfo==1); % carrier
        subjIdx2 = find(geneInfo==2); % non-carrier
        condiName = {'Carrier','noncarrier'};
    end
    subjIdx3 = 1:subjNoAll;
    subjIdxPool = {'subjIdx1','subjIdx2','subjIdx3'};

    tmpr = NaN(sum(geneInfo<3),2,2);
    for subjIdxNo = 1:2
        subjIdx = eval(subjIdxPool{subjIdxNo});
        dropErrorDis_itemMem = NaN(length(subjIdx),4,2);
        for conditionNo = 1:length(conditionLabel)
            tmp_label = conditionLabel{conditionNo};
            tmpConNo = conditionNoLabel(conditionNo);
            for subjNo = 1:length(subjIdx)
                tmpidx1 = find(dataAll(subjIdx(subjNo)).trialdata(:,4)==tmpConNo);
                tmpr(subjNo,conditionNo,subjIdxNo) = corr(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,2),...
                    dataAll(subjIdx(subjNo)).trialdata(tmpidx1,5),'type','Pearson');
            end
        end
    end
        
    subplot(1,2,nameIdx);
    meanD = squeeze(nanmean(tmpr));
    stdD=squeeze(nanstd(tmpr))/sqrt(length(subjIdx));
    ngroups = size(meanD, 1);
    nbars = size(meanD, 2);
    tmpcolor = colormap('jet');
    ncolors = tmpcolor(round(linspace(1,64,nbars)),:);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        bar(x,meanD(:,i),groupwidth/2,'FaceColor',ncolors(i,:),'EdgeColor',ncolors(i,:),'LineWidth',1.5); hold on
        errorbar(x, meanD(:,i), stdD(:,i), '.','LineWidth',3, 'color','k','HandleVisibility','off'); hold on
    end
    hold off
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',conditionLabel);
    legend(condiName)
    ylabel('drop error corr. dis to landmark')
    set(gca,'fontsize',14)
end



addpath('H:\matlabCode\with&withoutpic\privateFunctions');
subjIdx = find(geneInfo<3);
medianData = median(ageInfo);
ageInfo1 = ageInfo(subjIdx);
ageInfo1(ageInfo1<medianData) = 1;
ageInfo1(ageInfo1>=medianData) = 2;
betweenSubVar = cat(2,geneInfo(subjIdx),ageInfo1);
tmpr = NaN(sum(geneInfo<3),2,2);
dataANOVA = NaN(sum(geneInfo<3),2);
for subjNo = 1:length(subjIdx)
    for conditionNo = 1:length(conditionLabel)
        tmp_label = conditionLabel{conditionNo};
        tmpConNo = conditionNoLabel(conditionNo);
        tmpidx1 = find(dataAll(subjIdx(subjNo)).trialdata(:,4)==tmpConNo);
        dataANOVA(subjNo,conditionNo) = corr(dataAll(subjIdx(subjNo)).trialdata(tmpidx1,2),...
                    dataAll(subjIdx(subjNo)).trialdata(tmpidx1,5),'type','Pearson');
    end
end

anovaResult = simple_mixed_anova(dataANOVA, betweenSubVar, {'lure_tel'},{'carrierOrnot','youngOld'});
    
figure;
for nameIdx = 1:2
    if nameIdx == 1
        condiName = {'young','old'};
        condiInfo = ageInfo1;
    elseif nameIdx == 2
        condiName = {'Carrier','noncarrier'};
        condiInfo = geneInfo(geneInfo<3);
    end

    tmpData1 = dataANOVA(condiInfo==1,:);
    tmpData2 = dataANOVA(condiInfo==2,:);
    [~,p1,~,~] = ttest(tmpData1(:,1),tmpData1(:,2));
    [~,p2,~,~] = ttest(tmpData2(:,1),tmpData2(:,2));
    [~,p3,~,~] = ttest2(tmpData1(:,1),tmpData2(:,1));
    [~,p4,~,~] = ttest2(tmpData1(:,2),tmpData2(:,2));
    [~,p5,~,~] = ttest2(tmpData1(:,1)-tmpData1(:,2),...
        tmpData2(:,1)-tmpData2(:,2));

    meanD = zeros(2,2);
    meanD(:,1) = nanmean(tmpData1);
    meanD(:,2) = nanmean(tmpData2);
    stdD = zeros(2,2);
    stdD(:,1) = nanstd(tmpData1)/sqrt(size(tmpData1,1));
    stdD(:,2) = nanstd(tmpData2)/sqrt(size(tmpData2,1));
    
    subplot(1,2,nameIdx);
    ngroups = size(meanD, 1);
    nbars = size(meanD, 2);
    tmpcolor = colormap('jet');
    ncolors = tmpcolor(round(linspace(1,64,nbars)),:);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        bar(x,meanD(:,i),groupwidth/2,'FaceColor',ncolors(i,:),'EdgeColor',ncolors(i,:),'LineWidth',1.5); hold on
        errorbar(x, meanD(:,i), stdD(:,i), '.','LineWidth',3, 'color','k','HandleVisibility','off'); hold on
    end
    plot([0.8 1.8],[mean(meanD(:,1)) mean(meanD(:,1))],'k-'); hold on
    text(1.3,mean(meanD(:,1)),num2str(p1));
    plot([1.2 2.2],[mean(meanD(:,2)) mean(meanD(:,2))],'k-'); hold on
    text(1.7,mean(meanD(:,2)),num2str(p2));
    plot([0.8 1.2],[mean(meanD(1,:)) mean(meanD(1,:))],'k-'); hold on
    text(1,mean(meanD(1,:)),num2str(p3));
    plot([1.8 2.2],[mean(meanD(2,:)) mean(meanD(2,:))],'k-'); hold on
    text(2,mean(meanD(2,:)),num2str(p3));
    
    for i = 1:size(tmpData1,1)
    end
    
    hold off
    set(gca,'xtick',[1 2]);
    set(gca,'xticklabel',conditionLabel);
    legend(condiName);
    ylabel('drop error corr. dis to landmark');
    set(gca,'fontsize',14);
    title(num2str(p5))
end


%% calculate the navigatning pattern
mountainLoc = [-16 16]; %[-20 20];
largeShip = [-16 -16];%[-20 -20];
smallShip = [16 -16];%[20 -20];

dataAll = zeros(subjNoAll,2,6);
for subjNo = 1:subjNoAll
    disp(subjNo)
    trialNoAll = length(behInfo(subjNo).trialInfo);
    trialdata = zeros(trialNoAll,13);
    for trialNo = 1:trialNoAll
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;
        incomingXY = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
        incomingXYTime = behInfo(subjNo).trialInfo(trialNo).posXYTime_startToDeliver_dropTarget;
        incomingSpeed1 = arrayfun(@(x) pdist([incomingXY(x,:);incomingXY(x+1,:)],'euclidean'),1:size(incomingXY,1)-1);
        incomingSpeed2 = arrayfun(@(x) incomingXYTime(x+1)-incomingXYTime(x),1:size(incomingXY,1)-1);
        incomingSpeed3 = arrayfun(@(x) incomingSpeed1(x+1)+incomingSpeed1(x),1:size(incomingSpeed1,2)-1);
        incomingSpeed4 = arrayfun(@(x) incomingSpeed2(x+1)+incomingSpeed2(x),1:size(incomingSpeed2,2)-1);
        incomingSpeed5 = [0 incomingSpeed3./incomingSpeed4 0];
        
        trialdata(trialNo,8) = sum(abs(incomingXY(:,1))>10 | abs(incomingXY(:,2))>10);
        trialdata(trialNo,9) = sum((abs(incomingXY(:,1))<=10 & abs(incomingXY(:,1))>5) ...
            | (abs(incomingXY(:,2))<=10 & abs(incomingXY(:,2))>5));
        trialdata(trialNo,10) = sum(abs(incomingXY(:,1))<=5 | abs(incomingXY(:,2))<=5);
        
        incomingXY1 = incomingXY(incomingSpeed5>0,:);
        trialdata(trialNo,11) = sum(abs(incomingXY1(:,1))>10 | abs(incomingXY1(:,2))>10);
        trialdata(trialNo,12) = sum((abs(incomingXY1(:,1))<=10 & abs(incomingXY1(:,1))>5) ...
            | (abs(incomingXY1(:,2))<=10 & abs(incomingXY1(:,2))>5));
        trialdata(trialNo,13) = sum(abs(incomingXY1(:,1))<=5 | abs(incomingXY1(:,2))<=5);
        
        incomingTime = behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ... 
            behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime;

        dist_landmark = min([pdist([targetLoc;mountainLoc],'euclidean'),...
            pdist([targetLoc;largeShip],'euclidean') ...
            pdist([targetLoc;smallShip],'euclidean')]);
        pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
        dropError = pdist([targetLoc;dropLoc],'euclidean');
        travelLenth = pdist([startSearchLoc;dropLoc],'euclidean');
        incomingPath = sum(arrayfun(@(x) pdist([incomingXY(x,:);incomingXY(x+1,:)],'euclidean'),1:size(incomingXY,1)-1));
        trialdata(trialNo,1) = travelLenth;
        trialdata(trialNo,2) = dropError;%dropError;
        trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
        trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
        trialdata(trialNo,5) = dist_landmark;
        trialdata(trialNo,6) = incomingTime;
        trialdata(trialNo,7) = incomingPath;
    end
    for condiNo = 1:2
        trialIdx = find(trialdata(:,4)==(2-condiNo));
        dataAll(subjNo,condiNo,1:6) = sum(trialdata(trialIdx,8:13));
    end
end

dataAll_v1(:,:,1:3) = dataAll(:,:,1:3)./repmat(nansum(dataAll(:,:,1:3),3),[1 1 3]);
dataAll_v1(:,:,4:6) = dataAll(:,:,4:6)./repmat(nansum(dataAll(:,:,4:6),3),[1 1 3]);

tmp_index = 4:6;

addpath('H:\matlabCode\with&withoutpic\privateFunctions');
subjIdx = find(geneInfo<3);
medianData = median(ageInfo);
ageInfo1 = ageInfo(subjIdx);
ageInfo1(ageInfo1<medianData) = 1;
ageInfo1(ageInfo1>=medianData) = 2;
betweenSubVar = cat(2,geneInfo(subjIdx),ageInfo1);
tmpr = NaN(sum(geneInfo<3),2,2);
dataANOVA = dataAll_v1(subjIdx,:,tmp_index);
anovaResult = simple_mixed_anova(dataANOVA, betweenSubVar, {'lure_tel','edge'},{'carrierOrnot','youngOld'});
    
lureTelName = {'lure','teleport'};
figure;
for nameIdx = 1:2
    if nameIdx == 1
        condiName = {'young','old'};
        condiInfo = ageInfo1;
    elseif nameIdx == 2
        condiName = {'Carrier','noncarrier'};
        condiInfo = geneInfo(geneInfo<3);
    end

    
    for LureTel = 1:2
        tmpData1 = dataAll_v1(condiInfo==1,LureTel,tmp_index);
        tmpData2 = dataAll_v1(condiInfo==2,LureTel,tmp_index);

        [~,p1,~,~] = ttest2(tmpData1(:,1),tmpData2(:,1));
        [~,p2,~,~] = ttest2(tmpData1(:,3),tmpData2(:,3));
        [~,p3,~,~] = ttest2(tmpData1(:,1)-tmpData1(:,3),...
            tmpData2(:,1)-tmpData2(:,3));
        
        meanD = zeros(2,3);
        meanD(1,:) = nanmean(tmpData1);
        meanD(2,:) = nanmean(tmpData2);
        stdD = zeros(2,3);
        stdD(1,:) = nanstd(tmpData1)/sqrt(size(tmpData1,1));
        stdD(2,:) = nanstd(tmpData2)/sqrt(size(tmpData2,1));
    
        subplot(2,2,(nameIdx-1)*2+LureTel);
        ngroups = size(meanD, 1);
        nbars = size(meanD, 2);
        tmpcolor = colormap('jet');
        ncolors = tmpcolor(round(linspace(1,64,nbars)),:);
        % Calculating the width for each bar group
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            bar(x,meanD(:,i),groupwidth/2,'FaceColor',ncolors(i,:),'EdgeColor',ncolors(i,:),'LineWidth',1.5); hold on
            errorbar(x, meanD(:,i), stdD(:,i), '.','LineWidth',3, 'color','k','HandleVisibility','off'); hold on
        end

        hold off
        set(gca,'xtick',[1 2]);
        set(gca,'xticklabel',condiName);
        legend({'ege','mid','ctr'});
        ylabel('data points');
        set(gca,'fontsize',14);
        ylim([0 0.6])
        title({strcat(lureTelName{LureTel}) strcat('edgeP=',num2str(p1),' ctrP=',num2str(p2)) ...
            strcat('interactionP=',num2str(p3))})
    end
end








































































































%% distance to the nearest landmark
mountainLoc = [-16 16]; %[-20 20];
largeShip = [-16 -16];%[-20 -20];
smallShip = [16 -16];%[20 -20];

dataAll = struct;
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    if trialNoAll >70
        trialNoAll = 70;
        trialdata = zeros(trialNoAll,8);
        for trialNo = 1:trialNoAll
            startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
            targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
            dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;
            incomingXY = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
            incomingTime = behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ... 
                behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime;

            dist_landmark = min([pdist([targetLoc;mountainLoc],'euclidean'),...
                pdist([targetLoc;largeShip],'euclidean') ...
                pdist([targetLoc;smallShip],'euclidean')]);
            pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
            dropError = pdist([targetLoc;dropLoc],'euclidean');
            travelLenth = pdist([startSearchLoc;dropLoc],'euclidean');
            incomingPath = sum(arrayfun(@(x) pdist([incomingXY(x,:);incomingXY(x+1,:)],'euclidean'),1:size(incomingXY,1)-1));
            trialdata(trialNo,1) = travelLenth;
            trialdata(trialNo,2) = dropError;%dropError;
            trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
            trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
            trialdata(trialNo,5) = dist_landmark;
            trialdata(trialNo,6) = incomingTime;
            trialdata(trialNo,7) = incomingPath;
        end
        dataAll(subjNo).trialdata = trialdata;
    else
        dataAll(subjNo).trialdata = [];
    end
end


%% GCLR calculation
addpath('F:\project\matlabCode\hui_code');
dataAll = struct;
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    trialdata = zeros(trialNoAll,4);
    for trialNo = 1:trialNoAll
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;
        incomingXY = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
        incomingTime = behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ... 
            behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime;
        
        incomingAngle = atan2(targetLoc(2)-startSearchLoc(2),targetLoc(1)-startSearchLoc(1));
        dropError = pdist([targetLoc;dropLoc],'euclidean');
        trialdata(trialNo,1) = incomingAngle;
        trialdata(trialNo,2) = dropError;%dropError;
        trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
        trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
    end
    dataAll(subjNo).trialdata = trialdata;
end

foldNo = 6;
figure;
for subjNo = 1:subjNoAll
    tmpAngle = dataAll(subjNo).trialdata(:,1);
    dropError = dataAll(subjNo).trialdata(:,2);
    
    for rptNo = 1:rpt
        powerAll1 = dropError(randperm(length(dropError)));
        [betaValues, preferDir] = GCLRfunction( powerAll1, tmpAngle, 'rad', foldNo);
        tmpData(1:3,rptNo) = betaValues;
        tmpData(4,rptNo) = preferDir;
    end
    [betaValues, preferDir] = GCLRfunction( dropError, tmpAngle, 'rad', foldNo);
    subplot(3,4,subjNo);
    plot(1:rpt,squeeze(tmpData(2:3,:))); hold on
%     plot(1:rpt,betaValues(1)*ones(1,rpt),'r--'); hold on
    plot(1:rpt,betaValues(2)*ones(1,rpt),'g--'); hold on
    plot(1:rpt,betaValues(3)*ones(1,rpt),'b--'); hold on
%     plot(1:rpt,preferDir*ones(1,rpt),'k--'); hold on
    legend
end





foldNo = 4;
runningDir_deg = 0:5:360;
runningDir_rad = runningDir_deg*pi/180;
runningDir_rad_60 = runningDir_rad *foldNo;
preDirs = (5:(360/foldNo):360)*pi/180;
powerAll = zeros(size(runningDir_rad));
for idx = 1:length(powerAll)
    tmpAng = min(abs(runningDir_rad(idx) - preDirs));
    if tmpAng<30*pi/180
        tmpAng = 360/foldNo*pi/180-tmpAng;
    end
    powerAll(idx) =  5*tmpAng/(360/foldNo*pi/180);
end

figure;
bar(runningDir_deg,powerAll)

rpt = 100;
amp_tmp = 0:0.5:10;
tmpData = zeros(4,length(amp_tmp),rpt);
foldNo = 6;
for rptNo = 1:rpt
    for ampNo = 1:length(amp_tmp)
        powerAll1 = powerAll + rand(size(powerAll))*amp_tmp(ampNo);
        [betaValues, preferDir] = GCLRfunction( powerAll1, runningDir_rad, 'rad', foldNo);
        tmpData(1:3,ampNo,rptNo) = betaValues;
        tmpData(4,ampNo,rptNo) = preferDir;
    end
end

figure;
legendLabel = cell(1,length(amp_tmp));
for ampNo = 1:length(amp_tmp)
    plot(1:rpt,squeeze(tmpData(3,ampNo,:))); hold on
    legendLabel{ampNo} = num2str(amp_tmp(ampNo));
end
legend(legendLabel)


tmpData = zeros(4,rpt);
foldNo = 6;
for rptNo = 1:rpt
    powerAll1 = powerAll(randperm(length(powerAll)));
    [betaValues, preferDir] = GCLRfunction( powerAll1, runningDir_rad, 'rad', foldNo);
    tmpData(1:3,rptNo) = betaValues;
    tmpData(4,rptNo) = preferDir;
end
[betaValues, preferDir] = GCLRfunction( powerAll, runningDir_rad, 'rad', foldNo);
figure;
plot(1:rpt,squeeze(tmpData(:,:))); hold on
plot(1:rpt,betaValues(1)*ones(1,rpt),'r--'); hold on
plot(1:rpt,betaValues(2)*ones(1,rpt),'g--'); hold on
plot(1:rpt,betaValues(3)*ones(1,rpt),'b--'); hold on
plot(1:rpt,preferDir*ones(1,rpt),'k--'); hold on
legend
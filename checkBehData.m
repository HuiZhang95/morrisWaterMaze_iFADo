dataPath = 'F:\project\navigationParadigm\MorrisMaze\data\ifADo Data\';
saveDataPath = 'F:\project\navigationParadigm\MorrisMaze\analyzeBehData\';
dinfo = dir(dataPath);
dinfo(ismember( {dinfo.name}, {'.', '..'})) = [];
subjNoAll = length(dinfo);
behInfo = struct;
dataPath_exl = 'F:\project\navigationParadigm\MorrisMaze\data\';
[num,txt,raw] = xlsread(strcat(dataPath_exl,'Demographics_04092020_IFaDo.xlsx'));
for subjNo = 1:subjNoAll
    disp(subjNo);
    tic
    if contains(dinfo(subjNo).name,'_m_')
        behInfo(subjNo).sex = 'm';
    elseif contains(dinfo(subjNo).name,'_f_')
        behInfo(subjNo).sex = 'f';
    end
    subjIdx = find(ismember(raw(:,1),dinfo(subjNo).name(7:end)));
    behInfo(subjNo).age = raw{subjIdx,2};
    behInfo(subjNo).geneType = raw{subjIdx,4};
    behInfo(subjNo).geneType2 = raw{subjIdx,5};
    
    behInfo(subjNo).subjID = dinfo(subjNo).name;
    subPath = strcat(dataPath,dinfo(subjNo).name,filesep,'session_0',filesep);
    fid = fopen(strcat(subPath,'log.txt'));
    fid_text = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',{'\t'});
    % trigger information
    eventLineNo = find(ismember(fid_text{3},'bigEvent'));
    eventName = fid_text{4}(eventLineNo);
    eventTime = str2double(fid_text{1}(eventLineNo));
    
    startRealExperimentLineNo = find(ismember(fid_text{3},'bigEvent') & ismember(fid_text{4},'startRealExperimentTrigger'));
    trialStartLineNo = find(ismember(fid_text{3},'bigEvent') & ismember(fid_text{4},'trialStartTrigger'));
    trialStartLineNo(end) = [];
    trialEndLineNo = find(ismember(fid_text{3},'bigEvent') & contains(fid_text{4},'trialInfo:'));
    
    if ~isempty(startRealExperimentLineNo)
        trialStartLineNo(trialStartLineNo<startRealExperimentLineNo) = [];
        trialEndLineNo(trialEndLineNo<startRealExperimentLineNo) = [];
    end
        
    
    % avatar position information
    posXYLineNo = find(ismember(fid_text{4},'PandaEPL_avatar') & ismember(fid_text{3},'VROBJECT_POS'));
    posXYTime = str2double(fid_text{1}(posXYLineNo));
    [tmpidx1,tmpidx2] = unique(posXYTime);
    posXYTime = tmpidx1;
    posXYLineNo = posXYLineNo(tmpidx2);
    posXY_heading_LineNo = find(ismember(fid_text{4},'PandaEPL_avatar') & ismember(fid_text{3},'VROBJECT_HEADING'));
    posXYTime_heading = str2double(fid_text{1}(posXY_heading_LineNo));
    posXY_turningSpeedLineNo = find(ismember(fid_text{4},'PandaEPL_avatar') & ismember(fid_text{3},'MOVINGOBJECT_TURNINGSPEED'));
    posXYTime_turningSpeed = str2double(fid_text{1}(posXY_turningSpeedLineNo));
    posXY_linearSpeedLineNo = find(ismember(fid_text{4},'PandaEPL_avatar') & ismember(fid_text{3},'MOVINGOBJECT_LINEARSPEED'));
    posXYTime_linearSpeed = str2double(fid_text{1}(posXY_linearSpeedLineNo));
    
    posXY = zeros(length(posXYTime),2);
    posXY_heading = NaN(length(posXYTime),1);
    posXY_tuningSpeed = NaN(length(posXYTime),1);
    posXY_linearSpeed = NaN(length(posXYTime),1);
    for idx = 1:length(posXYTime)
        tmptext = fid_text{5}{posXYLineNo(idx)};
        tmpText1 = strsplit(tmptext,{',','(',')'});
        posXY(idx,:) = [str2double(tmpText1{2}) str2double(tmpText1{3})];
        tmpidx = find(posXYTime_heading==posXYTime(i));
        if ~isempty(tmpidx)
            tmptext = fid_text{5}{posXY_heading_LineNo(tmpidx(1))};
            posXY_heading(idx) = str2double(tmptext);
        end
        tmpidx = find(posXYTime_turningSpeed==posXYTime(i));
        if ~isempty(tmpidx)
            tmptext = fid_text{5}{posXY_turningSpeedLineNo(tmpidx(1))};
            posXY_tuningSpeed(idx) = str2double(tmptext);
        end
        tmpidx = find(posXYTime_linearSpeed==posXYTime(i));
        if ~isempty(tmpidx)
            tmptext = fid_text{5}{posXY_linearSpeedLineNo(tmpidx(1))};
            posXY_linearSpeed(idx) = str2double(tmptext);
        end
    end
%     posXY = cat(2,str2double(fid_text{6}(posXYLineNo)),str2double(fid_text{7}(posXYLineNo)));
    
    trialNoAll = length(trialStartLineNo);
    trialInfo = struct;
    for trialNo = 1:trialNoAll
        tmpidx = find(eventLineNo>trialStartLineNo(trialNo) & eventLineNo<trialEndLineNo(trialNo));
        eventNames = eventName(tmpidx);

        startToTarget = trialStartLineNo(trialNo);
        startSearchTime = str2double(fid_text{1}(startToTarget));
        trialInfo(trialNo).startTime = startSearchTime;
        
        collectTarget = eventLineNo(tmpidx(ismember(eventNames,'collectTargetTrigger')));
        collectTargetTime = str2double(fid_text{1}(collectTarget));
        trialInfo(trialNo).collectTargetTime = startSearchTime;
        
        tmpXY_idx = find(posXYLineNo>startToTarget & posXYLineNo<collectTarget);
        trialInfo(trialNo).posXYTime_start_collectTarget = posXYTime(tmpXY_idx);
        trialInfo(trialNo).posXY_start_collectTarget = posXY(tmpXY_idx,:);

        trialInfo(trialNo).posXY_heading_start_collectTarget = posXY_heading(tmpXY_idx);
        trialInfo(trialNo).posXY_tuningSpeed_start_collectTarget = posXY_tuningSpeed(tmpXY_idx);
        trialInfo(trialNo).posXY_linearSpeed_start_collectTarget = posXY_linearSpeed(tmpXY_idx);
        
        instrToFindTarget = eventLineNo(tmpidx(ismember(eventNames,'instrToFindTargetTrigger')));
        instrToFindTargetTime = str2double(fid_text{1}(instrToFindTarget));
        trialInfo(trialNo).instrToFindTargetTime = instrToFindTargetTime;
        
        if find(ismember(eventNames,'preinstrToFindTargetTrigger'))
            trialInfo(trialNo).trialType = 'lure';% with lure model; 0 is teleported
            tmpXY_idx = find(posXYLineNo>collectTarget & posXYLineNo<instrToFindTarget);
            trialInfo(trialNo).posXYTime_collectTarget_collectLure = posXYTime(tmpXY_idx);
            trialInfo(trialNo).posXY_collectTarget_collectLure = posXY(tmpXY_idx,:);
            
            trialInfo(trialNo).posXY_heading_collectTarget_collectLure = posXY_heading(tmpXY_idx);
            trialInfo(trialNo).posXY_tuningSpeed_collectTarget_collectLure = posXY_tuningSpeed(tmpXY_idx);
            trialInfo(trialNo).posXY_linearSpeed_collectTarget_collectLure = posXY_linearSpeed(tmpXY_idx);

        else
            trialInfo(trialNo).trialType = 'teleport';
        end
        
        startToDeliverTarget = eventLineNo(tmpidx(ismember(eventNames,'findTargetTrigger')));
        startToDeliverTargetTime = str2double(fid_text{1}(startToDeliverTarget));
        trialInfo(trialNo).startToDeliverTargetTime = startToDeliverTargetTime;
        
        dropTarget = eventLineNo(tmpidx(ismember(eventNames,'dropTargetTrigger')));
        dropTargetTime = str2double(fid_text{1}(dropTarget));
        trialInfo(trialNo).dropTargetTime = dropTargetTime;
        
        tmpXY_idx = find(posXYLineNo>startToDeliverTarget & posXYLineNo<dropTarget);
        trialInfo(trialNo).posXYTime_startToDeliver_dropTarget = posXYTime(tmpXY_idx);
        trialInfo(trialNo).posXY_startToDeliver_dropTarget = posXY(tmpXY_idx,:);
        trialInfo(trialNo).posXY_heading_startToDeliver_dropTarget = posXY_heading(tmpXY_idx);
        trialInfo(trialNo).posXY_tuningSpeed_startToDeliver_dropTarget = posXY_tuningSpeed(tmpXY_idx);
        trialInfo(trialNo).posXY_linearSpeed_startToDeliver_dropTarget = posXY_linearSpeed(tmpXY_idx);
        
        tmpText = fid_text{4}{trialEndLineNo(trialNo)};
        tmpText1 = strsplit(tmpText,',');
        tmpText2 = strsplit(tmpText1{2},':');
        trialInfo(trialNo).targetPic = tmpText2{2};
        
        tmpText2 = strsplit(tmpText1{3},'[');
        tmpx = str2double(tmpText2{2});
        tmpText2 = strsplit(tmpText1{4},']');
        tmpy = str2double(tmpText2{1});
        trialInfo(trialNo).targetLoc = [tmpx tmpy];
        trialInfo(trialNo).targetLoc_adjust = trialInfo(trialNo).posXY_start_collectTarget(end,:);
        
        tmpText2 = strsplit(tmpText1{5},':');
        targetPicPos = tmpText2{2};
        tmpText2 = strsplit(tmpText1{6},':');
        responsePicPos = tmpText2{2};
        if strcmp(targetPicPos, responsePicPos)
            trialInfo(trialNo).picMemory = 1;
        else
            trialInfo(trialNo).picMemory = 0;
        end
        
        tmpText2 = strsplit(tmpText1{7},'[');
        tmpx = str2double(tmpText2{2});
        tmpText2 = strsplit(tmpText1{8},']');
        tmpy = str2double(tmpText2{1});
        trialInfo(trialNo).startSearchLoc = [tmpx tmpy];
        trialInfo(trialNo).startSearchLoc_adjust = trialInfo(trialNo).posXY_startToDeliver_dropTarget(1,:);
        
        tmpText2 = strsplit(tmpText1{11},'(');
        tmpx = str2double(tmpText2{2});
        tmpy = str2double(tmpText1{12});
        trialInfo(trialNo).dropLoc = [tmpx tmpy];
    end
    fclose(fid);
    behInfo(subjNo).trialInfo = trialInfo;
    toc
end
% behInfo(3) = [];
save([saveDataPath,'behInfo_ifado'],'behInfo');


saveDataPath = 'F:\project\navigationParadigm\MorrisMaze\analyzeBehData\';
load([saveDataPath,'behInfo_ifado'])

for subjNo = 1:length(behInfo)
    fileID = fopen(strcat(saveDataPath,behInfo(subjNo).subjID,'.txt'),'w');
    fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
        'trialNo','trialType','startTime','collectTargetTime','targetLoc_x','targetLoc_y',...
        'instrToFindTargetTime','startSearchLoc_x','startSearchLoc_y',...
        'dropTargetTime','dropLoc_x','dropLoc_y','targetPic','picMemory')
    for trialNo = 1:length(behInfo(subjNo).trialInfo)
        trialType = behInfo(subjNo).trialInfo(trialNo).trialType;
        startTime = behInfo(subjNo).trialInfo(trialNo).startTime;
        
        collectTargetTime = behInfo(subjNo).trialInfo(trialNo).collectTargetTime;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        
        instrToFindTargetTime = behInfo(subjNo).trialInfo(trialNo).instrToFindTargetTime;
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        
        dropTargetTime = behInfo(subjNo).trialInfo(trialNo).dropTargetTime;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;
        
        targetPic = behInfo(subjNo).trialInfo(trialNo).targetPic;
        picMemory = behInfo(subjNo).trialInfo(trialNo).picMemory;
        fprintf(fileID,'%3d, %s, %s, %s, %2.5d, %2.5d, %s, %2.5d, %2.5d, %s, %2.5d, %2.5d, %s, %1d\n',...
                    trialNo,trialType,num2str(startTime),num2str(collectTargetTime),targetLoc(1),targetLoc(2),...
                    num2str(instrToFindTargetTime),startSearchLoc(1),startSearchLoc(2),num2str(dropTargetTime),...
                    dropLoc(1),dropLoc(2),targetPic,picMemory);
    end
    fclose(fileID);
end




subjNoAll = length(behInfo);
figure;
subjData = zeros(subjNoAll,5);
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    
    trialdata = zeros(trialNoAll,4);
    for trialNo = 1:70%trialNoAll
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;
        
        pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
        dropError = pdist([targetLoc;dropLoc],'euclidean');
        travelLenth = pdist([startSearchLoc;dropLoc],'euclidean');
        trialdata(trialNo,1) = pathLenth;
        trialdata(trialNo,2) = dropError;%dropError;
        trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
        trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
        trialdata(trialNo,5) = (behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ...
            behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime)*dropError/pathLenth;
    end
    tmpidx = find(trialdata(:,4)==1);
    subplot(3,ceil(subjNoAll/3),subjNo);
    plot(trialdata(tmpidx,1),trialdata(tmpidx,2),'.')
    
    tmpItem = trialdata(tmpidx,1);
    tmpSpac = trialdata(tmpidx,2);
    [r,pr] = corr(tmpItem,tmpSpac,'type','spearman');
    p = polyfit(tmpItem,tmpSpac,1);
    tmpx = linspace(min(tmpItem),max(tmpItem),10);
    tmpy = polyval(p,tmpx);
    plot(tmpItem,tmpSpac,'ks','MarkerFaceColor',[0 .5 .5],'MarkerEdgeColor',[0 .5 .5]); hold on
    plot(tmpx,tmpy,'k-','LineWidth',3);
    title(num2str(r))
    ylim([0 20])
end



subjNoAll = length(behInfo);
subjData = NaN(subjNoAll,5);
ageInfo = NaN(subjNoAll,1);
geneInfo = NaN(subjNoAll,1);
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    
    if trialNoAll>70
        ageInfo(subjNo) = behInfo(subjNo).age;
        geneInfo(subjNo) = behInfo(subjNo).geneType2;
        trialNoAll = 70;
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
    
%     subplot(3,ceil(subjNoAll/3),subjNo);
    tmpidx = find(trialdata(:,4)==1);
%     plot(trialdata(tmpidx,2),'g-');
%     
    tmpidx1 = find(trialdata(tmpidx,3)==1);
%     hold on
%     plot(tmpidx1,trialdata(tmpidx(tmpidx1),2),'r*');
    subjData(subjNo,1) = nanmean(trialdata(tmpidx,2));
    subjData(subjNo,2) = length(tmpidx1)/length(tmpidx);
%     
    tmpidx = find(trialdata(:,4)==0);
%     plot(trialdata(tmpidx,2),'b-');
    tmpidx1 = find(trialdata(tmpidx,3)==1);
%     hold on
%     plot(tmpidx1,trialdata(tmpidx(tmpidx1),2),'r*');
    subjData(subjNo,3) = nanmean(trialdata(tmpidx,2));
    subjData(subjNo,4) = length(tmpidx1)/length(tmpidx);
%     title(behInfo(subjNo).subjID);
    
%     ylim([0 25])
    end
end
% suptitle('blue: teleteport')


subjNoAll = length(behInfo);
subjData = NaN(subjNoAll,5);
ageInfo = NaN(subjNoAll,1);
geneInfo = NaN(subjNoAll,1);
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);

    ageInfo(subjNo) = behInfo(subjNo).age;
    geneInfo(subjNo) = behInfo(subjNo).geneType2;
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

ageInfo1 = ageInfo(~isnan(ageInfo));
medianData = median(ageInfo1);
subjIdx1 = find(ageInfo>=medianData); nameIdx = 1;
subjIdx2 = find(ageInfo<medianData);
% subjIdx1 = find(geneInfo==2); nameIdx = 2;
% subjIdx2 = find(geneInfo<2);
subjIdx3 = 1:subjNoAll;
subjIdxPool = {'subjIdx1','subjIdx2','subjIdx3'};
carrierAge = ageInfo(geneInfo<2);
nonCarrAge = ageInfo(geneInfo==2);
[~,p,~,stats]= ttest2(carrierAge,nonCarrAge);
figure;
bar(1:2,[nanmean(carrierAge) nanmean(nonCarrAge)],0.4,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5); hold on
errorbar(1:2,[nanmean(carrierAge) nanmean(nonCarrAge)],[nanstd(carrierAge)/sqrt(length(carrierAge)) nanstd(nonCarrAge)/sqrt(length(nonCarrAge))],'.','LineWidth',3, 'color',[0 .9 .9]); hold on
plot(ones(size(carrierAge)),carrierAge,'o'); hold on
plot(2*ones(size(nonCarrAge)),nonCarrAge,'o');
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'carrier','non-carrier'});
title(strcat('age p=',num2str(p),' carrierNo=',num2str(length(carrierAge))));

figure;
for i = 1:2
    subjIdx = eval(subjIdxPool{i});
    [~,p,~,stats] = ttest(subjData(subjIdx,1),subjData(subjIdx,3));
    meanD = squeeze(nanmean(subjData(subjIdx,:)));
    stdD = squeeze(nanstd(subjData(subjIdx,:)));
    subplot(1,2,i);
    bar(1:2,meanD([1 3]),0.4,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5); hold on
    errorbar(1:2,meanD([1 3]),stdD([1 3]),'.','LineWidth',3, 'color',[0 .9 .9]); hold on
    for subjNo = subjIdx
        plot(1:2,[subjData(subjNo,1) subjData(subjNo,3)],'-','LineWidth',2); hold on
    end
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'lure','teleport'});
    title({'drop error' strcat('p=',num2str(p))});
end
[~,p1,~,stats] = ttest2(subjData(subjIdx1,1),subjData(subjIdx2,1));
[~,p2,~,stats] = ttest2(subjData(subjIdx1,3),subjData(subjIdx2,3));
if nameIdx ==1
    suptitle(strcat('old vs. young; Lure p=',num2str(p1),'teleport p=',num2str(p2)));
elseif nameIdx ==2
    suptitle(strcat('noGene vs. Gene; Lure p=',num2str(p1),'teleport p=',num2str(p2)));
end


figure;
for i = 1:2
    subjIdx = eval(subjIdxPool{i});
    [~,p,~,stats] = ttest(subjData(subjIdx,2),subjData(subjIdx,4));
    meanD = squeeze(nanmean(subjData(subjIdx,:)));
    stdD = squeeze(nanstd(subjData(subjIdx,:)));
    subplot(1,2,i);
    bar(1:2,meanD([2 4]),0.4,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5);hold on
    errorbar(1:2,meanD([2 4]),stdD([2 4]),'.','LineWidth',3, 'color',[0 .9 .9]);hold on
    for subjNo = subjIdx
        plot(1:2,[subjData(subjNo,2) subjData(subjNo,4)],'-','LineWidth',2); hold on
    end
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'lure','teleport'});
    title({'item memory' strcat('p=',num2str(p))});
end
[~,p1,~,stats] = ttest2(subjData(subjIdx1,2),subjData(subjIdx2,2));
[~,p2,~,stats] = ttest2(subjData(subjIdx1,4),subjData(subjIdx2,4));
if nameIdx ==1
    suptitle(strcat('old vs. young; Lure p=',num2str(p1),'teleport p=',num2str(p2)));
elseif nameIdx ==2
    suptitle(strcat('noGene vs. Gene; Lure p=',num2str(p1),'teleport p=',num2str(p2)));
end
%% check the correlation between spatial memory and item memory across participants
% subjIdx1 = find(ageInfo>=medianData); nameIdx = 1;
% subjIdx2 = find(ageInfo<medianData);
subjIdx1 = find(geneInfo==2); nameIdx = 2;
subjIdx2 = find(geneInfo<2);
figure
for i = 1:2
    subjIdx = eval(subjIdxPool{i});
    tmpItem = subjData(subjIdx,2);
    tmpSpac = subjData(subjIdx,1);
    tmpItem(isnan(tmpItem)) = [];
    tmpSpac(isnan(tmpSpac)) = [];
    [r,pr] = corr(tmpItem,tmpSpac,'type','spearman');
    p = polyfit(tmpItem,tmpSpac,1);
    tmpx = linspace(min(tmpItem),max(tmpItem),10);
    tmpy = polyval(p,tmpx);
    subplot(2,2,(i-1)*2+1);
    plot(tmpItem,tmpSpac,'ks','MarkerFaceColor',[0 .5 .5],'MarkerEdgeColor',[0 .5 .5]); hold on
    plot(tmpx,tmpy,'k-','LineWidth',3);
    title(strcat('lure r=',num2str(r),' p=',num2str(pr)))
    xlabel('item memory')
    ylabel('drop error')

    tmpItem = subjData(subjIdx,4);
    tmpSpac = subjData(subjIdx,3);
    tmpItem(isnan(tmpItem)) = [];
    tmpSpac(isnan(tmpSpac)) = [];
    [r,pr] = corr(tmpItem,tmpSpac,'type','spearman');
    p = polyfit(tmpItem,tmpSpac,1);
    tmpx = linspace(min(tmpItem),max(tmpItem),10);
    tmpy = polyval(p,tmpx);
    subplot(2,2,(i-1)*2+2);
    plot(tmpItem,tmpSpac,'ks','MarkerFaceColor',[0 .5 .5],'MarkerEdgeColor',[0 .5 .5]); hold on
    plot(tmpx,tmpy,'k-','LineWidth',3);
    title(strcat('teleport r=',num2str(r),' p=',num2str(pr)))
    xlabel('item memory')
    ylabel('drop error')
end
suptitle('correlation across participants')

%% relation of item memory and spatial memory
subjData = NaN(subjNoAll,2,2);
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    trialdata = zeros(trialNoAll,4);
    for trialNo = 1:trialNoAll
        startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
        targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;

        pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
        dropError = pdist([targetLoc;dropLoc],'euclidean');
        trialdata(trialNo,1) = pathLenth;
        trialdata(trialNo,2) = dropError;
        trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
        trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
    end
    tmpidx = find(trialdata(:,3)==1 & trialdata(:,4)==1);
    subjData(subjNo,1,1) = mean(trialdata(tmpidx,2));
    tmpidx = find(trialdata(:,3)==1 & trialdata(:,4)==0);
    subjData(subjNo,1,2) = mean(trialdata(tmpidx,2));
    tmpidx = find(trialdata(:,3)==0 & trialdata(:,4)==1);
    subjData(subjNo,2,1) = mean(trialdata(tmpidx,2));
    tmpidx = find(trialdata(:,3)==0 & trialdata(:,4)==0);
    subjData(subjNo,2,2) = mean(trialdata(tmpidx,2));
end

meanD = squeeze(nanmean(subjData));
stdD = squeeze(nanstd(subjData));

figure;
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
    for subjNo = 1:subjNoAll
        plot(x,[subjData(subjNo,1,i) subjData(subjNo,2,i)],'-','LineWidth',2,'HandleVisibility','off'); hold on
    end
end
hold off
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'correct','incorrect'});
legend({'lure','teleport'})
ylabel('drop error')
set(gca,'fontsize',14)




%% distribution of trials
figure;
for subjNo = 1:subjNoAll
    subplot(3,ceil(subjNoAll/3),subjNo);
%     viscircles([0,0],17,'Color','k','LineWidth',3); hold on;
    trialNoAll = length(behInfo(subjNo).trialInfo);
    targetLocAll = zeros(trialNoAll,2);
    startSearchLocAll = zeros(trialNoAll,2);
    dropLocAll = zeros(trialNoAll,2);
    for trialNo = 1:trialNoAll
        posXY_start_collectTarget = behInfo(subjNo).trialInfo(trialNo).posXY_start_collectTarget;
        plot(posXY_start_collectTarget(:,1),posXY_start_collectTarget(:,2),'k-','lineWidth',0.5); hold on;
        if strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure')
            posXY_collectTarget_collectLure = behInfo(subjNo).trialInfo(trialNo).posXY_collectTarget_collectLure;
            plot(posXY_collectTarget_collectLure(:,1),posXY_collectTarget_collectLure(:,2),'g-','lineWidth',0.5); hold on;
        else
            tmpx = linspace(posXY_start_collectTarget(end,1),behInfo(subjNo).trialInfo(trialNo).startSearchLoc(1),10);
            tmpy = linspace(posXY_start_collectTarget(end,2),behInfo(subjNo).trialInfo(trialNo).startSearchLoc(2),10);
            plot(tmpx,tmpy,'g--'); hold on
        end
        
        posXY_startToDeliver_dropTarget = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
        plot(posXY_startToDeliver_dropTarget(:,1),posXY_startToDeliver_dropTarget(:,2),'b-','lineWidth',0.5); hold on;
        targetLocAll(trialNo,:) = behInfo(subjNo).trialInfo(trialNo).targetLoc;
        dropLocAll(trialNo,:) = behInfo(subjNo).trialInfo(trialNo).dropLoc;
        startSearchLocAll(trialNo,:) = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
    end
    plot(targetLocAll(:,1),targetLocAll(:,2),'rd','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r'); hold on
    plot(startSearchLocAll(:,1),startSearchLocAll(:,2),'g*','MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor','g'); hold on
    plot(dropLocAll(:,1),dropLocAll(:,2),'b^','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b'); hold on
    xlim([-20 20]);
    ylim([-20 20]);
    axis('square')
    
    title(behInfo(subjNo).subjID);
end

%% check facing direction of each trial
% relay test
figure;
for subjNo = 1:subjNoAll
    subplot(3,ceil(subjNoAll/3),subjNo);
    trialNoAll = length(behInfo(subjNo).trialInfo);
    
    posX = [];
    posY = [];
    for trialNo = 1:trialNoAll
        posXY_start_collectTarget = behInfo(subjNo).trialInfo(trialNo).posXY_start_collectTarget;
        posX = cat(1,posX,posXY_start_collectTarget(:,1));
        posY = cat(1,posY,posXY_start_collectTarget(:,2));
        if strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure')
            posXY_collectTarget_collectLure = behInfo(subjNo).trialInfo(trialNo).posXY_collectTarget_collectLure;
            posX = cat(1,posX,posXY_collectTarget_collectLure(:,1));
            posY = cat(1,posY,posXY_collectTarget_collectLure(:,2));
        end
        
        posXY_startToDeliver_dropTarget = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
        posX = cat(1,posX,posXY_startToDeliver_dropTarget(:,1));
        posY = cat(1,posY,posXY_startToDeliver_dropTarget(:,2));
    end
    
    tmpx = diff(posX);
    tmpy = diff(posY);
    idx = find(tmpx == 0 & tmpy == 0);
    tmpx(idx) = [];
    tmpy(idx) = [];
    tmpSlop = atan2(tmpy,tmpx);
    histogram(tmpSlop,linspace(-pi,pi,20));
    axis('square')
    title(behInfo(subjNo).subjID);
end


figure;
posX = [];
posY = [];
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    for trialNo = 1:trialNoAll
        posXY_start_collectTarget = behInfo(subjNo).trialInfo(trialNo).posXY_start_collectTarget;
        posX = cat(1,posX,posXY_start_collectTarget(:,1));
        posY = cat(1,posY,posXY_start_collectTarget(:,2));
        if strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure')
            posXY_collectTarget_collectLure = behInfo(subjNo).trialInfo(trialNo).posXY_collectTarget_collectLure;
            posX = cat(1,posX,posXY_collectTarget_collectLure(:,1));
            posY = cat(1,posY,posXY_collectTarget_collectLure(:,2));
        end
        
        posXY_startToDeliver_dropTarget = behInfo(subjNo).trialInfo(trialNo).posXY_startToDeliver_dropTarget;
        posX = cat(1,posX,posXY_startToDeliver_dropTarget(:,1));
        posY = cat(1,posY,posXY_startToDeliver_dropTarget(:,2));
    end
end
    
tmpx = diff(posX);
tmpy = diff(posY);
idx = find(tmpx == 0 & tmpy == 0);
tmpx(idx) = [];
tmpy(idx) = [];
tmpSlop = atan2(tmpy,tmpx);
histogram(tmpSlop,linspace(-pi,pi,20));
axis('square')
title('all patients');


%% performance of each picture
subjNoAll = length(behInfo);
pictureInfo = NaN(400,subjNoAll);
for subjNo = 1:subjNoAll
    trialNoAll = length(behInfo(subjNo).trialInfo);
    for trialNo = 1:trialNoAll
        picNo = str2double(behInfo(subjNo).trialInfo(trialNo).targetPic(1:end-1));
        pictureInfo(picNo,subjNo) = behInfo(subjNo).trialInfo(trialNo).picMemory;
    end
end
picAll = mean(pictureInfo,2);
tmpidx = find(picAll==1);
figure;
plot(1:length(picAll),picAll,'s','MarkerFaceColor',[0 .5 .5],'MarkerEdgeColor',[0 .5 .5])
xlabel('pictureNo')
ylabel('mean accuracy across participants')
title(strcat(num2str(length(tmpidx)),' pictures without wrong answers'));



%% gender effect
gendInfo = {behInfo.sex};
gendPool = {'f','m'};
figure;
for gendNo = 1:length(gendPool)
    subjIdx = find(ismember(gendInfo,gendPool{gendNo}));
    subjData = zeros(length(subjIdx),5);
    count = 0;
    for subjNo = subjIdx
        count = count + 1;
        trialNoAll = length(behInfo(subjNo).trialInfo);
        trialdata = zeros(trialNoAll,4);
        for trialNo = 1:trialNoAll
            startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
            targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
            dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;

            pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
            dropError = pdist([targetLoc;dropLoc],'euclidean');
            trialdata(trialNo,1) = pathLenth;
            trialdata(trialNo,2) = dropError;%/pathLenth;%dropError;
            trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
            trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
            trialdata(trialNo,5) = (behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ...
                behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime)*dropError/pathLenth;
        end
        tmpidx = find(trialdata(:,4)==1);
        tmpidx1 = find(trialdata(tmpidx,3)==1);
        subjData(count,1) = nanmean(trialdata(tmpidx,2));
        subjData(count,2) = length(tmpidx1)/length(tmpidx);

        tmpidx = find(trialdata(:,4)==0);
        tmpidx1 = find(trialdata(tmpidx,3)==1);
        subjData(count,3) = nanmean(trialdata(tmpidx,2));
        subjData(count,4) = length(tmpidx1)/length(tmpidx);
    end

    [~,p,~,stats] = ttest(subjData(:,1),subjData(:,3));
    meanD = squeeze(mean(subjData));
    stdD = squeeze(std(subjData));
    
    subplot(2,2,(gendNo-1)*2+1);
    bar(1:2,meanD([1 3]),0.4,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5); hold on
    errorbar(1:2,meanD([1 3]),stdD([1 3]),'.','LineWidth',3, 'color',[0 .9 .9]); hold on
    for subjNo = 1:length(subjIdx)
        plot(1:2,[subjData(subjNo,1) subjData(subjNo,3)],'-','LineWidth',2); hold on
    end
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'lure','teleport'});
    title({gendPool{gendNo}, 'drop error' strcat('p=',num2str(p))});


    [~,p,~,stats] = ttest(subjData(:,2),subjData(:,4));
    subplot(2,2,(gendNo-1)*2+2);
    bar(1:2,meanD([2 4]),0.4,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5);hold on
    errorbar(1:2,meanD([2 4]),stdD([2 4]),'.','LineWidth',3, 'color',[0 .9 .9]);hold on
    for subjNo = 1:length(subjIdx)
        plot(1:2,[subjData(subjNo,2) subjData(subjNo,4)],'-','LineWidth',2); hold on
    end
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',{'lure','teleport'});
    title({gendPool{gendNo} 'item memory' strcat('p=',num2str(p))});
end


gendInfo = {behInfo.sex};
gendPool = {'f','m'};
subjDataAll = cell(1,length(gendPool));
for gendNo = 1:length(gendPool)
    subjIdx = find(ismember(gendInfo,gendPool{gendNo}));
    subjData = zeros(length(subjIdx),5);
    count = 0;
    for subjNo = subjIdx
        count = count + 1;
        trialNoAll = length(behInfo(subjNo).trialInfo);
        trialdata = zeros(trialNoAll,4);
        for trialNo = 1:trialNoAll
            startSearchLoc = behInfo(subjNo).trialInfo(trialNo).startSearchLoc;
            targetLoc = behInfo(subjNo).trialInfo(trialNo).targetLoc;
            dropLoc = behInfo(subjNo).trialInfo(trialNo).dropLoc;

            pathLenth = pdist([targetLoc;startSearchLoc],'euclidean');
            dropError = pdist([targetLoc;dropLoc],'euclidean');
            trialdata(trialNo,1) = pathLenth;
            trialdata(trialNo,2) = dropError;%/pathLenth;%dropError;
            trialdata(trialNo,3) = behInfo(subjNo).trialInfo(trialNo).picMemory;
            trialdata(trialNo,4) = strcmp(behInfo(subjNo).trialInfo(trialNo).trialType,'lure');
            trialdata(trialNo,5) = (behInfo(subjNo).trialInfo(trialNo).dropTargetTime - ...
                behInfo(subjNo).trialInfo(trialNo).startToDeliverTargetTime)*dropError/pathLenth;
        end
        tmpidx = find(trialdata(:,4)==1);
        tmpidx1 = find(trialdata(tmpidx,3)==1);
        subjData(count,1) = nanmean(trialdata(tmpidx,2));
        subjData(count,2) = length(tmpidx1)/length(tmpidx);

        tmpidx = find(trialdata(:,4)==0);
        tmpidx1 = find(trialdata(tmpidx,3)==1);
        subjData(count,3) = nanmean(trialdata(tmpidx,2));
        subjData(count,4) = length(tmpidx1)/length(tmpidx);
    end
    subjDataAll{gendNo} = subjData(:,[1 3]);
end

meanD = zeros(2,2);
stdD = zeros(2,2);
meanD(1,:) = squeeze(mean(subjDataAll{1}));
meanD(2,:) = squeeze(mean(subjDataAll{2}));
stdD(1,:) = squeeze(std(subjDataAll{1}));
stdD(2,:) = squeeze(std(subjDataAll{2}));

figure;
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

i = 1;
x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
for i=1:ngroups
    x1 = [x(i) i+i-x(i)];
    for subjNo = 1:length(subjDataAll{i})
        plot(x1,[subjDataAll{i}(subjNo,1) subjDataAll{i}(subjNo,2)],'-','LineWidth',2,'HandleVisibility','off'); hold on
    end
end
hold off
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'famale','male'});
legend({'lure','teleport'})
ylabel('drop error')
set(gca,'fontsize',14)

tmpF = subjDataAll{1}(:,1) - subjDataAll{1}(:,2);
tmpM = subjDataAll{2}(:,1) - subjDataAll{2}(:,2);
[~,p,~,stat] = ttest2(tmpF,tmpM);
titleStr = strcat('interaction p=',num2str(p));

tmpF = mean(subjDataAll{1},2);
tmpM = mean(subjDataAll{2},2);
[~,p,~,stat] = ttest2(tmpF,tmpM);
titleStr = strcat(titleStr, ' gender p=',num2str(p));

title(titleStr);


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

for conditionNo = 1:length(conditionLabel)
    tmp_label = conditionLabel{conditionNo};
    tmpConNo = conditionNoLabel(conditionNo);
%     figure;
%     for subjNo = 1:subjNoAll
%         if ~isempty(dataAll(subjNo).trialdata)
%             tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo);
%             r = corr(dataAll(subjNo).trialdata(tmpidx1,5),dataAll(subjNo).trialdata(tmpidx1,6),'type','Pearson');
%             subplot(3,5,subjNo);
%             plot(dataAll(subjNo).trialdata(tmpidx1,5),dataAll(subjNo).trialdata(tmpidx1,6),'k.');
%             xlabel('incoming path_ideal');
%             ylabel('incoming time');
%             title(num2str(r))
%         end
%     end
%     suptitle(tmp_label)

%     figure;
%     for subjNo = 1:subjNoAll
%         tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo);
%         r = corr(dataAll(subjNo).trialdata(tmpidx1,5),dataAll(subjNo).trialdata(tmpidx1,2),'type','Pearson');
%         subplot(3,5,subjNo);
%         plot(dataAll(subjNo).trialdata(tmpidx1,5),dataAll(subjNo).trialdata(tmpidx1,2),'k.');
%         xlabel('incoming path_ideal');
%         ylabel('drop error');
%         title(num2str(r))
%     end
%     suptitle(tmp_label)

%     figure;
%     for subjNo = 1:subjNoAll
%         if ~isempty(dataAll(subjNo).trialdata)
%             tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo);
%             r = corr(dataAll(subjNo).trialdata(tmpidx1,6),dataAll(subjNo).trialdata(tmpidx1,2),'type','Pearson');
%             subplot(3,5,subjNo);
%             plot(dataAll(subjNo).trialdata(tmpidx1,6),dataAll(subjNo).trialdata(tmpidx1,2),'k.');
%             xlabel('incoming time');
%             ylabel('drop error');
%             title(num2str(r))
%         end
%     end
%     suptitle(tmp_label)

    
    tmpr = NaN(1,subjNoAll);
    for subjNo = 1:subjNoAll
        if ~isempty(dataAll(subjNo).trialdata)
            tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo);
            dropError = dataAll(subjNo).trialdata(tmpidx1,2);
            incomingDist = dataAll(subjNo).trialdata(tmpidx1,5);
            incomingTime = dataAll(subjNo).trialdata(tmpidx1,6);
            tmpr(subjNo) = partialcorr(dropError,incomingDist,incomingTime);
        end
    end
    subjIdx1 = find(geneInfo==2); nameIdx = 2;
    subjIdx2 = find(geneInfo<2);
    figure
    for i = 1:2
        subplot(2,2,i);
        subjIdx = eval(subjIdxPool{i});
        bar(1,nanmean(tmpr(subjIdx))); hold on
        plot(ones(1,length(subjIdx)),tmpr(subjIdx),'.','markerSize',30)
        [~,p,~,stat] = ttest(tmpr(subjIdx));
        title(strcat(tmp_label,'-partial correlation-',num2str(p)));
    end


    dropErrorDis_itemMem = NaN(subjNoAll,4,2);
    for subjNo = 1:subjNoAll   
        if ~isempty(dataAll(subjNo).trialdata)
    %         tmpidx1 = find(dataAll(subjNo).trialdata(:,3)==1);
            tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo & dataAll(subjNo).trialdata(:,3)==1);
            dropErrorDis_itemMem(subjNo,1,1) = nanmean(dataAll(subjNo).trialdata(tmpidx1,6));
            dropErrorDis_itemMem(subjNo,2,1) = nanmean(dataAll(subjNo).trialdata(tmpidx1,5));
            dropErrorDis_itemMem(subjNo,3,1) = nanmean(dataAll(subjNo).trialdata(tmpidx1,7));
            dropErrorDis_itemMem(subjNo,4,1) = nanmean(dataAll(subjNo).trialdata(tmpidx1,2));
            tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo & dataAll(subjNo).trialdata(:,3)==0);
    %         tmpidx1 = find(dataAll(subjNo).trialdata(:,3)==0);
            dropErrorDis_itemMem(subjNo,1,2) = nanmean(dataAll(subjNo).trialdata(tmpidx1,6));
            dropErrorDis_itemMem(subjNo,2,2) = nanmean(dataAll(subjNo).trialdata(tmpidx1,5));
            dropErrorDis_itemMem(subjNo,3,2) = nanmean(dataAll(subjNo).trialdata(tmpidx1,7));
            dropErrorDis_itemMem(subjNo,4,2) = nanmean(dataAll(subjNo).trialdata(tmpidx1,2));
        end
    end
    


    % incoming time
    indexPool = {'incoming time','incoming dist_ideal','incoming dist_prac','drop error'};
    figure;
    for idxNo = 1:length(indexPool)
        tmpData = squeeze(dropErrorDis_itemMem(:,idxNo,:));
        [~,p,~,stat] = ttest(tmpData(:,1),tmpData(:,2));
        meanData = nanmean(tmpData,1);
        stdData = nanstd(tmpData)/sqrt(subjNoAll);
        subplot(2,2,idxNo);
        bar(1:2,meanData); hold on
        errorbar(1:2,meanData,stdData,'k.');
        set(gca,'xtick',[1 2]);
        set(gca,'xticklabel',{'Remember','Forgotten'});
        xlabel('item memory');
        ylabel(indexPool{idxNo});
        title(num2str(p));
    end
    suptitle(tmp_label)
    
end


tmpr = NaN(subjNoAll,2);
for conditionNo = 1:length(conditionLabel)
    tmpConNo = conditionNoLabel(conditionNo);
    for subjNo = 1:subjNoAll
        if ~isempty(dataAll(subjNo).trialdata)
            tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo);
            tmpr(subjNo,conditionNo) = corr(dataAll(subjNo).trialdata(tmpidx1,2),dataAll(subjNo).trialdata(tmpidx1,5),'type','Pearson');
    %         tmpr(subjNo,conditionNo) = partialcorr(dataAll(subjNo).trialdata(tmpidx1,1),...
    %             dataAll(subjNo).trialdata(tmpidx1,2),...
    %             dataAll(subjNo).trialdata(tmpidx1,5),'type','Pearson');
        end
    end
end

[~,p,~,stat] = ttest(tmpr(:,1),tmpr(:,2));
figure;
meanData = nanmean(tmpr);
stdData=nanstd(tmpr)/sqrt(subjNoAll);
bar(1:2,meanData);
hold on
errorbar(1:2,meanData,stdData,'k.');
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',conditionLabel);
ylabel('corr between travel path and drop error partial out ideal incoming path length')
% ylabel('corr between travel path and drop error')
% ylabel('corr between travel path and incoming pathlen')
% ylabel('corr between drop error and incoming pathlen partial out travel length')
% ylabel('corr between drop error and incoming pathlen')
title(strcat('p=',num2str(p)))


% compare the ideal path length between conditions
tmpr = NaN(subjNoAll,2);
for conditionNo = 1:length(conditionLabel)
    tmpConNo = conditionNoLabel(conditionNo);
    for subjNo = 1:subjNoAll
        if ~isempty(dataAll(subjNo).trialdata)
            tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo);
            tmpr(subjNo,conditionNo) = nanmean(dataAll(subjNo).trialdata(tmpidx1,5));
        end
    end
end
[~,p,~,stat] = ttest(tmpr(:,1),tmpr(:,2));
figure;
meanData = nanmean(tmpr);
stdData=nanstd(tmpr)/sqrt(subjNoAll);
bar(1:2,meanData);
hold on
errorbar(1:2,meanData,stdData,'k.');
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',conditionLabel);
ylabel('ideal incoming path length')
title(strcat('p=',num2str(p)))

%% corelation between item memory performance and spatial memory performance across participants

for conditionNo = 1:length(conditionLabel)
    tmp_label = conditionLabel{conditionNo};
    tmpConNo = conditionNoLabel(conditionNo);
    
    dropErrorDis_itemMem = NaN(subjNoAll,5);
    for subjNo = 1:subjNoAll   
        if ~isempty(dataAll(subjNo).trialdata)
    %         tmpidx1 = find(dataAll(subjNo).trialdata(:,3)==1);
            tmpidx1 = find(dataAll(subjNo).trialdata(:,4)==tmpConNo);
            dropErrorDis_itemMem(subjNo,1) = nanmean(dataAll(subjNo).trialdata(tmpidx1,6));
            dropErrorDis_itemMem(subjNo,2) = nanmean(dataAll(subjNo).trialdata(tmpidx1,5));
            dropErrorDis_itemMem(subjNo,3) = nanmean(dataAll(subjNo).trialdata(tmpidx1,7));
            dropErrorDis_itemMem(subjNo,4) = nanmean(dataAll(subjNo).trialdata(tmpidx1,2));
            dropErrorDis_itemMem(subjNo,5) = (nansum(dataAll(subjNo).trialdata(tmpidx1,3)))/length(tmpidx1);
        end
    end
    
    % incoming time
    indexPool = {'incoming time','incoming dist_ideal','incoming dist_prac','drop error'};
    
    
%     subjIdx1 = find(ageInfo>=medianData); nameIdx = 1;
%     subjIdx2 = find(ageInfo<medianData);
    subjIdx1 = find(geneInfo==2); nameIdx = 2;
    subjIdx2 = find(geneInfo<2);
    if nameIdx == 2
        condiName = {'non-carrier','carrier'};
    elseif nameIdx == 1
        condiName = {'old','young'};
    end
    for i = 1:2
        subjIdx = eval(subjIdxPool{i});
        itemMemory = dropErrorDis_itemMem(:,5);
        itemMemory(isnan(itemMemory)) = [];
        itemMemory1 = itemMemory(subjIdx);
        figure;
        for idxNo = 1:length(indexPool)
            tmpData = squeeze(dropErrorDis_itemMem(subjIdx,idxNo));
            tmpData(isnan(tmpData)) = [];
            [r,p] = corr(tmpData,itemMemory1,'type','Pearson');
            subplot(2,2,idxNo);
            plot(itemMemory1,tmpData,'k.','markerSize',30);
            xlabel('item memory');
            ylabel(indexPool{idxNo});
            title(strcat('r=',num2str(r),' p=',num2str(p)));
        end
        suptitle(strcat('corr between ## and ## across participants ',tmp_label,condiName{i}));


        % partial out time
        indexPool = {'incoming time','incoming dist_ideal','incoming dist_prac','drop error'};
        itemMemory = dropErrorDis_itemMem(subjIdx,5);
        incomingTime =  dropErrorDis_itemMem(subjIdx,1);
        itemMemory(isnan(itemMemory)) = [];
        incomingTime(isnan(incomingTime)) = [];
        figure;
        for idxNo = 2:length(indexPool)
            tmpData = squeeze(dropErrorDis_itemMem(subjIdx,idxNo));
            tmpData(isnan(tmpData)) = [];
            [r,p] = partialcorr(tmpData,itemMemory,incomingTime);
            subplot(2,2,idxNo);
            plot(itemMemory,tmpData,'k.','markerSize',30);
            xlabel('item memory');
            ylabel(indexPool{idxNo});
            title(strcat('r=',num2str(r),' p=',num2str(p)));
        end
        suptitle(strcat('partial out incoming time across participants  ',tmp_label,condiName{i}))
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

figure;
for subjNo = 1:subjNoAll
    r = corr(dataAll(subjNo).trialdata(:,5),dataAll(subjNo).trialdata(:,2),'type','Pearson');
    subplot(4,4,subjNo);
    plot(dataAll(subjNo).trialdata(:,5),dataAll(subjNo).trialdata(:,2),'k.');
    xlabel('distance to the nearest landmark');
    ylabel('drop error');
    title(num2str(r))
end


% subjIdx1 = find(ageInfo>=medianData); nameIdx = 1;
% subjIdx2 = find(ageInfo<medianData);
subjIdx1 = find(geneInfo==2); nameIdx = 2;
subjIdx2 = find(geneInfo==1);
if nameIdx == 2
    condiName = {'non-carrier','carrier'};
elseif nameIdx == 1
    condiName = {'old','young'};
end
figure;
dataTemp = cell(2,2);
for condiNo = 1:2
    subplot(1,2,condiNo);
    for i = 1:2
        subjIdx = eval(subjIdxPool{i});
        tmpCorr = NaN(1,length(subjIdx));
        for subjNo = 1:length(subjIdx)
            if ~isempty(dataAll(subjIdx(subjNo)).trialdata)
                tmpTrialidx = find(dataAll(subjIdx(subjNo)).trialdata(:,4)==((condiNo-1)*1));
                tmpCorr(subjNo) = corr(dataAll(subjIdx(subjNo)).trialdata(tmpTrialidx,5),dataAll(subjIdx(subjNo)).trialdata(tmpTrialidx,2),'type','Pearson');
            end
        end
        dataTemp{condiNo,i} = tmpCorr;
        [~,p,~,stat] = ttest(tmpCorr);
        bar(i,nanmean(tmpCorr)); hold on
        errorbar(i,nanmean(tmpCorr),nanstd(tmpCorr)/sqrt(length(subjIdx)),'k.'); hold on
        plot(ones(size(tmpCorr))*i,tmpCorr,'o');
        ylabel('drop error corr. dis to Nearest landmark');
        text(i,nanmean(tmpCorr),strcat('p=',num2str(p),' t=',num2str(stat.tstat)));
    end
    [~,p,~,stat] = ttest2(dataTemp{condiNo,1},dataTemp{condiNo,2});
    set(gca,'xtick',[1 2]);
    set(gca,'xticklabel',condiName);
    title(strcat(conditionLabel{condiNo},' p=',num2str(p),' t=',num2str(stat.tstat)));
end
[~,p1,~,stat1] = ttest(dataTemp{1,1},dataTemp{2,1});
[~,p2,~,stat2] = ttest(dataTemp{1,2},dataTemp{2,2});
suptitle({strcat(condiName{1},conditionLabel{1},' vs ',conditionLabel{2},' p=',num2str(p1),' t=',num2str(stat1.tstat)) ...
    strcat(condiName{2},conditionLabel{1},' vs ',conditionLabel{2},' p=',num2str(p2),' t=',num2str(stat2.tstat))})


%% calculate the navigatning pattern
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
            incomingXYTime = behInfo(subjNo).trialInfo(trialNo).posXYTime_startToDeliver_dropTarget;
            incomingSpeed1 = arrayfun(@(x) pdist([incomingXY(x,:);incomingXY(x+1,:)],'euclidean'),1:size(incomingXY,1)-1);
            incomingSpeed2 = arrayfun(@(x) incomingXYTime(x+1)-incomingXYTime(x),1:size(incomingXY,1)-1);
            incomingSpeed3 = arrayfun(@(x) incomingSpeed1(x+1)+incomingSpeed1(x),1:size(incomingSpeed1,2)-1);
            incomingSpeed4 = arrayfun(@(x) incomingSpeed2(x+1)+incomingSpeed2(x),1:size(incomingSpeed2,2)-1);
            incomingSpeed5 = incomingSpeed3./incomingSpeed4;
            
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
    if tmpAng<30*pi/180;
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
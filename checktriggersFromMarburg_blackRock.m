% extract photosensor triggers
clear all
addpath('H:\matlabCode\with&withoutpic\privateFunctions\');
addpath('F:\project\matlabCode\shadedErrorBar\');
addpath('F:\project\matlabCode\hui_code\');
addpath F:\project\VirtualMaze\matlabCode\nrlx2edf\;
addpath('F:\project\VirtualMaze\matlabCode\NPMK\');
% targetFolder = 'F:\project\navigationParadigm\MorrisMaze\data\Marburg\22Feb2021\iEEG\';
targetFolder = 'F:\project\navigationParadigm\MorrisMaze\data\Marburg\Test_Marburg_032021\blackrock\20210325-040034_island_trigg_both_cond_one_block\';

loadFileName = '20210325-040034-001';%'20210222-091919-001';
loadfileName = strcat(targetFolder,loadFileName,'.ns3');
dataNs = openNSx(loadfileName);
sr = 2000;%dataNs.MetaTags.TimeRes;
loadfileName = strcat(targetFolder,loadFileName,'.nev');
dataNev = openNEV(loadfileName,'nosave');

triggerData = double(dataNs.Data(129,:));

figure;plot(triggerData,'.')

minTrigger = min(triggerData);
maxTrigger = max(triggerData);
triggerData1 = diff(triggerData);
triggerData2 = zeros(size(triggerData1))';
% triggerData2(abs(triggerData1)<10) = 1;
triggerData2(triggerData<1000) = 1;
[clusterInfo,clusterNo] = extractCluster_1d(triggerData2);
triggerOnset_photosensor = NaN(clusterNo,1);
for tmpi = 1:clusterNo
    if min(triggerData(clusterInfo==tmpi))<840
        triggerOnset_photosensor(tmpi) = find(clusterInfo==tmpi,1,'last');
    end
end
triggerOnset_photosensor(isnan(triggerOnset_photosensor))=[];
timeDiff_eeg = round(diff(triggerOnset_photosensor)*1000/sr);
%% behavial data
% targetFolder = 'F:\project\navigationParadigm\MorrisMaze\data\Marburg\22Feb2021\beh\m_22022021_ID001_3\session_0\';
targetFolder = 'F:\project\navigationParadigm\MorrisMaze\data\Marburg\Test_Marburg_032021\behavior\morrisMaze_m_25032021_subjTest\session_0\';
loadFileName = 'log.txt';
fid = fopen(strcat(targetFolder,loadFileName));
fid_text = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',{'\t',',','(',')'});
% start and end time stamp
startTime = str2double(fid_text{1}(1));
endTime = str2double(fid_text{1}(end));
% trigger information
eventLineNo = find(ismember(fid_text{3},'EEGTrigger'));
eventName = fid_text{4}(eventLineNo);
eventTime = str2double(fid_text{1}(eventLineNo));

eventLineNo_bigEvent = find(ismember(fid_text{3},'bigEvent'));
eventName_bigEvent = fid_text{4}(eventLineNo_bigEvent);
eventTime_bigEvent = str2double(fid_text{1}(eventLineNo_bigEvent));


bigEvent_eventName = cell(size(eventName));
bigEvent_eventTime = NaN(size(eventTime));
for i = 1:length(eventName)
    lineNo = eventLineNo(i);
    eventLineNo_bigEvent1 = eventLineNo_bigEvent;
    eventLineNo_bigEvent1(eventLineNo_bigEvent1>lineNo) = [];
    bigEvent_eventName{i} = eventName_bigEvent{length(eventLineNo_bigEvent1)};
    bigEvent_eventTime(i) = eventTime_bigEvent(length(eventLineNo_bigEvent1));
end

% eventTime2 = unique(eventTime);
eventTime1 = round(diff(eventTime));
eventTime2 = round(diff(bigEvent_eventTime)*1000);

% find startTrial idx 
tmpIdx = find(ismember(bigEvent_eventName,'trialStartTrigger'));
tmpTime = diff(bigEvent_eventTime(tmpIdx))*1000;

figure; plot(timeDiff_eeg,'k.');
hold on; plot(tmpTime,'r.');
% 
% 
% tmpIdx(tmpIdx>length(eventTime2)-1)=[];
% tmpTime = eventTime2(tmpIdx);
% tmpTime_eegIdx = NaN(length(tmpTime),2);
% for i = 1:length(tmpTime)
%     [minData,minIdx] = min(abs(timeDiff_eeg-tmpTime(i)));
%     tmpTime_eegIdx(i,:) = [minData minIdx];
% end
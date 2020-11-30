mouse = 'cd019';
rootpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config';
[nFrames_oneplane, ~] = func_readnFrames([rootpath '\mouse\' mouse]);
behavpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\cd019\behavior';
sep = '\';
% behavioral variables
SESSION = 1;
TRIAL = 2;
TONE = 3;
RESPONSE = 4;
H = 1; M = 2; FA = 3; CR = 4;
NOLICK_PERIOD = 5;
RESP_TIME = 6;
DELAY_AFTER_RESP = 7;
TOTAL_TRIAL_DUR_MINUS_RESP_TIME = 8;
LICKF = 9;
REWARDF = 10;
TOTAL_TRIAL_DUR = 11; 
TONEF = 12;
CONTEXT=13;
% Default values
nPlanes = 2;
acq = 30.98/nPlanes;
pretone = 0.5; % for PSTH tone, in s
posttone = 1.5; % for PSTH tone, in s
nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];
nFrames_add = sum(nFrames_oneplane(2:end,:),2);
nFrames = diff([0;nFrames_add]);

% Get frames of 
imgData = func_loadMouseConfig(mouse,'root',rootpath);
allDay = unique(imgData.Day);
behavDay = unique(imgData.Day(strcmp(imgData.BehavType,'Behavior')));
tuningDayIdx = cellfun(@(x)any(x==behavDay), num2cell(allDay));
tuningDay = allDay(~tuningDayIdx);

nDays = length(allDay);
avg_day = cell(nDays,1);
foil_fr = [];target_fr = [];

b_count = 0;
for i=1:length(allDay)
    this_day = allDay(i);
    sessionIdx_thisday_behav = find(imgData.Day == this_day & strcmp(imgData.BehavType,'Behavior'));
    sessionIdx_thisday = find(imgData.Day==this_day);
    list_b{i} = b_count+(1:length(sessionIdx_thisday_behav)); % update the number of behavioral files today
    list_f{i} = sessionIdx_thisday;
    b_count = b_count + length(sessionIdx_thisday_behav);
end

%---------TONE FRAME FOR EACH SESSION-----------
for i = 1:length(behavDay)
    this_day = behavDay(i);
    sessionIdx_thisday_behav = find(imgData.Day == this_day & strcmp(imgData.BehavType,'Behavior'));
    sessionIdx_thisday = find(imgData.Day==this_day);
    for j = 1:length(sessionIdx_thisday_behav)
        % REVISE THIS LINE
        behavMatrix = load([behavpath sep imgData.BehavFile{sessionIdx_thisday_behav(j)}]);
        nTrials = size(behavMatrix,1);

        foil = unique([behavMatrix(behavMatrix(:,RESPONSE)==FA,TONE);behavMatrix(behavMatrix(:,RESPONSE)==CR,TONE)]);
        target = unique([behavMatrix(behavMatrix(:,RESPONSE)==H,TONE);behavMatrix(behavMatrix(:,RESPONSE)==M,TONE)]);
                
        foil_frames = behavMatrix(ismember(behavMatrix(:,TONE),foil),TONEF);
        target_frames = behavMatrix(ismember(behavMatrix(:,TONE),target),TONEF);
        nFoil = size(foil_frames,1);
        nTarget = size(target_frames,1);
        
        foilFramePlane = [foil_frames foil_frames];
        foilFramePlane(logical(mod(foil_frames,2)),:) = [round(foil_frames(logical(mod(foil_frames,2)))/nPlanes) ...
            round(foil_frames(logical(mod(foil_frames,2)))/nPlanes)-1];
        foilFramePlane(~mod(foil_frames,2),:) = [foil_frames(~mod(foil_frames,2))/nPlanes ...
            foil_frames(~mod(foil_frames,2))/nPlanes];

        targetFramePlane = [target_frames target_frames];
        targetFramePlane(logical(mod(target_frames,2)),:) = [round(target_frames(logical(mod(target_frames,2)))/nPlanes) ...
            round(target_frames(logical(mod(target_frames,2)))/nPlanes)-1];
        targetFramePlane(~mod(target_frames,2),:) = [target_frames(~mod(target_frames,2))/nPlanes ...
            target_frames(~mod(target_frames,2))/nPlanes];

        %if ~(j==1 && i==1)
        foilFramePlane = foilFramePlane+repelem(nFrames_oneplane(sessionIdx_thisday_behav(j),:),nFoil,1);
        targetFramePlane = targetFramePlane+repelem(nFrames_oneplane(sessionIdx_thisday_behav(j),:),nTarget,1);
        %end
        foil_fr = [foil_fr;foilFramePlane];
        target_fr = [target_fr;targetFramePlane];
    end
    %list_b{i} = sessionIdx_thisday_behav;
    %list_f{i} = sessionIdx_thisday;
end

toneorder = [45255 8000 13454 4757 5657,...
    22627 64000 53817 4000 9514,...
    16000 6727 19027 26909 32000,...
    11314 38055];


for i=1:length(tuningDay)
    this_day = tuningDay(i);
    sessionIdx_thisday_tuning = find(imgData.Day == this_day & (nFrames==16999) );
    sessionIdx_thisday = find(imgData.Day==this_day);
    for j = 1:length(sessionIdx_thisday_tuning)
        targIdx = find(target==toneorder);
        foilIdx = find(foil==toneorder);
        
        foil_frames = ((1700:1700:16999) + (foilIdx-1)*100)';
        target_frames = ((1700:1700:16999) + (targIdx-1)*100)';
        nFoil = size(foil_frames,1);
        nTarget = size(target_frames,1);
        
        foilFramePlane = [foil_frames foil_frames];
        foilFramePlane(logical(mod(foil_frames,2)),:) = [round(foil_frames(logical(mod(foil_frames,2)))/nPlanes) ...
            round(foil_frames(logical(mod(foil_frames,2)))/nPlanes)-1];
        foilFramePlane(~mod(foil_frames,2),:) = [foil_frames(~mod(foil_frames,2))/nPlanes ...
            foil_frames(~mod(foil_frames,2))/nPlanes];

        targetFramePlane = [target_frames target_frames];
        targetFramePlane(logical(mod(target_frames,2)),:) = [round(target_frames(logical(mod(target_frames,2)))/nPlanes) ...
            round(target_frames(logical(mod(target_frames,2)))/nPlanes)-1];
        targetFramePlane(~mod(target_frames,2),:) = [target_frames(~mod(target_frames,2))/nPlanes ...
            target_frames(~mod(target_frames,2))/nPlanes];

        %if ~(j==1 && i==1)
        foilFramePlane = foilFramePlane+repelem(nFrames_oneplane(sessionIdx_thisday_tuning(j),:),nFoil,1);
        targetFramePlane = targetFramePlane+repelem(nFrames_oneplane(sessionIdx_thisday_tuning(j),:),nTarget,1);
        %end
        foil_fr = [foil_fr;foilFramePlane];
        target_fr = [target_fr;targetFramePlane];
    
    end
    
end

start_foil = foil_fr - round(pretone*acq); % - pretone sec before tone onset
start_target = target_fr - round(pretone*acq); 
nframes_psth = round(pretone*acq) + round(posttone*acq); 

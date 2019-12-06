function go_nogo_reinforced_probe_ctrl_ZZ
%(newfile,trainingday,mousename,numtrials,session,switchtopassive,npassivetrials,numoptotuneplanes)

% Open GUI
prompt = {'Mouse name:','Training day:'};
titleWindow = 'Go/No-go task';
dims = [1 40];
if exist('D:\2Pdata\Previous_info.mat','file')
    load('D:\2Pdata\Previous_info.mat');    
    definput = {previous.name,num2str(previous.day)};
    answer = inputdlg(prompt,titleWindow,dims,definput);
else
    answer = inputdlg(prompt,titleWindow,dims);
end
mousename = answer{1};trainingday = str2double(answer{2});
path = ['D:\2Pdata\' mousename '\behavior\'];
if ~exist(path,'dir'), mkdir(['D:\2Pdata\',mousename,'\behavior\']); end
cd(path);
filelist = dir([mousename '_' num2str(trainingday) 'v*.txt']);
nFiles = size(filelist,1);
if isempty(filelist)
    session = 1;
else
    lastSession = str2double(filelist(nFiles).name(9));
    session = lastSession+1;        
end

% Default values
toneduration = 0.1; % in s
tonedelay = 0.2;
toneAmplitude = 1.5; % amplitude of tone function in TDT
solenoidOpen = 15; % number of ms to keep solenoid open
tonedurationms = toneduration*1000; % number of ms to play tone
numlicksamples = 2000000; %this is the buffer size for lick and running file, (sampling rate of 200Khz * 1/200 * 1000 seconds = 1,000,000 samples)
lickTubeOut = 0; lickTubeIn = 1;
hitdelay = 4; % ******* default should be 4, 2 for just licking
missdelay = 2;
falsealarmdelay = 7; 
nogodelay = 2;
responsewindow = 2.75;
scanningMode = 'bi';
framerate = 31.25;%30.98; changed for new scanbox, 9/13/2019
responsewindow = round(responsewindow*framerate);
nToneCounterbalanced = 10; % nb of counterbalancedtones in one bloc
nDaysFiveTarget = 5; % nb of trainingday with 5 target tones at the beginning of session
if trainingday <= nDaysFiveTarget
    nFirstTargetTones = 5;
    nToneCounterbalanced_firstbloc = 20; % nb of counterbalanced tones in the first bloc when 5 1st tones are TARGET 
    dolicktraining = 1;
else 
    nFirstTargetTones = 2;
    nToneCounterbalanced_firstbloc = 10; % nb of counterbalanced tones in the first bloc when 2 1st tones are TARGET 
    dolicktraining = 0;
end
nSwitchTargetTones = 2; % nb of target tones play after a switch (i.e. from Reinforced to Probe or vice versa)
pretone_scan = 1000; %2500; % nb of frames before starting the session (because of brightness)
npassivetrials = 20;
numoptotuneplanes = 2;
numtrials = 100;
switchtopassive = 61;

% path = ['X:\LabData1\celine\cd012\behavior'];
% Open GUI for chosing parameters
prompt = {'Session:','Nb of trials:','Switch to Passive:','Nb of Passive trials:','Lick training (0/1):','Nb of frames before beginning of session:','Nb of optotune planes:','Scanning mode (uni/di):'};
titleWindow = 'Go/No-go parameters';
dims = [1 50];
if exist(['D:\2Pdata\' mousename '\' mousename '_log.mat'],'file')
    load(['D:\2Pdata\' mousename '\' mousename '_log.mat']);    
    numtrials = prev.numtrials;
    switchtopassive = prev.switchtopassive;
    npassivetrials = prev.npassivetrials;
    dolicktraining = prev.dolicktraining;
    pretone_scan = prev.pretonescan;
    numoptotuneplanes = prev.numoptotuneplanes;
    scanningMode = prev.scanningMode;
end    
definput = {num2str(session),num2str(numtrials),num2str(switchtopassive),num2str(npassivetrials),...
    num2str(dolicktraining),num2str(pretone_scan),num2str(numoptotuneplanes),scanningMode};
answer = inputdlg(prompt,titleWindow,dims,definput);
session = str2double(answer{1});
compare = nan(nFiles,1);
for i=1:nFiles
    compare(i,1) = strcmp([mousename '_' num2str(trainingday) 'v' num2str(session) '.txt'],filelist(i).name);
end
if sum(compare)~=0
    while sum(compare)~=0
        answ = input('This session already exists. Do you want to write over it? yes/no','s');
        if strcmp(answ,'no')
            definput = {num2str(session+1),'100','',num2str(npassivetrials),num2str(dolicktraining),'',scanningMode,num2str(pretone_scan)};
            answer = inputdlg(prompt,titleWindow,dims,definput);
            session = str2double(answer{1});
            compare = nan(nFiles,1);
            for i=1:nFiles
                compare(i,1) = strcmp([mousename '_' num2str(trainingday) 'v' num2str(4) '.txt'],filelist(i).name);
            end
        elseif strcmp(answ,'yes'), break; 
        end            
    end
    session = str2double(answer{1});
end
numtrials = str2double(answer{2});
while mod(numtrials,nToneCounterbalanced)~=0
    disp(['The number of trials has to be a multiple of ' num2str(nToneCounterbalanced) '. Try again!']);
    definput = {num2str(session),'60','',num2str(npassivetrials),num2str(dolicktraining),'',scanningMode,num2str(pretone_scan)};
    answer = inputdlg(prompt,titleWindow,dims,definput);
    numtrials = str2double(answer{2});
end
switchtopassive = str2double(answer{3});
while mod(switchtopassive-1,nToneCounterbalanced)~=0
    disp(['You want the switch to passive to occur at trial n*' num2str(nToneCounterbalanced) '+1. Try again!']);
    definput = {num2str(session),'60','',num2str(npassivetrials),num2str(dolicktraining),'',scanningMode,num2str(pretone_scan)};
    answer = inputdlg(prompt,titleWindow,dims,definput);
    numtrials = str2double(answer{2});
end
npassivetrials = str2double(answer{4});
while mod(npassivetrials,nToneCounterbalanced)~=0
    disp(['The number of passive trials has to be a multiple of ' num2str(nToneCounterbalanced) '. Try again!']);
    definput = {num2str(session),'60','',num2str(npassivetrials),num2str(dolicktraining),'',scanningMode,num2str(pretone_scan)};
    answer = inputdlg(prompt,titleWindow,dims,definput);
    numtrials = str2double(answer{2});
end
switchtoactive = switchtopassive+npassivetrials;
dolicktraining = str2double(answer{5});  
pretone_scan = str2double(answer{6});
numoptotuneplanes = str2double(answer{7});   
scanningMode = answer{8};if ~strcmp(scanningMode,'bi'), framerate = 15.63; end % scanningMode = 'uni'; % 15.49;before new scanbox

% Save current session info
clear previous;
previous.name = mousename;previous.day = trainingday;
save('D:\2Pdata\Previous_info.mat','previous');
clear prev;
prev.numtrials = numtrials;
prev.switchtopassive = switchtopassive;
prev.npassivetrials = npassivetrials;
prev.dolicktraining = dolicktraining;
prev.pretonescan = pretone_scan;
prev.numoptotuneplanes = numoptotuneplanes;
prev.scanningMode = scanningMode;
save(['D:\2Pdata\' mousename '\' mousename '_log.mat'],'prev');
    
% Create behavior file
outfname = [path mousename '_' num2str(trainingday) 'v' num2str(session)];
outfname = [char(outfname),'.txt'];
dlmwrite(outfname,[]);

% Initialize and load TDT circuit, set baseline params
RP = TDTRP('C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_NLW\MatlabSyncTone_LickDetect_Frame.rcx', 'RZ6', 'INTERFACE', 'UZ3');
RP.write('BaseAmp',toneAmplitude); 
RP.write('SolenoidOn',solenoidOpen);
RP.write('ToneDur',tonedurationms); 
RP.write('numlicksamples',numlicksamples);
RP.write('ToneOn',0);
RP.write('LickPosition',lickTubeOut);

licktubeornot = questdlg('Do you want to leave the lick tube present while you search for the good plane?',...
    'Wait for User','Yes','No','Yes');
switch licktubeornot
    case 'Yes'
        RP.write('LickPosition',lickTubeIn);
        RP.write('LickTrngOn2',1);
end

moveforward = questdlg('Are you ready to start the bloc?','Wait for User','Go','Exit','Go');
RP.write('LickTrngOn2',0);
switch moveforward
  case 'Exit'
    disp('Exiting program');
  case 'Go'
    outfnamelick = [char(outfname),'licks.txt'];
    dlmwrite(outfnamelick,[]);
    outfnamewater = [char(outfname),'waterdelivery.txt'];
    dlmwrite(outfnamewater,[]);
    outfnameframe = [char(outfname),'frames.txt'];
    dlmwrite(outfnameframe,[]);
    
    % Generate tone vector for the session
    toneset = select_mouseKK(lower(mousename));%get target and foil tones
    foil = toneset(1); % foil
    target = toneset(2); % target
    initialvector = [foil*ones(nToneCounterbalanced/2,1);target*ones(nToneCounterbalanced/2,1)];
    tones = zeros(nToneCounterbalanced,150);
    for x = 1:150
        tones(:,x) = initialvector(randperm(10));
    end
    tones = reshape(tones,1500,1);
    % Now apply restrictions:
    % 1) first 5 (or 2) tones are target (i.e counterbalanced the first 20-trial (or 10-trial) bloc accordingly)
    tone_vector = [target*ones(nToneCounterbalanced_firstbloc/2-nFirstTargetTones,1);foil*ones(nToneCounterbalanced_firstbloc/2,1)];
    first_bloc_tones = [target*ones(nFirstTargetTones,1);tone_vector(randperm(nToneCounterbalanced_firstbloc-nFirstTargetTones))];
    tones(1:nToneCounterbalanced_firstbloc)= first_bloc_tones;
   
    % 2) Make sure that the first two tones after a switch are target tones
    tone_vector_swich = [target*ones(nToneCounterbalanced/2-nSwitchTargetTones,1);foil*ones(nToneCounterbalanced/2,1)];
    first_bloc_probe_tones = [target*ones(nSwitchTargetTones,1);tone_vector_swich(randperm(nToneCounterbalanced-nSwitchTargetTones))];
    tones(switchtopassive:switchtopassive+nToneCounterbalanced-1) = first_bloc_probe_tones; 
    swich_active_tones = [target*ones(nSwitchTargetTones,1);tone_vector_swich(randperm(nToneCounterbalanced-nSwitchTargetTones))];
    tones(switchtoactive:switchtoactive+nToneCounterbalanced-1) = swich_active_tones; 
   
    % Get the appropriate attenuation values from txt file
    attenuation = load('C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_tones\0.25 octave 2P behavior.txt'); 
    allfreq = attenuation(:,1);
    temp = abs(allfreq(:,1)-foil);[~,tempidx] = min(temp);
    amps_foil = attenuation(tempidx,2);
    temp = abs(allfreq(:,1)-target);[~,tempidx] = min(temp);
    amps_target = attenuation(tempidx,2);
   
    disp(['Mouse name: ' mousename ' - target = ' num2str(target) ' Hz, attenuation = ' num2str(amps_target)...
        ' dB - foil = ' num2str(foil) ' Hz,  attenuation = ' num2str(amps_foil) ' dB']);    
    amps = tones; 
    amps(tones==foil) = amps_foil;
    amps(tones==target)= amps_target;
            
    % Create summary figure
    figure('Name','Summary data','Position',[11    61   707   288]);
    active_fig=gcf;
    subplot(2,2,1);title('Reinforced Context','HandleVisibility','Off');
    xlim([0 5]);
    set(gca,'XTickLabel',{'','H' 'M' 'FA' 'CR'},'NextPlot','replacechildren');
    subplot(2,2,3);
    hold on;
    ylim([0 1]);
    xlabel('Trial Number');
    ylabel('% Correct (Hit + Correct Reject)');
    
    subplot(2,2,2);title('Probe/Passive Context','HandleVisibility','Off');
    xlim([0 5]);
    set(gca,'XTickLabel',{'' 'H' 'M' 'FA' 'CR'},'NextPlot','replacechildren');
    subplot(2,2,4);
    hold on;
    ylim([0 1]);
    xlabel('Trial Number');
    ylabel('% Correct (Hit + Correct Reject)');

    % Trigger timing and scanning for scanning
    disp('Trials will start now');
    c = RP.write('ResetFrCount',1);
    c = RP.write('ResetFrCount',1);
    c = RP.write('ResetFrCount',0);

    RP.trg(1); %trigger scanning using Soft1 in RP circuit
    
    z = RP.write('DetectLicks',1); %trigger lick detection / writing
   
    zz = RP.write('ResetData',1);
    zz = RP.write('ResetData',0);
    
    disp(['Waiting for ',num2str(pretone_scan),' frames before playing first tone...']);
    solenoidvalue=1;
    context = nan(numtrials,1);
    response = nan(numtrials,1);
    responsetime = nan(numtrials,1);
    %% ------------- trial structure, revised from original GNG task
    for trialnum = 1:numtrials
        if trialnum==1
            while RP.read('FrameNumber') ~= pretone_scan
%                 RP.write('LickPosition',lickTubeOut);
                RP.write('LickPosition',lickTubeIn);
                if dolicktraining == 1
%                     RP.write('LickPosition',lickTubeIn);
                    RP.write('LickTrngOn',1);
                end
            end
            RP.write('LickTrngOn',0);
            if trialnum==switchtopassive
                RP.write('LickPosition',lickTubeOut);
            else
                RP.write('LickPosition',lickTubeIn);
            end
        end
        
        if trialnum == switchtopassive
            RP.write('LickPosition',lickTubeOut);
            solenoidvalue = 0;
            pause(5);
        end
        
        if trialnum == switchtoactive
            RP.write('LickPosition',lickTubeIn);
            solenoidvalue = 1;
            pause(5);
        end
        lickposition = RP.read('LickPosition');
        context(trialnum) = lickposition;  
        if lickposition==1
            ind=find(context==1); plotidx=0;
        elseif lickposition==0
            ind=find(context==0); plotidx=1;
        end
        
        currentframe = RP.read('FrameNumber');
        trialStartFrame = currentframe;
        % Reset values in TDT circuit
        a = RP.write('ToneOn',0);
        d = RP.write('Solenoid',0);
        disp(['Trial ',num2str(trialnum),' starts at frame ',num2str(currentframe)]);
        tic
%         disp(['The solenoid is now set to: ',num2str(RP.read('Solenoid'))]);
%         disp(['Lickometer is now set to: ',num2str(RP.read('LickVal'))]);
        
        % Check for licking before tone onset
        lick = RP.read('LickVal');
        while lick == 1
            currentframe = RP.read('FrameNumber');
            disp(['New check frame is FIRST',num2str(currentframe),' and wait until no licking until ',num2str(currentframe+round(framerate))]);
            lick = RP.read('LickVal');
        end
        
        while RP.read('FrameNumber') < (currentframe+round(framerate)) 
            lick = RP.read('LickVal');
            if lick == 1
                currentframe = RP.read('FrameNumber');
                disp('Resetting timer...');
                disp(['New check frame is ',num2str(currentframe),' and wait until no licking until ',num2str(currentframe+round(framerate))]);
            end
        end
        
        %-----------NEW generate tone timing
        meanToneITI = 8;
        toneITI = exprnd(meanToneITI,1,10);
        toneTiming = cumsum(toneITI);
        toneFrames = round(toneTiming*framerate) + trialStartFrame; % NEED TO REVISE TO FRAME
        toneID = randi(2,1,10);
        tempToneFreq = [target, foil];
        toneFreqThisTrial = tempToneFreq(toneID);
        toneAmpThisTrial = tempToneAmp(toneID);
        toneFramePlay = [];
        toneIDPlay = [];
        toneTime = [];
        %-----------
        prelick_temp = toc;
        % Trigger the tone in TDT
%         pretone_temp = toc;
%         disp(['The tone is ',num2str(tones(trialnum)),' and attenuation is ',num2str(amps(trialnum))]);
%         RP.write('Frequency',tones(trialnum));
%         RP.write('Attenuation',amps(trialnum));
%         while mod(RP.read('FrameNumber'),numoptotuneplanes) > 0, end % ensure that tone always plays at the same optotune plane               
%         a = RP.write('ToneOn',1);
%         temp_toneon_framenumber = RP.read('FrameNumber');
%         disp(['Tone starts at frame ',num2str(temp_toneon_framenumber)]);
%         
%         pause(toneduration + tonedelay);

        % Detect lick during response period 
        lickWindowFrame = RP.read('FrameNumber');
        lick = RP.read('LickVal');
        while lick < 1 && RP.read('FrameNumber')<(temp_toneon_framenumber + responsewindow)
            lick = RP.read('LickVal');
            % ---------NEW addded detection of tone timing
            if any(toneFrames == RP.read('FrameNumber'))
               % play tone
                tempFrame = RP.read('FrameNumber');
                toneIdx = find(tempFrame == toneFrames);
                disp(['The tone is ',num2str(toneFreqThisTrial(toneIdx)),' and attenuation is ',num2str(toneAmpThisTrial(toneIdx))]);
                RP.write('Frequency',toneFreqThisTrial(toneIdx));
                RP.write('Attenuation',toneAmpThisTrial(toneIdx));
                while mod(RP.read('FrameNumber'),numoptotuneplanes) > 0, end % ensure that tone always plays at the same optotune plane               
                a = RP.write('ToneOn',1);
                temp_toneon_framenumber = RP.read('FrameNumber');
                disp(['Tone starts at frame ',num2str(temp_toneon_framenumber)]);
                tempToneTime = toc;
                
                toneTime = [toneTime tempToneTime];
                toneFramePlay = [toneFramePlay temp_toneon_framenumber];
                toneIDPlay = [toneIDPlay toneIdx];
                
                pause(toneduration + tonedelay);
            end
            %---------
        end
        
        responsetime_temp = toc - prelick_temp;
        
        if tones(trialnum) == target
            if lick == 1 %hit
                response(trialnum) = 1; %a response of 1 = Hit
                responsetime(trialnum) = responsetime_temp;
                lickframe = RP.read('FrameNumber');
                while RP.read('FrameNumber') < (lickframe + 2)
                    
                end %wait 2 frames before releasing water
                while mod(RP.read('FrameNumber'),numoptotuneplanes) > 0, end  %wait in this while loop until currentframe is divisible by the number of optotune planes
                d = RP.write('Solenoid',solenoidvalue); %Trigger and reset solenoid
                d = RP.write('Solenoid',0);
                waterframe = RP.read('FrameNumber');
                delaytime_temp = hitdelay;
                plotval='g*';
            else %miss
                response(trialnum) = 2; %a response of 2 = Miss
                responsetime(trialnum) = 999999;
                delaytime_temp = missdelay;
                plotval='r*';
                lickframe = 999999;
                waterframe = 999999;
            end
        else
            if lick == 1 %falsealarm
                response(trialnum) = 3; %a response of 3 = False Alarm
                responsetime(trialnum) = responsetime_temp;
                delaytime_temp = falsealarmdelay;
                plotval='ro';
                lickframe = RP.read('FrameNumber');
                waterframe = 999999;
            else %nogo
                response(trialnum) = 4; 
                responsetime(trialnum) = responsetime_temp;
                delaytime_temp = nogodelay;
                plotval='go';
                lickframe = 999999;
                waterframe = 999999;
            end
        end
        
        %closes the if statement for if mouse licks DURING the tone
        %end %closes the if statement for if mouse licks BEFORE the tone
        
        % Plot current status
        figure(active_fig);
        [n,x]=hist(response(ind),[1 2 3 4]);
        subplot(2,2,plotidx+1);
        bar(x,n);
        ylim([0 max(n)+1]);
        text(x,n',num2str(n','%0.0f'),... 
        'HorizontalAlignment','center',... 
        'VerticalAlignment','bottom');
        set(gca,'XTickLabel',{'H' 'M' 'FA' 'CR'});
        
        subplot(2,2,plotidx+3);
        plot(trialnum,((n(1)+n(4))/length(ind)),plotval);
        xlim([0 trialnum+1]);
        
        predelay=toc;
        %pause(delaytime_temp);
        beforeDelayFrame = RP.read('FrameNumber');
        delayFrames = beforeDelayFrame + round(delaytime_temp * framerate);
        % -------- NEW change pause into a waiting period, check tone presentation
        while RP.read('FrameNumber') < delayFrames
            % ---------NEW addded detection of tone timing
            if any(toneFrames == RP.read('FrameNumber'))
                % play tone
                tempFrame = RP.read('FrameNumber');
                toneIdx = find(tempFrame == toneFrames);
                disp(['The tone is ',num2str(toneFreqThisTrial(toneIdx)),' and attenuation is ',num2str(toneAmpThisTrial(toneIdx))]);
                RP.write('Frequency',toneFreqThisTrial(toneIdx));
                RP.write('Attenuation',toneAmpThisTrial(toneIdx));
                while mod(RP.read('FrameNumber'),numoptotuneplanes) > 0, end % ensure that tone always plays at the same optotune plane               
                a = RP.write('ToneOn',1);
                temp_toneon_framenumber = RP.read('FrameNumber');
                disp(['Tone starts at frame ',num2str(temp_toneon_framenumber)]);
                tempToneTime = toc;
                
                toneTime = [toneTime tempToneTime];
                toneFramePlay = [toneFramePlay temp_toneon_framenumber];
                toneIDPlay = [toneIDPlay toneIdx];
                
                pause(toneduration + tonedelay);
            end
            %---------
        end
        %---------
        postdelay=toc-predelay;
        disp(['The delay is ',num2str(postdelay),' second']);
        trialEndFrame = RP.read('FrameNumber');
        % Save data
        if isempty(toneIDPlay)
            data = [trainingday trialnum nan response(trialnum) nan responsetime_temp delaytime_temp (toc - responsetime_temp) lickframe waterframe toc nan context(trialnum) trialStartFrame trialEndFrame lickWindowFrame];
            dlmwrite(outfname,data,'-append');
        else
            for i = 1:length(toneIDPlay)
                data = [trainingday trialnum tempToneFreq(toneIDPlay(i)) response(trialnum) toneTime(i) responsetime_temp delaytime_temp (toc - responsetime_temp) lickframe waterframe toc toneFramePlay(i) context(trialnum) trialStartFrame trialEndFrame lickWindowFrame];
                dlmwrite(outfname,data,'-append');
            end
            
        end
    end
    %% ---------------------
    RP.trg(1);
    templickdata = RP.read('LickData',0,numlicksamples);
    tempwaterdata = RP.read('WaterDeliveryData',0,numlicksamples);
    tempframedata = RP.read('FrameData',0,numlicksamples);
    dlmwrite(outfnamelick,templickdata,' ');
    dlmwrite(outfnamewater,tempwaterdata,' ');
    dlmwrite(outfnameframe,tempframedata,' ');
    disp(['Total number of frames collected is: ',num2str(sum(diff(tempframedata)>0))]);
end



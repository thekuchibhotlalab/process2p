dataPath = 'D:\labData\excitatory\tuning\masterData\';
%mouseGNG = { 'cd017', 'cd036', 'cd037', 'cd042','cd044'};
mouseGNG = { 'cd044'};
tStay =[];tSwitch =[];tAway =[];
fStay =[];fSwitch =[];fAway =[];
midF =[];midT =[];midStay =[];
tSideStay =[];tSideT =[];tSideF =[];
fSideStay =[];fSideT =[];fSideF =[];
tUp =[];tDown =[];tStayD =[];tDownSw =[];
fUp =[];fDown =[];fStayD =[];fUpSw =[];
nUp =[];nDown =[];nStay =[];
nUpT =[];nUpF =[];nDownT =[];nDownF =[];

nRespCell = 0;
nDecodeCell = 0;

for i = 1:length(mouseGNG)
    load([dataPath mouseGNG{i} '_prePostTuning.mat']);

    tStay = [tStay prepostTuning.peak.tPre.stay];
    tSwitch = [tSwitch prepostTuning.peak.tPre.switch];
    tAway = [tAway prepostTuning.peak.tPre.away];

    fStay = [fStay prepostTuning.peak.fPre.stay];
    fSwitch = [fSwitch prepostTuning.peak.fPre.switch];
    fAway = [fAway prepostTuning.peak.fPre.away];

    midF = [midF prepostTuning.peak.middlePre.stay];
    midT = [midT prepostTuning.peak.middlePre.T];
    midStay =  [midStay prepostTuning.peak.middlePre.F];

    tSideStay = [tSideStay prepostTuning.peak.tSidePre.stay];
    tSideT = [tSideT prepostTuning.peak.tSidePre.closeT];
    tSideF = [tSideF prepostTuning.peak.tSidePre.closeF];

    fSideStay =[fSideStay prepostTuning.peak.fSidePre.stay];
    fSideT = [fSideT prepostTuning.peak.fSidePre.closeT];
    fSideF = [fSideF prepostTuning.peak.fSidePre.closeF];

    tUp = [tUp prepostTuning.decode.tChange.up];
    tDown = [tDown prepostTuning.decode.tChange.down];
    tStayD = [tStayD prepostTuning.decode.tChange.stay];
    tDownSw = [tDownSw prepostTuning.decode.tChange.downSwitch];

    fUp = [fUp prepostTuning.decode.fChange.up];
    fDown = [fDown prepostTuning.decode.fChange.down];
    fStayD = [fStayD prepostTuning.decode.fChange.stay];
    fUpSw = [fUpSw prepostTuning.decode.fChange.upSwitch];

    nUp = [nUp prepostTuning.decode.nChange.up];
    nDown = [nDown prepostTuning.decode.nChange.down];
    nStay = [nStay prepostTuning.decode.nChange.stay];

    nUpT = [nUpT prepostTuning.decode.nChange.upT];
    nUpF = [nUpF prepostTuning.decode.nChange.upF];
    nDownT = [nDownT prepostTuning.decode.nChange.downT];
    nDownF = [nDownF prepostTuning.decode.nChange.downF];
    
    nRespCell = nRespCell + prepostTuning.peak.tuningUp + ...
        prepostTuning.peak.tuningDown + prepostTuning.peak.tuningStay;
    nDecodeCell = nDecodeCell + prepostTuning.decode.nDecode;
end


figure; 
subplot(1,5,1)
tPre = [sum(tSwitch) sum(tAway) sum(tStay)];
tPre(tPre==0) = 0.1;
pie(tPre / sum(tPre),[1 1 1])
legend({'Switch2F', 'Away', 'Stay'},'Location','Best'); title ('Peak at Target')
subplot(1,5,2)
fPre = [sum(fSwitch) sum(fAway) sum(fStay)];
fPre(fPre==0) = 0.1;
pie(fPre/sum(fPre),[1 1 1])
legend({'Switch2T', 'Away', 'Stay'},'Location','Best'); title ('Peak at Foil')
subplot(1,5,3)
midPre = [sum(midT) sum(midF) sum(midStay)];
midPre(midPre==0) = 0.1;
pie(midPre / sum(midPre),[1 1 1])
legend({'move2T', 'move2F', 'Stay'},'Location','Best'); title ('Peak in Middle')
subplot(1,5,4)
tSidePre = [sum(tSideT) sum(tSideF) sum(tSideStay)];
tSidePre(tSidePre==0)= 0.1;
pie(tSidePre / sum(tSidePre),[1 1 1])
legend({'Close2T', 'Close2F', 'Stay'},'Location','Best'); title ('Peak on Tside')
subplot(1,5,5)
fSidePre = [sum(fSideT) sum(fSideF) sum(fSideStay)];
fSidePre(fSidePre==0) = 0.1;
pie(fSidePre/sum(fSidePre),[1 1 1])
legend({'Close2T', 'Close2F', 'Stay'},'Location','Best'); title ('Peak on Fside')


figure; 
subplot(1,3,1)
tDecode = [sum(tUp) sum(tDown-tDownSw) sum(tDownSw) sum(tStayD)];
tDecode(tDecode==0) = 0.1;
pie(tDecode / sum(tDecode),[1 1 1 1])
legend({'Incease', 'Decrease&No-Switch','Decrease&Switch', 'Stay'},'Location','Best'); title ('Target Preferring')
subplot(1,3,2)
fDecode = [sum(fDown) sum(fUp-fUpSw) sum(fUpSw) sum(fStayD)];
fDecode(fDecode==0) = 0.1;
pie(fDecode/sum(fDecode),[1 1 1 1])
legend({'Incease', 'Decrease&No-Switch','Decrease&Switch', 'Stay'}); title ('Foil Preferring')
subplot(1,3,3)
nDecode = [ sum(nUpT) sum(nUpF) sum(nDownT) sum(nDownF)  sum(nStay) ];
nDecode (nDecode == 0) = 0.1;
pie(nDecode / sum(nDecode),[1 1 1 1 1])
legend({'Move-To-T&End-State-T', 'Move-To-T&End-State-F','Move-To-T&End-State-F',...
    'Move-To-F&End-State-F','Stay'},'Location','Best'); title ('No Preference')

    
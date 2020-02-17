load('cd017_TC_restrict_traces_remove_neuropil.mat')

for i = 1:50
    subplot_tight(20,5,floor((i-1)/5)*10+mod(i,5)+(mod(i,5)==0)*5,[0.005 0.005])
    frame = randi(size(s_deconvolve_all_session{1},1)-2000);
    plot(s_deconvolve_all_session{i}(frame:frame+2000));
    xlim([0 2000])
    axis off
    subplot_tight(20,5,floor((i-1)/5)*10+mod(i,5)+5+(mod(i,5)==0)*5,[0.005 0.005]); 
    plot(restrict_traces{i}(frame:frame+2000));
    xlim([0 2000])
    axis off
end

%%
neuron = 16;
rawfluor = signals{1}(neuron,:,1);
dff = signals{1}(neuron,:,2);
zscore = signals{1}(neuron,:,3);

rawfluor_removed = signals_baselineremoved{1}(neuron,:,1);
dff_removed = signals_baselineremoved{1}(neuron,:,2);
zscore_removed = signals_baselineremoved{1}(neuron,:,3);

%%
range = 180000 + (36500:38000);
figure;
subplot(3,1,1)
plot(rawfluor(range))
subplot(3,1,2)
plot(dff(range))
subplot(3,1,3)
plot(zscore(range))

figure;
subplot(3,1,1)
plot(rawfluor_removed(range))
subplot(3,1,2)
plot(dff_removed(range))
subplot(3,1,3)
plot(zscore_removed(range))

%% plot neuropil extraction 
%range = (36500:38000);
%neuron = 1;
%dff{1} = signals_nothereremoved{1}(:,:);dff{2} = signals_nothereremoved{2}(:,:);
%dff{1} = signals{1}(:,:,2);dff{2} = signals{2}(:,:,2);

figure;
for i = 1:50
    frame = randi(size(dff{1},2)-2000);
    subplot_tight(20,5,floor((i-1)/5)*10+mod(i,5)+(mod(i,5)==0)*5,[0.005 0.005])
    plot(dff{1}(i,frame:frame+2000))
    xlim([0 2000])
    ylim([-1 3])
    set(gca,'xticklabels',[])
    set(gca,'yticklabels',[])
    %subplot(3,1,2)
    %plot(c_deconvolve_by_trial{1}(neuron,range))
    subplot_tight(20,5,floor((i-1)/5)*10+mod(i,5)+5+(mod(i,5)==0)*5,[0.005 0.005]); 
    plot(s_deconvolve_by_trial{1}(i,frame:frame+2000))
    xlim([0 2000])
    ylim([0 0.8])
    set(gca,'xticklabels',[])
    set(gca,'yticklabels',[])
    
end

%% plot neuropil extraction
neuron = 2;
figure;
neuron_per_row = 8;
neuron_start = 64;
for i = 1:4*neuron_per_row
    subplot_tight(16,neuron_per_row,floor((i-1)/neuron_per_row)*4*...
        neuron_per_row+mod(i,neuron_per_row)+(mod(i,neuron_per_row)==0)...
        *neuron_per_row,[0.02 0.005])
    plot(TC{1}(neuron_start+i,:));
    hold on;
    plot(neuroPil{1}(neuron_start+i,:))
    ylimm = ylim;
    ylim([ylimm(1)-4000 ylimm(2)])
    xlim([0 size(neuroPil{1},2)])
    title(['raw fluor neuron' int2str(i+neuron_start)])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    subplot_tight(16,neuron_per_row,floor((i-1)/neuron_per_row)*4*...
        neuron_per_row+mod(i,neuron_per_row)+neuron_per_row+...
        (mod(i,neuron_per_row)==0)*neuron_per_row,[0.02 0.005])
    plot(TC{1}(neuron_start+i,:) - neuroPil{1}(neuron_start+i,:))
    ylim([ylimm(1)-4000 ylimm(2)])
    title('true fluor')
    xlim([0 size(neuroPil{1},2)])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    subplot_tight(16,neuron_per_row,floor((i-1)/neuron_per_row)*4*...
        neuron_per_row+mod(i,neuron_per_row)+neuron_per_row*2+...
        (mod(i,neuron_per_row)==0)*neuron_per_row,[0.02 0.005])
    plot(signals{1}(neuron_start+i,:,2))
    ylim([-3 5])
    title('df/f, all session')
    xlim([0 size(neuroPil{1},2)])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    subplot_tight(16,neuron_per_row,floor((i-1)/neuron_per_row)*4*...
        neuron_per_row+mod(i,neuron_per_row)+neuron_per_row*3+...
        (mod(i,neuron_per_row)==0)*neuron_per_row,[0.02 0.005])
    temp = signals_nothereremoved{1}(neuron_start+i,:);
    temp(isnan(temp)) = 0;
    plot(temp)
    ylim([-3 5])
    title('df/f, session present')
    xlim([0 size(neuroPil{1},2)])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
end
function plotArtifact(TC, xoff, yoff, sessionName)
% Centralizes variables with respect to the median
nPlanes = 2;
for j = 1:nPlanes
    nTC = TC{j};
    nxoff = xoff{j};
    nyoff = yoff{j};
    
    nxoff = nxoff - median(nxoff);
    nyoff = nyoff - median(nyoff);
    
    xdiff = [diff(nxoff) 0];
    ydiff = [diff(nyoff) 0];
    
    % Finds correlation between xoff and yoff
    
    % Creates figure and first subplot of x-y shift 
    figure;
    subplot(4,2,1);plot(nxoff); hold on; plot(nyoff)
    xlabel('frames'); ylabel('shift')
    title(['x and y shift, corr. = ' num2str(corr(nxoff', nyoff'), '%0.3f')])
    legend('x shift', 'y shift', 'Location', 'Best')
    
    subplot(4,2,2);plot(xdiff); hold on; plot(ydiff)
    xlabel('frames'); ylabel('shift')
    title(['x and y shift derivatives, corr. = ' num2str(corr(xdiff', ydiff'), '%0.3f')])
    legend('x shift', 'y shift', 'Location', 'Best')

    % Creates histogram of x-y shifts and calculates mean/standard deviation
    subplot(4,2,3)
    histogram(nxoff, 'FaceAlpha', 0.5);hold on;histogram(nyoff, 'FaceAlpha', 0.5)
    ylabel('frames');xlabel('shift')
    title(['x and y shift, x avg = ' num2str(mean(nxoff), '%0.3f'),...
        ' (' num2str(std(nxoff), '%0.3f'), ')', ', y avg = '...
        num2str(mean(nyoff), '%0.3f'), ' (' num2str(std(nyoff), '%0.3f'), ')'])
    legend('x shift', 'y shift', 'Location', 'Best')


    % Computes the correlation between x and y shift with fluorescence and
    % plots histograms
    subplot(4,2,4)
    row = size(nTC,1);
    xfloCorr = zeros(1,row);
    for i = 1:row
        xfloCorr(1,i) = corr(nxoff',nTC(i,:)');
    end
    yfloCorr = zeros(1,row);
    for i = 1:row
        yfloCorr(1,i) = corr(nyoff',nTC(i,:)');
    end
    histogram(xfloCorr, 'FaceAlpha', 0.5); hold on; histogram(yfloCorr, 'FaceAlpha', 0.5);
    xlabel('shift-fluorescence corr'); ylabel('frames')
    title(['corr(dff,x) = ' num2str(mean(xfloCorr), '%0.3f'),...
        ' corr(dff,y) = ' num2str(mean(yfloCorr), '%0.3f'),])
    legend('x shift', 'y shift', 'Location', 'Best')

    % Uses x/y shifts to predict fluorescence level of each neuron and plot the
    % histogram of variance explained (r squared)
    subplot(4,2,5)
    xvar = zeros(1,row);
    for i = 1:row
        resultsx = fitlm(nxoff',nTC(i,:)');
        xvar(1,i) = resultsx.Rsquared.Ordinary;
    end
    yvar= zeros(1,row);
    for i = 1:row
        resultsy = fitlm(nyoff',nTC(i,:)');
        yvar(1,i) = resultsy.Rsquared.Ordinary;
    end
    xyvar = zeros(1,row);
    for i = 1:row
        resultsxy = fitlm([nxoff',nyoff'], nTC(i,:)');
        xyvar(1,i) = resultsxy.Rsquared.Ordinary;
    end
    histogram(xvar, 'FaceAlpha', 0.5); hold on; 
    histogram(yvar, 'FaceAlpha', 0.5); hold on; histogram(xyvar, 'FaceAlpha', 0.5);
    xlabel('r-squared value'); ylabel('frames')
    title(['variance explained between x/y shift and cell signal'])
    legend('x shift', 'y shift', 'x-y shift', 'Location', 'Best')
    
    % Uses x/y shifts derivative to predict fluorescence level of each neuron and plot the
    % histogram of variance explained (r squared)
    subplot(4,2,6)
    xvar = zeros(1,row);
    for i = 1:row
        resultsx = fitlm(xdiff',nTC(i,:)');
        xvar(1,i) = resultsx.Rsquared.Ordinary;
    end
    yvar= zeros(1,row);
    for i = 1:row
        resultsy = fitlm(ydiff',nTC(i,:)');
        yvar(1,i) = resultsy.Rsquared.Ordinary;
    end
    xyvar = zeros(1,row);
    for i = 1:row
        resultsxy = fitlm([xdiff',ydiff'], nTC(i,:)');
        xyvar(1,i) = resultsxy.Rsquared.Ordinary;
    end
    histogram(xvar, 'FaceAlpha', 0.5); hold on; 
    histogram(yvar, 'FaceAlpha', 0.5); hold on; histogram(xyvar, 'FaceAlpha', 0.5);
    xlabel('r-squared value'); ylabel('frames')
    title(['variance explained between x/y shift diff and cell signal'])
    legend('x shift', 'y shift', 'x-y shift', 'Location', 'Best')
    
    % Computes p value for correlation between x/y shift and signal per
    % neuron and plots pie chart of significant/insignificant correlations
    subplot(4,2,7)
    sig = 0;
    for i = 1:row
        f = fitlm([xdiff',ydiff'], nTC(i,:)');
        px = f.Coefficients.pValue(2);
        py = f.Coefficients.pValue(3);
        if px < 0.05 || py < 0.05
            sig = sig+1;
        end
    end
    X = [sig, (row-sig)]; labels = {'p < 0.05', 'p > 0.05'}; 
    pie(X, labels);
    title(['percent significant correlations'])
    
    % Breaks session down into 4 time segments and computes correlation
    % between cell signal and shift per segment. 
    subplot(4,2,8)
    col = size(nTC,2); 
    Q = zeros(4,row); c = 1;
    for n = 1:4
        p = (fix(col/4))*n; 
        for i = 1:row 
            f = fitlm([nxoff(c:p)',nyoff(c:p)'], nTC(i,c:p)'); 
            Q(n,i) = sqrt(f.Rsquared.Ordinary);
        end
        c = 1+p;
    end
    Q = sort(Q,2);
    maxDiff = [Q(1,1) Q(2,row) Q(3,1) Q(4,row)];
    avgCorr = [mean(Q(1,:)) mean(Q(2,:)) mean(Q(3,:)) mean(Q(4,:))];
    plot(maxDiff); hold on; plot(avgCorr);
    title(['corr(dff,shift)']); xlabel('time block'); ylabel('r value');
    legend('max diff. between corr', 'mean corr', 'Location', 'Best')
    
    % Saves figure as pdf and data as .mat file
    fileName = strcat(sessionName, num2str(j), '.png');
    saveas(gcf, fileName)
    save(fileName)
    
    % Creates new figure and computes p values from permutation test.
    % Returns histogram of p values.
    figure;
    subplot(1,2,1)
    pVal = zeros(1,row); 
    nPerm = 1000; corrValue = zeros(1,nPerm);
    for i = 1:row
        t = nTC(i,:); c = corr(ydiff',nTC(i,:)');
        for q = 1:nPerm
            permTC = t(randperm(length(t)));
            corrValue(q) = corr(ydiff',permTC');
        end
        nless = sum(corrValue < c); nequal = sum(corrValue == c);
        pVal(i) = (nless + nequal)/ length(corrValue);
    end
    nsig = 0;
    for i = 1:row
        if pVal(i) < 0.05
            nsig = nsig+1;
        end
    end
    mask = pVal < 0.05; bin_edges = 0:0.01:1;
    histogram(pVal(mask),bin_edges); hold on; histogram(pVal(~mask),...
        bin_edges,'FaceColor','red');
    xlabel('p value'); ylabel('n cells');
    title(['p values, ' num2str((nsig/length(pVal))*100), '% < 0.05']);
    legend('p < 0.05', 'p > 0.05', 'Location', 'Best')
    
    % Creates new figure and computes p values from permutation test for each neuron.
    % Returns histogram of p values with significant correlation marked.
    figure;
    subplot(2,2,2)
    pVal = zeros(1,row); nsig = 0;
    nPerm = 500; corrValue = zeros(1,nPerm); sigCorr = [];
    for i = 1:row
        t = nTC(i,:); cy = corr(ydiff',nTC(i,:)'); 
        for q = 1:nPerm
            permTC = t(randperm(length(t)));
            corrValue(q) = corr(ydiff',permTC');
        end
        nless = sum(corrValue < cy); nequal = sum(corrValue == cy);
        pVal(i) = (nless + nequal)/ length(corrValue);
        if pVal(i) < 0.025 || pVal(i) > 0.95
            nsig = nsig+1;
            sigCorr(i) = cy;
            sigCorr_YN(i) = i;
        else
            insigCorr_YN(i) = i;
        end
    end
    v = nonzeros(sigCorr'); newSigCorr = reshape(v,nsig,1)';
    diffCorr = zeros(1,row); 
    for i = 1:row
        diffCorr(1,i) = corr(ydiff',nTC(i,:)');
    end
    histogram(diffCorr, 'NumBins', 23); hold on; histogram(newSigCorr, 'FaceColor','red', 'NumBins', 23);
    xlabel('corr'); ylabel('frames');
    title(['corr(dff,ydiff), ' num2str((nsig/length(pVal))*100), '% < 0.05']);
    legend('p > 0.05', 'p < 0.05', 'Location', 'Best')
    
    % Computes correlations using permutation test with ABS VALUE of y diff
    subplot(2,2,4)
    pVal = zeros(1,row); nsig = 0; y_AB = abs(ydiff);
    nPerm = 500; corrValue = zeros(1,nPerm); sigCorr = [];
    for i = 1:row
        t = nTC(i,:); cy = corr(y_AB',nTC(i,:)'); 
        for q = 1:nPerm
            permTC = t(randperm(length(t)));
            corrValue(q) = corr(y_AB',permTC');
        end
        nless = sum(corrValue < cy); nequal = sum(corrValue == cy);
        pVal(i) = (nless + nequal)/ length(corrValue);
        if pVal(i) < 0.025 || pVal(i) > 0.95
            nsig = nsig+1;
            sigCorr(i) = cy;
            sigCorr_YAB(i) = i;
        else
            insigCorr_YAB(i) = i;
        end
    end
    v = nonzeros(sigCorr'); newSigCorr = reshape(v,nsig,1)';
    diffCorr = zeros(1,row); 
    for i = 1:row
        diffCorr(1,i) = corr(y_AB',nTC(i,:)');
    end
    histogram(diffCorr, 'NumBins', 23); hold on; histogram(newSigCorr, 'FaceColor','red', 'NumBins', 23);
    xlabel('corr'); ylabel('frames');
    title(['corr(dff,|ydiff|), ' num2str((nsig/length(pVal))*100), '% < 0.05']);
    legend('p > 0.05', 'p < 0.05', 'Location', 'Best')
    
    % Computes correlations using permutation test with x diff
    subplot(2,2,1)
    pVal = zeros(1,row); nsig = 0;
    nPerm = 500; corrValue = zeros(1,nPerm); sigCorr = [];
    for i = 1:row
        t = nTC(i,:); cx = corr(xdiff',nTC(i,:)'); 
        for q = 1:nPerm
            permTC = t(randperm(length(t)));
            corrValue(q) = corr(xdiff',permTC');
        end
        nless = sum(corrValue < cx); nequal = sum(corrValue == cx);
        pVal(i) = (nless + nequal)/ length(corrValue);
        if pVal(i) < 0.025 || pVal(i) > 0.95
            nsig = nsig+1;
            sigCorr(i) = cy;
            sigCorr_XN(i) = i;
        else
            insigCorr_XN(i) = i;
        end
    end
    v = nonzeros(sigCorr'); newSigCorr = reshape(v,nsig,1)';
    diffCorr = zeros(1,row); 
    for i = 1:row
        diffCorr(1,i) = corr(xdiff',nTC(i,:)');
    end
    histogram(diffCorr,'NumBins',23); hold on; histogram(newSigCorr, 'FaceColor','red', 'NumBins', 23);
    xlabel('corr'); ylabel('frames');
    title(['corr(dff,xdiff), ' num2str((nsig/length(pVal))*100), '% < 0.05']);
    legend('p > 0.05', 'p < 0.05', 'Location', 'Best')
    
    % Computes correlations using permutation test with ABS VALUE of x diff
    subplot(2,2,3)
    pVal = zeros(1,row); nsig = 0; x_AB = abs(xdiff);
    nPerm = 500; corrValue = zeros(1,nPerm); sigCorr = [];
    for i = 1:row
        t = nTC(i,:); cx = corr(x_AB',nTC(i,:)'); 
        for q = 1:nPerm
            permTC = t(randperm(length(t)));
            corrValue(q) = corr(x_AB',permTC');
        end
        nless = sum(corrValue < cx); nequal = sum(corrValue == cx);
        pVal(i) = (nless + nequal)/ length(corrValue);
        if pVal(i) < 0.025 || pVal(i) > 0.95
            nsig = nsig+1;
            sigCorr(i) = cy;
            sigCorr_XAB(i) = i;
        else
            insigCorr_XAB(i) = i;
        end
    end
    v = nonzeros(sigCorr'); newSigCorr = reshape(v,nsig,1)';
    diffCorr = zeros(1,row); 
    for i = 1:row
        diffCorr(1,i) = corr(x_AB',nTC(i,:)');
    end
    histogram(diffCorr,'NumBins',23); hold on; histogram(newSigCorr, 'FaceColor','red', 'NumBins', 23);
    xlabel('corr'); ylabel('frames');
    title(['corr(dff,|xdiff|), ' num2str((nsig/length(pVal))*100), '% < 0.05']);
    legend('p > 0.05', 'p < 0.05', 'Location', 'Best')
    
    % Computes cross correlation for each neuron with x/y diff. 
    % Averages correlogram for x/y diff with all neurons and plot avg. line.
    figure;
    subplot(1,2,1)
    for i = 1:row
        [cx, lagsx] = xcorr(xdiff',nTC(i,:)',20, 'normalize'); 
        [cy, lagsy] = xcorr(ydiff',nTC(i,:)',20,'normalize');
        for q = 1:41
            CCx(q,i) = cx(q); CCy(q,i) = cy(q);
        end
    end
    for i = 1:41
        meanCCx(i) = mean(CCx(i,:)); meanCCy(i) = mean(CCy(i,:));
        SEM_x(i) = std(meanCCx)/sqrt(length(meanCCx)); SEM_y(i) = std(meanCCy)/sqrt(length(meanCCy)); 
    end
    x1 = (meanCCx + SEM_x); x2 = (meanCCx - SEM_x); x = -20:20; 
    fill([x fliplr(x)], [x1 fliplr(x2)], 'b', 'FaceAlpha', 0.3, 'LineStyle', 'none'); hold on;
    plot(x, meanCCx, 'b'); 
    xlabel('lag'); ylabel('xcorr');
    title(['Xcorr between x diff and TC, avg. lag = ' num2str(mean(meanCCx),'%0.4f')]);
    subplot(1,2,2)
    y1 = (meanCCy + SEM_y); y2 = (meanCCy - SEM_y); 
    fill([x fliplr(x)], [y1 fliplr(y2)], 'r', 'FaceAlpha', 0.3, 'LineStyle', 'none'); hold on;
    plot(x, meanCCy, 'r'); 
    xlabel('lag'); ylabel('xcorr');
    title(['Xcorr between y diff and TC, avg. lag = ' num2str(mean(meanCCy),'%0.5f')]);
    
    % Creates figure with sample cross correlograms of individual neurons.
    figure;
    p = (fix(row/4));
    for i = 1:2:7
        subplot(4,2,i);
        [c, lags] = xcorr(xdiff',nTC(p,:),20,'normalize');
        plot(lags, c);
        xlabel('lag'); ylabel('xcorr');
        title(['Xcorr between xdiff and signal, cell ' num2str(p)]);
        p = p + 50;
    end
    p = (fix(row/4));
    for i = 2:2:8
        subplot(4,2,i);
        [c, lags] = xcorr(ydiff',nTC(p,:)',20,'normalize');
        plot(lags, c);
        xlabel('lag'); ylabel('xcorr');
        title(['Xcorr between ydiff and signal, cell ' num2str(p)]);
        p = p + 50;
    end
    
    % Separates individual calcium traces by significance of correlation
    % with x/y shift and absolute value of x/y shift. Returns average
    % calcium trace for each group.
    sigCorr_YN = nonzeros(sigCorr_YN)'; sigCorr_XN = nonzeros(sigCorr_XN)';
    sigCorr_YAB = nonzeros(sigCorr_YAB)'; sigCorr_XAB = nonzeros(sigCorr_XAB)';
    insigCorr_YN = nonzeros(insigCorr_YN)'; insigCorr_XN = nonzeros(insigCorr_XN)';
    insigCorr_YAB = nonzeros(insigCorr_YN)'; insigCorr_XAB = nonzeros(insigCorr_XN)';
    for i = 1:length(sigCorr_XN); n = sigCorr_XN(i); A = nTC(n,:); sigTC_XN(i,:) = A; end;
    for i = 1:length(sigCorr_YN); n = sigCorr_YN(i); A = nTC(n,:); sigTC_YN(i,:) = A; end;
    for i = 1:length(sigCorr_XAB); n = sigCorr_XAB(i); A = nTC(n,:); sigTC_XAB(i,:) = A; end;
    for i = 1:length(sigCorr_YAB); n = sigCorr_YAB(i); A = nTC(n,:); sigTC_YAB(i,:) = A; end;
    for i = 1:length(insigCorr_XN); n = insigCorr_XN(i); A = nTC(n,:); insigTC_XN(i,:) = A; end;
    for i = 1:length(insigCorr_YN); n = insigCorr_YN(i); A = nTC(n,:); insigTC_YN(i,:) = A; end;
    for i = 1:length(insigCorr_XAB); n = insigCorr_XAB(i); A = nTC(n,:); insigTC_XAB(i,:) = A; end;
    for i = 1:length(insigCorr_YAB); n = insigCorr_YAB(i); A = nTC(n,:); insigTC_YAB(i,:) = A; end;
    figure; subplot(2,2,1); plot(mean(sigTC_XN)); hold on; plot(mean(insigTC_XN)); 
    legend('sig','insig'); title(['x corr']);
    subplot(2,2,2); plot(mean(sigTC_YN)); hold on; plot(mean(insigTC_YN)); 
    legend('sig','insig'); title(['y corr']);
    subplot(2,2,3); plot(mean(sigTC_XAB)); hold on; plot(mean(insigTC_XAB)); 
    legend('sig','insig'); title(['x ABS corr']);
    subplot(2,2,4); plot(mean(sigTC_YAB)); hold on; plot(mean(insigTC_YAB)); 
    legend('sig','insig'); title(['y ABS corr']);
    
end


end
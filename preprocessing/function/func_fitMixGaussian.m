function param = func_fitMixGaussian (targDistribution)
    [targpdf,tempxbins] = histcounts(targDistribution);
    for i = 1:length(tempxbins)-1
       xbins(i) = mean(tempxbins(i:i+1));
    end
    targpdf = targpdf ./ sum(targpdf);
    param0 = [mean(targDistribution)-300, mean(targDistribution)+300,...
        std(targDistribution), std(targDistribution)];
    param = fmincon(@(x)func_getMixGaussianCost...
        (xbins,targpdf, x),param0, [],[]);
end


function mixgauss = func_getMixGaussian(xbins, param)
    mu1 = param(1);
    mu2 = param(2);
    std1 = param(3);
    std2 = param(4);
    gauss1 = normpdf(xbins, mu1, std1);
    gauss2 = normpdf(xbins, mu2, std2);
    mixgauss = gauss1 + gauss2;
    mixgauss = mixgauss ./ sum(mixgauss);
end

function cost = func_getMixGaussianCost(xbins,targpdf,param)
    mixgauss = func_getMixGaussian(xbins, param);
    cost = sum((mixgauss - targpdf).^2);
end
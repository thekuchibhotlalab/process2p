function filenames = fn_fileListAutoSelect(configTable,selectedOpt)

cateAllFiles = {'extractTC','extractRedraw','deconvolve','roiTracking','roiRevision','roiRevisionPlot','saveMeanImg'};
cateSelectFiles = {'getTuning','none'};

if contains(selectedOpt,cateAllFiles); filenames = configTable.('ImagingFile'); 
elseif contains(selectedOpt,cateSelectFiles); filenames = 'SelectFile'; 
else; filenames = 'SelectFile'; 
end
    

end
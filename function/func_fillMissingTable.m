function T = func_fillMissingTable(T)

missingCol = T.Properties.VariableNames(all(ismissing(T),1));
nRow = size(T,1);
tempCol = cell(nRow,1);
tempCol(:) = {' '};

for i = 1:length(missingCol)
    T.(missingCol{i}) = tempCol;
end

end
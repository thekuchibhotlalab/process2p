function newdir = getDirNAbove(mydir,n)
tempStr   = strsplit(mydir,'\');
newdir = strjoin(tempStr(1:end-n),'\');
end
function sbx2tifKK(fname,varargin)

% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is converted

z = sbxread(fname,1,1);
global info;

if(nargin>1)
    N=varargin{1};%N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end

k = 0;
done = 0;
while(~done && k<=N)
    try
        q = sbxread(fname,k,1);
        green = squeeze(q(1,:,:));
        red = squeeze(q(2,:,:));
        if(k==1)
            imwrite(green,[fname '_green.tif'],'tif');
            imwrite(red,[fname '_red.tif'],'tif');
        else
            imwrite(green,[fname '_green.tif'],'tif','writemode','append');
            imwrite(red,[fname '_red.tif'],'tif','writemode','append');
        end
    catch
        done = 1;
    end
    k = k+1;
end
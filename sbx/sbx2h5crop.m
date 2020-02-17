function sbx2h5crop(fname,varargin)

% sbx2h5
% Generates h5 file from sbx files
%cropping paramaters for optotune and bidi

x_crop_start = 110; %x-value to start strips
x_crop_end = 0;
y_crop_start = 100; %y-value to start strips
y_crop_end = 0;

fnh = [fname ,'.h5'];

z = sbxread(fname,1,1);
global info;
 
if(nargin>1)
N = min(varargin{1},info.max_idx);
else
N = info.max_idx;
end
 
k = 0;
done = 0;
 
blksize = 200; % block size
 
to_read = min(blksize,N-k);
 
while(~done && to_read>0)
    try
        q = sbxread(fname,k,to_read);
        q = squeeze(q(1,:,:,:)); % extract green channel only
        q = permute(q,[2 1 3]);
        q = q(y_crop_start:(end-y_crop_end),x_crop_start:(end-x_crop_end),:); %crop optotune and bidi
        if(k==0)
            h5create(fnh,'/data',[697 403 Inf],'DataType','uint16','ChunkSize',[697 403 to_read]);
            h5write(fnh,'/data',q,[1 1 1],[697 403 to_read]);
%             f = waitbar(0,'Converting to hdf5');
        else
            h5write(fnh,'/data',q,[1 1 k+1],[697 403 to_read]);
        end
    catch
        done = 1;
%         delete(f);
    end
    k = k+to_read;
    to_read = min(blksize,N-k);
% 	waitbar(k/N,f);
end
 
% delete(f);
end


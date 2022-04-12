function sbx2h5cropdir(fname,targetDir,channel,varargin)

% sbx2h5
% Generates h5 file from sbx files

if nargin <3; channel = 1; end
x_crop_start = 110; %x-value to start strips
x_crop_end = 0;
y_crop_start = 100; %y-value to start strips
y_crop_end = 0; % take y from 100:755
 
fnh = [targetDir '\' fname ,'.h5'];

z = sbxread(fname,1,1);
global info;
 
if(nargin>3)
N = min(varargin{1},info.max_idx);
else
N = info.max_idx;
end
 
k = 0;
done = 0;
 
blksize = 2000; % block size
 
to_read = min(blksize,N-k);
 
while(~done && to_read>0)
try
q = sbxread(fname,k,to_read);
q = squeeze(q(channel,:,:,:)); % extract green channel only
q = permute(q,[2 1 3]);
q = q(y_crop_start:(end-y_crop_end),x_crop_start:(end-x_crop_end),:); %crop optotune and bidi
if(k==0)
h5create(fnh,'/data',[size(q,1) size(q,2) Inf],'DataType','uint16','ChunkSize',[size(q,1) size(q,2) to_read]);
h5write(fnh,'/data',q,[1 1 1],[size(q,1) size(q,2) to_read]);
f = waitbar(0,'Converting to hdf5');
else
h5write(fnh,'/data',q,[1 1 k+1],[size(q,1) size(q,2) to_read]);
end
catch
done = 1;
delete(f);
end
k = k+to_read;
to_read = min(blksize,N-k); 
waitbar(k/N,f);
end
 
delete(f);
end


function func_attachRefImg(h5path, h5name,refImg,varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('rep', 20000)
p.addParameter('root', pwd)

p.parse(varargin{:});
sep = '\';

frames = h5read([h5path sep h5name],'/data');

repImg = repmat(refImg,1,1,20000);

end
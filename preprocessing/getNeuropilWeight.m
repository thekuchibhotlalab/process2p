function weight = getNeuropilWeight(imgSize, rois, varargin)
    warning off
    p = inputParser;
    p.addParameter('minRadii', 0); % radius around the cells to exclude
    p.addParameter('maxRadii', 20);
    p.addParameter('avgMethod', 'mean');
    p.addParameter('avgParams', []);
    p.addParameter('coefficient', 0.7);
    p.parse(varargin{:});

    imgSizeX = imgSize(1);
    imgSizeY = imgSize(2);
    %currentRoiPos = allRoi{currentRoi};
    
    neuroPilMask = false(imgSizeX, imgSizeY);
    for i = 1:length(rois)
        x = rois{1,i}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
        y = rois{1,i}.mnCoordinates(:,2);
        

        bw = roipoly(zeros(imgSizeX, imgSizeY),y,x); % y=col; x=rows
        [xt,yt]= find(bw>0);
        
        mask = createCirclesMask([imgSizeX, imgSizeY], [y x], ones(size(x))* p.Results.minRadii);
        neuroPilMask = neuroPilMask | mask | bw;
        %subplot(121)
        %plot(roiPolyShape)
        %subplot(122)
        %imagesc(neuroPilMask)
        %close all
    end

    xaxis = 1:imgSizeX;
    yaxis = 1:imgSizeY;

    [xgrid, ygrid] = meshgrid(xaxis, yaxis);
    
    weight = [];
    for i = 1:length(rois)
        
        x = rois{1,i}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
        y = rois{1,i}.mnCoordinates(:,2);
        roiPolyShape = polyshape(x,y);
        [roiCenterX, roiCenterY] = centroid(roiPolyShape);
        
        
        switch p.Results.avgMethod
            
            case 'gaussian'  
                xGaussProb = normpdf(xaxis, roiCenterX, p.Results.avgParams(1));
                yGaussProb = normpdf(yaxis, roiCenterY, p.Results.avgParams(1));
                [xGaussGrid, yGaussGrid] = meshgrid(yGaussProb, xGaussProb);
                gaussGrid = xGaussGrid .* yGaussGrid;
                gaussGrid(neuroPilMask==1) = 0; 
                gaussGrid = gaussGrid ./ sum(gaussGrid(:));
                weight{i} = gaussGrid;
            case 'mean'
                mask = createCirclesMask([imgSizeX, imgSizeY], [roiCenterY roiCenterX], p.Results.maxRadii);
                mask(neuroPilMask==1) = 0;
                mask = mask ./ sum(mask(:));
                weight{i} = mask;
        end
        weight{i} = weight{i} * p.Results.coefficient;
       
    end
    warning on
    
end

function mask = createCirclesMask(varargin)
%xDim,yDim,centers,radii)
% Create a binary mask from circle centers and radii
%
% SYNTAX:
% mask = createCirclesMask([xDim,yDim],centers,radii);
% OR
% mask = createCirclesMask(I,centers,radii);
%
% INPUTS: 
% [XDIM, YDIM]   A 1x2 vector indicating the size of the desired
%                mask, as returned by [xDim,yDim,~] = size(img);
%  
% I              As an alternate to specifying the size of the mask 
%                (as above), you may specify an input image, I,  from which
%                size metrics are to be determined.
% 
% CENTERS        An m x 2 vector of [x, y] coordinates of circle centers
%
% RADII          An m x 1 vector of circle radii
%
% OUTPUTS:
% MASK           A logical mask of size [xDim,yDim], true where the circles
%                are indicated, false elsewhere.
%
%%% EXAMPLE 1:
%   img = imread('coins.png');
%   [centers,radii] = imfindcircles(img,[20 30],...
%      'Sensitivity',0.8500,...
%      'EdgeThreshold',0.30,...
%      'Method','PhaseCode',...
%      'ObjectPolarity','Bright');
%   figure
%   subplot(1,2,1);
%   imshow(img)
%   mask = createCirclesMask(img,centers,radii);
%   subplot(1,2,2);
%   imshow(mask)
%
%%% EXAMPLE 2:
%   % Note: Mask creation is the same as in Example 1, but the image is
%   % passed in, rather than the size of the image.
%
%   img = imread('coins.png');
%   [centers,radii] = imfindcircles(img,[20 30],...
%      'Sensitivity',0.8500,...
%      'EdgeThreshold',0.30,...
%      'Method','PhaseCode',...
%      'ObjectPolarity','Bright');
%   mask = createCirclesMask(size(img),centers,radii);
%
% See Also: imfindcircles, viscircles, CircleFinder
%
% Brett Shoelson, PhD
% 9/22/2014
% Comments, suggestions welcome: brett.shoelson@mathworks.com

% Copyright 2014 The MathWorks, Inc.

narginchk(3,3)
if numel(varargin{1}) == 2
	% SIZE specified
	xDim = varargin{1}(1);
	yDim = varargin{1}(2);
else
	% IMAGE specified
	[xDim,yDim] = size(varargin{1});
end
centers = varargin{2};
radii = varargin{3};
xc = centers(:,1);
yc = centers(:,2);
[xx,yy] = meshgrid(1:yDim,1:xDim);
xx = double(xx);
yy = double(yy);
mask = false(xDim,yDim);
for ii = 1:numel(radii)
	mask = mask | hypot(xx - xc(ii), yy - yc(ii)) <= radii(ii);
end
end

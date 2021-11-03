c = cell2mat(sessionMeanImg);
c = reshape(c,403,64,697);
c = permute(c,[1 3 2]);
implay(int16(c))
meanImg = mean(c(:,:,1:60),3);
imwrite(uint16(meanImg),'meanImg_plane1.tiff','tiff')
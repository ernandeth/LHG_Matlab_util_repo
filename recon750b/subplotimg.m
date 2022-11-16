function subplotimg(img,plotLabels,strTitle)

imgMin = min(img(:));
imgMax = max(img(:));
img = img - imgMin;
img = img./(imgMax - imgMin).*64;

numImgs = size(img,3);
numInRow = ceil(sqrt(numImgs));
numWholeRows = floor(numImgs/numInRow);
numExtra = numImgs - (numWholeRows*numInRow);
if numExtra > 0
    numRows = numWholeRows + 1;
else
    numRows = numWholeRows;
end

H = figure;
set(H,'pos',[337 93 876 702]);
if exist('strTitle','var')
    set(H,'name',strTitle);
end

for ind = 1:numImgs
    h = subplot(numRows,numInRow,ind);
    image(img(:,:,ind));
    set(h,'xticklabel','');
    set(h,'yticklabel','');
    xlabel(num2str(plotLabels(ind)));
    colormap(gray);
end
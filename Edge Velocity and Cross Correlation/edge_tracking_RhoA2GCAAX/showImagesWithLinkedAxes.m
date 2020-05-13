function showImagesWithLinkedAxes(imArr, style, lim, titles)
%Displays a set of images with linked axes

if nargin<2
    style='imagesc';
end
if nargin<4
    titles=[];
end
if (iscell(imArr) && length(imArr)==1 && size(imArr{1},3)>1) || (isnumeric(imArr) && size(imArr,3)>1)
    if iscell(imArr)
        temp=imArr{1};
    else
        temp=imArr;
    end
    clear imArr;
    for i=1:size(temp,3)
        imArr{i}=temp(:,:,i);
    end
end
scrsz=get(0,'ScreenSize');
figure('Position',[1 50 scrsz(3)-2 scrsz(4)-150]);
num=length(imArr);
for i=1:num
    ax(i)=subplot(1,num,i);
    if strcmpi(style,'imagesc')
        if nargin>=3
            imagesc(imArr{i},lim);
        else
            imagesc(imArr{i});
        end
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
    else
        if nargin<3
            lim=[];
        end
        imshow(imArr{i},lim);
    end
    if iscell(titles) & ~isempty(titles)
        title(titles{i});
    end
end
linkaxes(ax);

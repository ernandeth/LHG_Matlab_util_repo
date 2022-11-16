close all

a = zeros(100,100);
[x,y] = meshgrid(-50:50,-50:50);


a = exp(-((x-20).^2)/120 - ((y-20).^2)/80) + ...
    exp(-((x+20).^2)/120 - ((y-20).^2)/80) + ...
    exp(-((x-20).^2)/120 - ((y+20).^2)/80) + ...
    exp(-((x+18).^2)/120 - ((y+18).^2)/80);

%a = 1./sqrt((x-30).^2+ (y-20).^2 + eps);

mask = exp(-((x).^2)/180 - ((y).^2)/150);
imagesc(a)
colormap(gray)

mask(mask>0.03)=1;
mask(mask<=0.03)=0;

%[dx, dy] = gradient(mask);
%edge = circshift(dx,[0,1]) + circshift(dy,[1,0]);
%edge = abs(dx) + abs(dy);
%edge = find(edge ~= 0);
%[ex ey ]= ind2sub( size(dx) , edge);
%edge = [ex' ; ey'];

% Use this to read actual data
load sensitivities_0503_1
mask = bodymask;
a = squeeze(abs(sens64(4,:,:)));
%%%%


a = a.*mask;
maskContour = contour(mask,1);
ex = maskContour(1,2:end);
ey = maskContour(2,2:end);

el = fit_ellipse(ex,ey);
edge = make_ellipse(el.X0, el.Y0, el.a, el.b, el.phi);

figure
imagesc(a)
colormap gray
hold on
plot(edge(1,:), edge(2,:))
axis tight

figure
[c, h] = contour(a,20);
guess0 = [1 ,1,1,1];
optvar = optimset('lsqnonlin');
p=1;

c2 =[]
allellipses = [];

while p< size(c,2)
    % figure out how many points for the current contour:
    N = c(2,p);  % number of pairs
    L = c(1,p);  % intensity level

    % extract the coordinates of each contour
    tmp = c(:, p+1:p+N);
    plot(tmp(1,:), tmp(2,:));

    hold on

    % throw out half of each contour
    %incmplt = tmp(:,1:end/2);
    %plot(incmplt(1,:), incmplt(2,:),'r');

    %tx = incmplt(1,:);
    %ty = incmplt(2,:);

    % thrown out those points in the contour that are on the edge of the
    % brain
    tmp2 = [];
    fprintf('\rremoving edges from contour..%d . ', p)
    for cpoint=1:size(tmp,2)
        
        % check all points on the edge relative to this point
        removePt=0;
        for epoint = 1:size(edge,2)
            if ( sqrt(sum((tmp(:,cpoint) - edge(:,epoint)).^2)) <= 1)
                %fprintf ('found an edge point... ')
                removePt = 1;
                epoint = size(edge,2);
            end
        end
        if removePt==0
            tmp2 = [tmp2 tmp(:,cpoint)];
        end
    end
    if ~isempty(tmp2)
        tx = tmp2(1,:);
        ty = tmp2(2,:);
        % fit the incomplete ellipse:
        plot(tx,ty,'r')

        fprintf('...fitting ellipse')

        el = fit_ellipse(tx,ty);
        allellipses = [allellipses el];
        
       pause
    end
    % advance pointer to the next contour
    p = p + N +1;

end

N = size(allellipses,2);
for e = 1:N
    if  ~isempty(el)
        if isempty(el.status)
            my_ellipse = make_ellipse(el.X0, el.Y0, el.a, el.b, el.phi);
            plot(my_ellipse(1,:), my_ellipse(2,:), 'g');
            hold on
            drawnow


            % stick the results into a new contour (in matlab format)
            newContour = zeros(3,length(my_ellipse));
            %    newCountour(1,1) = L;
            %    newContour(2,1) = length(my_ellipse);
            %    newContour(1:2,2:end) = my_ellipse;
            newContour(1:2,:) = my_ellipse;
            newContour(3,:) = L;
            c2= [c2 newContour];
        end
    end
end
[x,y] = meshgrid(1:size(mask,1),1:size(mask,2));
result = griddata(c2(1,:), c2(2,:), c2(3,:), x,y);
figure
imagesc(result)
colormap gray
hold on
plot(edge(1,:), edge(2,:),'.')

function Intensity =  array_beam_path_sim;
%%
Npts = 50;
Nelements = 24;
attenuation = 0.015;

plane = 'sagittal'

switch (plane)
    case 'axial'
        % define a barrier (skull) using the tanh function
        bx = linspace(0,10,Npts)';
        by = 4*tanh(0.8*bx) ;
        by(end/2+1:end) = by(end/2:-1:1);
        bx = bx -5;
        
        B = [bx by];
        
        subplot(221)
        cla
        plot(bx,by)
        hold on
        
        % make an array of  sources parallel to the barrier with Nelements
        sx = linspace(0,10,Nelements)';
        sy = 4.5*tanh(0.8*sx) ;
        sy(end/2+1:end) = sy(end/2:-1:1);
        sx = 1.5*(sx -5);
        
        source = [sx sy];
        source(end/2-6:end/2+7,:) = [];
        %source = source(4:end-3,:);
        
        plot(source(:,1), source(:,2),'ro')
        axis square
        
    case 'sagittal'
        % define a barrier (skull) using the tanh function
        bx = linspace(0,10,Npts)';
        by = 4*bx.^0.3 ;
        bx = bx;
        B = [bx by];
        
        subplot(221)
        cla
        plot(bx,by)
        hold on
        
        % make an array of  sources parallel to the barrier with Nelements
        sy = linspace(0,max(by),Nelements)';
        sx = ((sy) .^ 3);
        sx = (sx)*8/max(sx);
        sx = sx - 2;
        source = [sx sy];
        %source(end/2-6:end/2+7,:) = [];
        %source = source(4:end-3,:);
        source = source(1:2:end,:) ;
        
        plot(source(:,1), source(:,2),'ro')
        axis square
end
drawnow
% define a set of targets under the barrier:
R =  linspace(min(bx), max(bx),Npts);
targets = [];
%%
Intensity = zeros(Npts);

% compute beam for each source:
for q = 1:size(source,1)
    S = source(q,:);
    
    % effect of beam on each target of the grid
    % Npts x Npts
    for m= 1:Npts
        for n = 1:Npts
            
            % if target is below the barrier:
            if R(n) <= by(m)
                T =  [R(m), R(n)];
                %plot(R(m), R(n),'k.')
                
                targets = [
                    targets;
                    R(m), R(n);
                    ];
                
                
                % line equation from target to source
                beam_path = get_line_equation(T,S);
                
                
                % Figure out where the beam hits the barrier by testing all
                % the elements in the barrier (B)
                min_dist_source=inf;
                x_seg = [];
                for p=1:size(B,1)-1
                    % line for each segment of the barrier
                    seg = get_line_equation( B(p,:) , B(p+1,:));
                    
                    % location where the beam_path interesects with the
                    % line describing each segment.
                    tmp_crossing = find_intersection(beam_path, seg);
                    %plot(tmp_crossing(1), tmp_crossing(2),'x')
                    
                    % length of each segment:
                    segLength = norm(B(p,:) - B(p+1,:));
                    % center of the segment
                    Bcenter =  (B(p,:) + B(p+1,:))/2;
                    % plot(Bcenter(1), Bcenter(2),'o')
                    % distance from segment to where lines interesect
                    seg2crossing = norm(tmp_crossing - Bcenter);
                    % distance from segment to source
                    seg2source = norm(Bcenter - S);
                    
                    % if the crossing is near the middle of the segment,
                    % and this segment is the closest one to the source
                    % then this is the right place: x_point
                    if (seg2crossing < segLength) && (seg2source < min_dist_source)
                        
                        %plot(tmp_crossing(1), tmp_crossing(2),'kx')
                        
                        min_dist_source = seg2source;
                        x_point = tmp_crossing;
                        x_seg = seg; % This is the segment where it crosses
                    end
                end
                
                if ~isempty(x_seg)
                    % calculate  angle between the beam_path and the segment?
                    x_theta = calc_intersection_angle(x_seg, beam_path);
                    
                    target2barrier = norm(x_point - T);
                    % accumulate the intensity from each beam.  They are
                    % weighted by the angle of the intersection with the
                    % barrier
                    Intensity(m,n) = Intensity(m,n) + sin(x_theta)*exp(-attenuation*target2barrier);
                    %plot(x_point(1), x_point(2),'rx')
                end
            end
        end
    end
end
axis square

subplot(222)
imagesc(Intensity')
axis square
axis xy
axis image

PP = Intensity(Intensity>0);
subplot(223)
hist(PP(:),50);

Threshold = max(PP(:)) - std(PP(:))
PP = Intensity;
PP(PP<Threshold) = 0;

subplot(224)
imagesc(PP');
axis square
axis xy
axis image

return

function theta = calc_intersection_angle(line1, line2)
% angle of incidence
% make two vectors using the slopes for the beam_path and the crossing segment
v1 =[1 line1(1)*1 ];
v2 =[1 line2(1)*1 ];
dotprod = v1*v2';
theta = (acos( dotprod / norm(v1)/norm(v2)));
return

function result = get_line_equation( s,t)
% return the slope and intercept of a line through two points
% s=[sx, sy] and t = [tx ty]
%
m = (s(2)-t(2)) / (s(1) - t(1) +eps);
b = s(2) - m*s(1);

result = [m, b];
return

function crossing = find_intersection(line1, line2)
% find the find_intersectionion between two lines mx + b ---> [m b]
%
m1 = line1(1); b1 = line1(2);

m2 = line2(1); b2 = line2(2);

x = (b2-b1) / (m1 - m2 );
y = m2*x + b2;

crossing = [x y];

return



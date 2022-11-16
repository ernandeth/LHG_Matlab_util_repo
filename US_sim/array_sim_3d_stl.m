function P = array_sim_3d_stl(fname)
Npts = 50;
attenuation = 0.15;

B = [];
if nargin==0
    fname = 'average305_t1_tal_lin.stl';
end
[f v n] = stlread(fname);
v = v/10; % convert to cm
% look only at the top
tmp = v(:,3);
v(find(tmp < 5.0),:)=[];
% shift down so origin is in nasion
v(:,3) = v(:,3)-5;

% shift x so origin is in nasion
v(:,1) = v(:,1)- max(v(:,1))/2;


B = v(1:500:end,:);

% define a barrier (skull) using the tanh function
W = max(v(:,1)) - min(v(:,1));
D = max(v(:,2)) - min(v(:,2));
H = max(v(:,3)) - min(v(:,3));

% Make a mesh inside the bone
xcoords = 1.2*linspace(-W,W,Npts)/2;
ycoords = 1.2*linspace(-D,2,Npts);
zcoords = 1.2*linspace(0,H, Npts);

[xmesh ymesh zmesh] = meshgrid(xcoords, ycoords, zcoords);


cla
plot3(B(:,1), B(:,2), B(:,3), '.')
axis([-10 10 -10 10 -10 10])
xlabel('x');
ylabel('y')
zlabel('z')
hold on



%% Define the array of transducer elements
Nelementsx = 16;
Nelementsz = 10;
S = [];
W2 = W+4;
H2 = H+2;
D2 = D+2;
sclz = 1;
crv = 0.4;

% model as a polynomial based on cadaver data
%load surf_data.mat
% model as a polynomial based on MNI data
load surf_data_MNI.mat

y_offset = max(axial(:,2));
axial(:,2) = axial(:,2)-y_offset;

y_offset = abs(max(sagittal(:,1)));
sagittal(:,1) = sagittal(:,1)-y_offset;


p_axial = polyfit(axial(:,1), axial(:,2), 6);
p_sagittal = polyfit(sagittal(:,1), sagittal(:,2), 6);

% The rings get smaller as we move up:
sz_array = linspace(0,H2, Nelementsz);  % slice location (z axis)

% slice shrinkage factor along y
slfactory = abs(polyval(p_sagittal, sz_array));
slfactory = slfactory/max(slfactory);
slfactory = 1- slfactory;

% slice shrinkage factor along x
slfactorx =  (H2-sz_array).^0.2 / (H2^0.2) ;

sx = linspace(-W2,W2,Nelementsx)/2;
sy = (D2/2)*tanh(crv*(W2/2 - abs(sx))) - W/2 ;
sy_offset = min(sy);

for z = 1:length(sz_array)
    % ring gets narrower as we go up:    
    tmpx = sx * slfactorx(z) ;
    tmpy = (sy - sy_offset) * slfactory(z) + sy_offset ;
    tmpz = ones(size(sx)) * sz_array(z);
    
    S = [S;
            tmpx' tmpy' tmpz'];
    
end

%single source test:
% S = [9,9,2];

plot3(S(:,1), S(:,2), S(:,3), 'r*')
axis(1.2*[-W2 W2 -W2 W2 -W W2])
xlabel('x');
ylabel('y')
zlabel('z')
hold off
drawnow

% compute the contributions from each source to the target field
% loop over the entire volume
P = zeros(Npts, Npts, Npts);
bone_x_coords = B(:,1);
bone_y_coords = B(:,2);
bone_z_coords = B(:,3);
mm=1;
deltay = D/Npts/2;

parfor z = 1:Npts
    for x = 1:Npts
        for y = 1:Npts
            target = [x y z];
            target_coords = [xcoords(x) ycoords(y) zcoords(z)];
            
            % only need the pressure at points behind the skull
            % test to see if the points are inside:
            %
            isinside = 0;
            ncrossing = 0;
            testpoint = target_coords;
            for n=1:Npts*4
                testpoint(2) = testpoint(2) + deltay;
                for m=1:size(B,1)
                    tmp = norm(testpoint - B(m,:));
                    if (tmp < deltay/2)
                        ncrossing=ncrossing+1;
                        %{
                        hold on
                        plot3(testpoint(1), testpoint(2), testpoint(3),'ko')
                        drawnow
                        plot3(B(m,1), B(m,2), B(m,3),'ro')
                        %}
                    end
                end
            end
            if (ncrossing==1)
                isinside = 1;
            end
            %}
            if (1)
            %if (isinside)
                %{
                hold on
                plot3(target_coords(1), target_coords(2), target_coords(3),'g*')
                drawnow
                hold on
                %}
                isinside=0;
                for n=1:size(S,1)  % loop over the sources
                    % find the line between source and target
                    source = S(n,:);
                    
                    % get line equation point-vector form
                    L = get_line_equation3d(source,target_coords);
                    
                    % find the nearest points on the skull nearest to the line
                    % check only the fron half of the skull
                    tmp_dists = inf(size(B,1),1);
                    for m=1:size(B,1)/2
                        x0 = B(m,:);
                        tmp_dists(m) = distance_point2line(x0, target_coords, source);
                        % weigh the distances by the distance to the source
                        tmp_dists(m) = tmp_dists(m)*norm(x0-source);
                    end
                    [tmp_dists inds] = sort(tmp_dists, 'ascend');
                    r=B(inds(1),:); % this is the closest point to the "entry"
                    s=B(inds(10),:);
                    t=B(inds(20),:);
                    entry_point = r;
                    tissue_path_length = norm(entry_point-target_coords);
                        
                    arr = B(inds(1:Npts/2),:);
                    % plane = get_plane_equation(r,s,t);
                    plane = get_plane_equation(arr);
                    incident_angle = calc_line_plane_angle(L, plane);
                                          
                    P(x,y,z) = P(x,y,z) + abs(sin(incident_angle))* exp(-attenuation*tissue_path_length);

                    %{
                    mm = mm+1;
                    if mm==500
                        hold off
                        plot3(B(1:3:end,1), B(1:3:end,2), B(1:3:end,3), '.')
                        hold on
                        plot3(arr(:,1), arr(:,2), arr(:,3),'rx');
                        plot3(target_coords(1), target_coords(2), target_coords(3),'g*')
                        plot3(source(1), source(2), source(3),'r*')
                        ray = [target_coords; source];
                        plot3(ray(:,1), ray(:,2), ray(:,3))
                        
                        axis(1.2*[-W2 W2 -W2 W2 -W2 W2])
                        xlabel('x');
                        ylabel('y')
                        zlabel('z')
                        drawnow
                        %junk(:)=0;
                        mm=1;
                        
                    end
                    %}
                    
                end
                %}
            end
        end
    end
    
end
return
end




function theta = calc_line_plane_angle(L, plane)
% get the angle between a line (between points x1 and x2)
% and a plane specified by: [a b c d]
% in the plane equation: ax + by + cz = d;
% the line comes in point-vector form
% [x0,  a
%  y0   b
%  z0   c]

% first find the normal to the plane from the equation
Nplane = plane(1:3);
Nline = L(:,2);

theta = acos( (Nplane*Nline) / norm(Nline) / norm(Nplane));

return
end

function theta = calc_intersection_angle(line1, line2)
% angle of incidence
% make two vectors using the slopes for the beam_path and the crossing segment
v1 =[1 line1(1)*1 ];
v2 =[1 line2(1)*1 ];
dotprod = v1*v2';
theta = (acos( dotprod / norm(v1)/norm(v2)));
return
end

function d = distance_point2line(x0, x1, x2);
% computes the distance between x0 and the line L
% L is specified as two points on the line (x1 and x2)
%
x12 = x1-x2;
x10 = x0-x1;
x2 = x2-x1;

d = norm(cross(x10, x12)) / norm(x12);

return
end

function d = distance_point2plane(p, S);
% distance from point p to a plane S
% specified by the equation
%       [a b c d]*[x y z 1]' = 0
% point p = [x0, y0, z0]
% plane S = [a b c d]

pp = [p 1];
d = abs(pp * S(:)) / norm(S(1:3));

return
end


function result = get_line_equation( s,t)
% return the slope and intercept of a line through two points
% s=[sx, sy] and t = [tx ty]
%
m = (s(2)-t(2)) / (s(1) - t(1) +eps);
b = s(2) - m*s(1);

result = [m, b];
return
end

function coeffs = get_line_equation3d(s,t);
% generate the equation for a line defined by two points in in R3
% points are s,t
% returns the coefficients in the equation
% x = x0 + nvec*t

x0 = s;
nvec = t-s;
nvec = nvec/norm(nvec);
coeffs = [x0(:)   nvec(:)];

return
end


function coeffs = get_plane_equation(arr);
%function coeffs = get_plane_equation(s,t,r);
% generate the equation for a plane defined by severa points in that plane
% points are s,t,r
% returns the coefficients in
% a*x(n) + b*y(n) + c*z(n) +d = 0
% -a*x(n) =  b*y(n) + c*z(n) + d
% in matrix form:
% X = A*coeffs

n = max(size(arr));
X = arr(:,3);
A = [arr(:, 2:3)  ones(n,1)];

coeffs = pinv(A)*X;
coeffs = [1; coeffs]';
return
end

function crossing = find_intersection(line1, line2)
% find the find_intersectionion between two lines mx + b ---> [m b]
%
m1 = line1(1); b1 = line1(2);

m2 = line2(1); b2 = line2(2);

x = (b2-b1) / (m1 - m2 );
y = m2*x + b2;

crossing = [x y];

return
end
%%
function gen_template_data
% from the MNI template:
axial = [-78 0;
    -76 12;
    -74 24;
    -80 -4 ;
    -72 36;
    -66 40;
    -68 34;
    -48 70;
    -60 56;
    -44 72;
    -34 76;
    -28 80;
    -20 82;
    -12 84;
    -2 84;
    8 84;
    16 84;
    16 84;
    20 80;
    30 80;
    36 76;
    42 74;
    48 70;
    56 62;
    66 42;
    70 32;
    74 24;
    76 14;
    80 4;
    ]/10;

sagittal = [
    86 -16;
    86 0;
    86 14;
    82 22;
    82 32;
    80 34;
    76 44;
    70 50;
    66 60;
    58 66;
    50 70;
    44 76;
    36 80
    30 82]/10;

sagittal(:,2) = sagittal(:,2) + 1.6;
sagittal(:,1) = sagittal(:,1) - min(sagittal(:,1));

save surf_data_MNI.mat axial sagittal

return


end


% Create Regressors from input stimulus and HRF
% Input Stimulus:  delta functions centered in the middle of the reaction time
% HRF Gamma Variate function
%
% the file times.mat must be present.  It contains:
%  	allcontrolends         90x1            720  double array
%  allcontrolstarts       90x1            720  double array
%  allends               144x4           4608  double array
%  allstarts             144x4           4608  double array
%  endcontroltimes        15x1            120  double array
%  lengths               666x1           5328  double array
%  onsets                666x1           5328  double array
%  startcontroltimes      15x1            120  double array



clear resp
load times

% Constants for GAmma variate function
delta=2500/2000;
tau = 1250/2000;

allstarts=[allstarts zeros(size(allstarts,1)  , 1 ) ];
allends=[allends zeros(size(allends,1)  , 1 ) ];
allstarts(1:size(allcontrolstarts,1), 5) = allcontrolstarts;
allends(1:size(allcontrolends,1), 5) = allcontrolends;

allmids = (allstarts + allends)/2;

% Create the input stimulus
input=zeros(1200, 5);
for col=1:5
  for row=1:size(allmids,1)
     input(floor(allmids(row,col))+1 , col ) = 1;
  end
end
 
plot(input)

%pause

% Create the HRF from a gamma variate function
gam=zeros(200,1);
dummy=gam;

for t=1:size(gam,1)
   gam(t)=((t-delta)./tau).^2.*exp(-(t-delta)./tau).*((t-delta) > 0);
end
plot(gam);
%pause


% Convolve the HRF with the input
for col=1:5
   resp(:,col)=conv(input(:,col),gam);
end
resp=resp(1:1200,:);

noswitch = resp(:,1);
counter = resp(:,2);
operation=resp(:,3);
dual = resp(:,4);
baseline=resp(:,5);

allmidpoints = allmids;


save regressors2 resp allmidpoints baseline operation counter operation dual

colormap('gray');
imagesc(resp);

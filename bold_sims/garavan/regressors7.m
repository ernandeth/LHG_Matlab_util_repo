% Create Regressors from input stimulus and HRF
% Input Stimulus:  square wave starting at 180 msec and ending 
% at 550 msec. after the onset of the stimulus 
% (NB: the lags between trials have been removed)
%
% HRF: Gamma Variate function
%
% the file 'times.mat' must be present.  It contains:
%  allcontrolends         90x1            720  double array
%  allcontrolstarts       90x1            720  double array
%  allends               144x4           4608  double array
%  allstarts             144x4           4608  double array
%  endcontroltimes        15x1            120  double array
%  lengths               666x1           5328  double array
%  onsets                666x1           5328  double array
%  startcontroltimes      15x1            120  double array
% 
%  The file 'lagstart.mat' must also be present.  It contains the lagtimes


clear resp
load times
load lagstart


% Constants for GAmma variate function
delta=2500/2000;
tau = 1250/2000;

allstarts=[allstarts zeros(size(allstarts,1)  , 1 ) ];
allends=[allends zeros(size(allends,1)  , 1 ) ];
allstarts(1:size(allcontrolstarts,1), 5) = allcontrolstarts;
allends(1:size(allcontrolends,1), 5) = allcontrolends;

allstarts = allstarts + 180/2000;
allends = allstarts + 370/2000;


% Shift the data to account for the removal of the scans
% where nothing was happening (lag between trials)
cleanstarts = allstarts;
cleanends = allends;
shift = 0;
for run = 2:7
	shift = shift + 200 * (run -1 ) - lagstart(run -1)
	start_index = find( allstarts > (run-1)*200  &  allstarts < run*200 );
	end_index = find( allends > (run-1)*200  &  allends < run*200 );
   cleanstarts(start_index) = allstarts(start_index) - shift;
   cleanends(end_index) = allends(end_index) - shift;
end

allstarts = cleanstarts;
allends = cleanends;

% Create the input stimulus
input=zeros(1200, 5);
for col=1:5
  for row=1:size(allstarts,1)
     input(floor(allstarts(row,col))+1 : floor(allends(row,col))+1, col ) = 1;
  end
end
colormap('gray') 
imagesc(input)

pause


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
resp=resp(1:1200 - shift,:);

noswitch = resp(:,1);
counter = resp(:,2);
operation=resp(:,3);
dual = resp(:,4);
baseline=resp(:,5);

save regressors7 resp  baseline operation counter operation dual noswitch cleanstarts cleanends

colormap('gray');
imagesc(resp);

% arrowtimes.m: 
% 
% Extract the start and end times of the arrow displays in the 
% garavan study.  Compute the durations, midppoints of the events, and write 
% a set of .mat files with the times for every particular kind of event.
% the variables generated can then be used for SPM99.

for filenum=1:6  %

   % create the filename string
   filename = strcat(basename,num2str(filenum),'-r.txt')
   block = read_mat(filename,14);

   % extract the event times for the regular trial types
   for i=3:6
      endcol = i*2;  
      startcol = i*2 - 1 ;
      starttimes(:,i-2) = block(find(block(:,startcol) > 0.1),startcol);
      endtimes(:,i-2) = block(find(block(:,endcol) > 0.1),endcol);
   end
   
   midtimes = (endtimes+starttimes)/2;
   durations = endtimes-starttimes;
	
   % extract the event times for the baseline data (there are less control events)
   startcontroltimes = block(find(block(:,13) > 0.1),13);
   endcontroltimes =  block(find(block(:,14) > 0.1),14);
   
   midcontroltimes = (endcontroltimes+startcontroltimes)/2;
   controldurations = endcontroltimes-startcontroltimes;
   
   % Figure out the time in scans units (2 sec./scan)
   starttimes = starttimes/2000 + 200*(filenum-1);
   startcontroltimes = startcontroltimes/2000 + 200*(filenum-1);
   endtimes = endtimes/2000 + 200*(filenum-1);
   endcontroltimes = endcontroltimes/2000 + 200*(filenum-1);
   midtimes = midtimes/2000 + 200*(filenum-1);
   midcontroltimes = midcontroltimes/2000 + 200*(filenum-1);
   durations = durations/2000;
   controldurations = controldurations/2000;
   
   
 

	if filenum >= 2
      alltrials=[alltrials; midtimes];
      allcontrols=[allcontrols; midcontroltimes];
      
      allstarts=[allstarts; starttimes];
      allcontrolstarts = [allcontrolstarts; startcontroltimes];
      allends=[allends; endtimes];
      allcontrolends=[allcontrolends; endcontroltimes];
      
      alldurations =[alldurations; durations];
      allcontroldurations = [allcontroldurations; controldurations];
	
	else
   	alltrials=midtimes;
     allcontrols = midcontroltimes;
     
     allstarts=starttimes;
     allcontrolstarts = startcontroltimes;
     allends=endtimes;
     allcontrolends= endcontroltimes;

     alldurations = durations;
     allcontroldurations = controldurations;
      
  end 
  
  
end

% Arrange start times and durations into a single column.
onsets=0;
lengths=0;
for col=1:4
    onsets = [onsets; allstarts(:,col)];
   lengths = [lengths ; alldurations(:,col)];
  end
onsets=onsets(2:size(onsets,1));
lengths=lengths(2:size(lengths,1));
  
onsets = [onsets; allcontrolstarts];
lengths =[lengths; allcontroldurations];

%create variables to be saved

noswitch=alltrials(:,1);
counter=alltrials(:,2);
operation=alltrials(:,3);
dual=alltrials(:,4);
baseline = allcontrols;


save times onsets lengths allstarts allends allcontrolstarts allcontrolends startcontroltimes endcontroltimes
save midpoints noswitch counter operation dual baseline
save alltimes allstarts alldurations allcontrolstarts allcontroldurations

write_mat(alltrials, 'alltrialtimes.txt');
write_mat(allcontrols,'allcontroltimes.txt');


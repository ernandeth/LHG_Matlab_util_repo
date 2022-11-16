
for filenum=1:6

   % create the filename string
   filename = strcat('dyblock',num2str(filenum),'-r.txt')
   block = read_mat(filename,14);

   
   for i=3:6
      endcol = i*2;  
      startcol = i*2 - 1 ;
      starttimes(:,i-2) = block(find(block(:,startcol) > 0.1),startcol);
      endtimes(:,i-2) = block(find(block(:,endcol) > 0.1),endcol);
      midtimes = (endtimes+starttimes)/2;
	end
   
   startcontroltimes = block(find(block(:,13) > 0.1),13);
   endcontroltimes =  block(find(block(:,14) > 0.1),14);
   midcontroltimes = (endcontroltimes+startcontroltimes)/2;
   
   % Figure out the time in terms of scans 
   midtimes = midtimes/2000 + 200*(filenum-1);
   midcontroltimes = midcontroltimes/2000 + 200*(filenum-1);
 

	if filenum >= 2
      alltrials=[alltrials; midtimes];
      allcontrols=[allcontrols; midcontroltimes];
	
	else
   	alltrials=midtimes;
   	allcontrols = midcontroltimes;
	end 


end

%create variables to be saved

noswitch=alltrials(:,1);
counter=alltrials(:,2);
operation=alltrials(:,3);
dual=alltrials(:,4);
baseline = allcontrols;

save midpoints noswitch counter operation dual baseline

write_mat(alltrials, 'alltrialtimes.txt');
write_mat(allcontrols,'allcontroltimes.txt');

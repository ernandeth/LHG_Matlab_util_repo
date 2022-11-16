function tmscall(action)

%     tmscall(action)
% Callback Switchyard for tmswin.m application
% Makes all the calls for the buttons and options of tmswin.m
% Luis Hernandez
% Last edit: 1-30-98

% Create an array of 5 tms structures with a file name and the data 
% proper.

global active tms

         
switch action
   
  case 'create'
     %variable initialization 
     active = 0;
     

  case 'load_xyz'
     % Create array with position data set from a file containing polhemus 
     % coordinates.  (x,y,z,pitch,roll,yaw)
     
      if active < 5
         [n p]= uigetfile('', 'Select Matrix Data File');
         cd (p);
         tms.name = strcat(p,n);
         tms.data = read_mat(tms.name,6);
     disp(tms.data) 
         hold on;  
         tms.colorstr = getcolor(active+1);  
         plot3( tms.data(:,1) ,...
            tms.data(:,2),...
            tms.data(:,3),...
            tms.colorstr);
         
         % Next: Update the active file slider
         active = active + 1;
         add_files(tms.name);
         labels;
     end
      

   case 'load_tms'
      % Create array with position and intensity data set from a file
      if active < 5
         [n p] = uigetfile('', 'Select Matrix Data File');  
         cd (p);
         tms.name = n;
         tms.data = read_mat(tms.name,4);
         
         hold on;
         tms.colorstr = getcolor(active +1);
         plot3( tms.data(:,1) ,...
            tms.data(:,2),...
            tms.data(:,3),...
            tms.colorstr);
         % Next: Update the active file slider
         active = active + 1;
         add_files(tms.name);
      end
      labels;

   case 'save_file'
      % Stores data as a four column matrix.
      [n p] = uiputfile('*.dat','Enter Output File name');
      name = strcat(p,n);
      write_mat( tms.data, name);
      
      

   case 'exit'
      % Exit the Program: clear global variables and close window
      main_window = findobj('Tag', 'Fig1')
      close (main_window)
      clear active Tmswin tms

      

   case 'ChangeActiveFile'
      % Change the active file number and name to those in the file
      % listbox.  
      % Reload the file into the tms data structure
      %lbox = findobj('Tag','Listbox1');
      filestring  = get(gcbo, 'String');
      index = get(gcbo,'Value');      
      
      % Not Working!!!!!!!!!!!!   
      

   case 'ClearFiles'
      % Remove the current plot and clear data variables
      clear tms
      active = 0;
      cla;
      set(gca,'Color',[0 0 0]);
      labels;
      
      lbox = findobj('Tag','Listbox1');
      set(lbox, 'String','');
      resultbox = findobj('Tag','ResultBox');
      set(resultbox,'String', '' );


   case 'remove_junk'
      % Keep Every Nth point in the data set.  Throw all others out.
             
      % clear plots and redraw the original set
      hold off;
      plot3( tms.data(:,1) ,...
         tms.data(:,2),...
         tms.data(:,3),...
         tms.colorstr);
      
      % Remove points that are out of range
      tms.data = remove(tms.data, -30, 30);
      
      % Draw the transformed data in the same color, but with different markers
      hold on ;
      tms.colorstr(2) = 'o';
      
      plot3( tms.data(:,1) ,...
         tms.data(:,2),...
         tms.data(:,3),...
         tms.colorstr);
      
      [n p] = uiputfile('*.dat','Enter Output File name');
      name = strcat(p,n);
      write_mat(tms.data, name);

   case 'fourths'
      
      % keep only every 4th point:  !
       tms.data = reduce(tms.data,4)

	% Draw the extracted data
      hold on ;
      
      plot3( tms.data(:,1) ,...
         tms.data(:,2),...
         tms.data(:,3),...
         'mo');
      
      [n p] = uiputfile('*.dat','Enter Output File name');
      name = strcat(p,n);
      write_mat(tms.data, name);

      
   case 'correct_motion'
      % Motion correction (requires a six column position data set)
      % Assumes that this is from a Polhemus data set (6 columns)
      
      % Redraw original data set
      hold off;
      plot3( tms.data(:,1) ,...
         tms.data(:,2),...
         tms.data(:,3),...
         tms.colorstr);
    
      tms.data = translation(tms.data);
      tms.data = rotation(tms.data);
      
      % remove the rotation angles and add a column of zeros for intensities
      tms.data = tms.data(:,1:3);
      sz = size(tms.data);
      tms.data = [tms.data(:,1:3) zeros(sz(1),1)];
      
      hold on;
      tms.colorstr(2) = 'o';
      plot3( tms.data(:,1) ,...
         tms.data(:,2),...
         tms.data(:,3),...
         tms.colorstr);
      labels;
      
      [n p] = uiputfile('*.dat','Enter Output File name');
      name = strcat(p,n)

      write_mat(tms.data, name);
   
   case 'fit'
      % Fit data from one coordinate frame to another.  Frames
      % are defined by two sets of fiducials

      transformed = tmsfit(tms.data);
      [n p] = uiputfile('*.dat','Enter Output File name');
      name = strcat(p,n);
      write_mat(transformed, name);

      plot3( transformed(:,1) ,...
         transformed(:,2),...
         transformed(:,3),...
         'wo');
      labels;

      
   case 'get_intensities'
      % Open a file containing a single columns with TMS intensities
      % and copy its contents into the fourth column of the data file
      
      [n p]= uigetfile('', 'Select Intensities File');
      cd (p);
         name = strcat(p,n);
         intensity_data = read_mat(name,1);
         szi = size(intensity_data);
         szp = size(tms.data);
         if szi(1) ~=szp(1)
            errormesg('Failure. Dimensions must agree !');
         else
            tms.data(:,4) = intensity_data(:);
            [n p] = uiputfile('*.dat','Enter Output File name');
            name = strcat(p,n);
            write_mat(name, tms.data);
         end
         

   case 'make_analyze'
   % Convert active tms data into an Analyze forma image
	
       tms2analyze(tms.data);

   case 'get_area'
      % Compute the area defined by a three dimensional data set
      area = area3(tms.data);
      resultbox = findobj('Tag','ResultBox');
      set(resultbox,'String', num2str(area) );
      
      
   case 'get_center'
      % Compute weighted average position of the data.
      % The weight of each point is determined by its 
      % corresponding intensity
      ctrpoint = stim_ctr(tms.data)
      resultbox = findobj('Tag','ResultBox');
      set(resultbox,'String', num2str(ctrpoint) );
      

   case 'project'
      % Calculate the projection of a data set to
      % a surface.  Requires an Analyze format volume.

	tms.data = tmsproject(tms.data);
	  
	plot3( tms.data(:,1) ,...
            tms.data(:,2),...
            tms.data(:,3),...
            tms.colorstr);
         
     tms.colorstr = getcolor(active+1);   
	write_mat(tms.data, strcat('p-',tms.name));

      
   end
   
return


function labels
h = gca;

xlabel('X');
ylabel('Y');
zlabel('Z');

return



function colorstring = getcolor(active)
% This function decides the plot's color and type of symbol
% to be plotted.

  switch(active)
   
  case 1
     colorstring = 'g*';
  case 2
     colorstring = 'b*';
  case 3
     colorstring = 'c*';
  case 4
     colorstring = 'w*';
  case 5
     colorstring = 'y*';
  end
  
return
   


function add_files(name)
% This function updates the active file list box
% by listing the contents of the tms structure
    global active 

    lbox = findobj('Tag','Listbox1');
    filestring  = get(lbox, 'String');
   
    if ~isempty(filestring)
       filestring = char(filestring, name);
    else
       filestring = name;
    end
    
    set(lbox, 'String', filestring);
    set(lbox, 'Value',active);
    
return













function dmalign2(action);
%
% 	This file contains a comprehensive form of distance based alignment algorithm
%	please read the accompanying text
%
% 	Cengizhan
% 	11.1.1996
%
if nargin<1,
	action= 'initialize';
end;
%
if strcmp(action,'initialize'),
	clf;
	fig = gcf;
	
	h1=text(.1,.9,'ALIGN (v1.1)');
	axis off;
	set(h1,'FontSize', 24);
	h2=text(.1,.8,'by  Cengizhan Ozturk, MD, PhD');
	set(h2,'FontSize', 18);
	
	h3=text(.1,.4,'Imaging and Computer Vision Center');
 	h4=text(.1,.33,'Drexel University, Philadelphia');
 	
 	% h5=text(.1,.45,'Director             :   Oleh J. Tretiak, Ph.D.');
	% set(h5,'Color',[1 0 0],'FontSize', 10);
	% h6=text(.1,.38,'Asst. Director   :   Jonathan Nissanov, Ph.D.');
	% set(h6,'Color',[1 0 0],'FontSize', 10);
	
	h7=text(.1,.60,'Medical Imaging Laboratory');
 	h8=text(.1,.53,'Johns Hopkins University School of Medicine');
 	
  	%h10=text(.1,.38,' Director   :  Elliot R. McVeigh, Ph.D.');
	%set(h10,'Color',[1 0 0],'FontSize', 9);
	
	pb1 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .05 .2 .1],...
		'String','ABOUT',...
		'Units','normalized',...
		'CallBack','dmalign2(''about'')');
	pb2 =uicontrol(fig,...
		'Style','push',...
		'Position',[.7 .05 .2 .1],...
		'String','ALIGN',...
		'Units','normalized',...
		'CallBack','dmalign2(''start'')');
	pb3 =uicontrol(fig,...
		'Style','push',...
		'Position',[.45 .05 .2 .1],...
		'Units','normalized',...
		'String','DEMOS',...
		'CallBack','dmalign2(''demos'')');
	pb3x =uicontrol(fig,...
		'Style','push',...
		'Position',[.45 .17 .2 .1],...
		'Units','normalized',...
		'String','GRAYSCALE',...
		'CallBack','dmalign2(''grayscale'')');
	pb4 =uicontrol(fig,...
		'Style','push',...
		'Position',[.7 .17 .2 .1],...
		'Units','normalized',...
		'String','QUIT',...
		'CallBack','close; clear; clear global;');
		

elseif strcmp(action,'about'),
	clf;
	fig = gcf;
 
	h1=text(.1,.9,'ALIGN (v1.1)');
	axis off;
   	h2=text(.1,.8,  '  for distance based alignment ');
    	h4=text(.1,.6, ' Transformations: Rigid and Affine ');
 	h5=text(.1,.5, ' Search: Iterative, modified gradient ');
 	h6=text(.1,.4, '              and euclidean distance map based');
 	h3=text(.1,.3, ' Contact: Dr. Cengizhan Ozturk (cozturk@mri.jhu.edu)');
 	h2=text(.1,.2, '  11.1.1996  (last mod : 7.6.1997) ');
 
 	pb3 =uicontrol(fig,...
		'Style','push',...
		'Position',[.7 .1 .2 .1],...
		'String','RETURN',...
		'Units','normalized',...
		'CallBack','dmalign2(''initialize'')');

elseif strcmp(action,'demos'),
	clf;
	fig = gcf;
	txt_1 = uicontrol(fig,...
		'BackgroundColor',[0 0 0],...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left',...
		'Style','text',...
		'String','Ellipse, Rigid Body',...
		'Position',[.4 .7 .6 .05]',...
		'Units','normalized');
		
	pb_d1 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .7 .15 .05],...
		'Units','normalized',...
		'String','Demo 1',...
		'CallBack','dmalign2(''demo1'')');

	txt_2 = uicontrol(fig,...
		'BackgroundColor',[0 0 0],...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left',...
		'Style','text',...
		'String','Cookie, Affine',...
		'Position',[.4 .6 .6 .05],...
		'Units','normalized');
	pb_d2 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .6 .15 .05],...
		'Units','normalized',...
		'String','Demo 2',...
		'CallBack','dmalign2(''demo2'')');
	
	txt_3 = uicontrol(fig,...
		'BackgroundColor',[0 0 0],...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left',...
		'Style','text',...
		'String','2D Cookie with occlusion, Rigid Body',...
		'Position',[.4 .5 .6 .05],...
		'Units','normalized');
	pb_d3 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .5 .15 .05],...
		'Units','normalized',...
		'String','Demo 3',...
		'CallBack','dmalign2(''demo3'')');
	
	txt_4 = uicontrol(fig,...
		'BackgroundColor',[0 0 0],...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left',...
		'Style','text',...
		'String','2D Occluded Cookie with outlier removal',...
		'Position',[.4 .4 .6 .05],...
		'Units','normalized');
	pb_d4 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .4 .15 .05],...
		'Units','normalized',...
		'String','Demo 4',...
		'CallBack','dmalign2(''demo4'')');
	
	txt_5 = uicontrol(fig,...
		'BackgroundColor',[0 0 0],...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left',...
		'Style','text',...
		'String','3D ellipsoid, Rigid Body',...
		'Position',[.4 .3 .6 .05],...
		'Units','normalized');
	pb_d5 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .3 .15 .05],...
		'Units','normalized',...
		'String','Demo 5',...
		'CallBack','dmalign2(''demo5'')');
		
	pb3 =uicontrol(fig,...
		'Style','push',...
		'Position',[.7 .1 .2 .1],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','dmalign2(''initialize'')');

elseif strcmp(action,'demo1'),
	clf;
	disp('DEMO 1');
	m2D(9001,9001,300,300,0,20,1,1,1,0,0.000001,0.000001,200,2,2,0,20);
	pb3 =uicontrol(gcf,...
		'Style','push',...
		'Position',[.8 0 .2 .07],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','close all;dmalign2(''demos'')');
		
elseif strcmp(action,'demo2'),
	clf;
	disp('DEMO 2');
	m2Daff(9002,9002,100,100,0,40,1,1,1,0,0.0000001,200,2,2,0,20);
	pb3 =uicontrol(gcf,...
		'Style','push',...
		'Position',[.8 0 .2 .07],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','close all;dmalign2(''demos'')');
		
elseif strcmp(action,'demo3'),
	clf;
	disp('DEMO 3');
	m2D(9003,9003,100,100,0,20,1,1,1,0,0.000001,0.000001,200,2,2,0,20);
	pb3 =uicontrol(gcf,...
		'Style','push',...
		'Position',[.8 .05 .1 .05],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','close all;dmalign2(''demos'')');
		
elseif strcmp(action,'demo4'),
	clf;
	disp('DEMO 4');
	m2D(9004,9004,100,100,0,20,1,1,1,1,0.000001,0.000001,200,2,2,0,20);
	pb3 =uicontrol(gcf,...
		'Style','push',...
		'Position',[.8 0 .2 .07],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','close all;dmalign2(''demos'')');
		
elseif strcmp(action,'demo5'),
	clf;
	disp('DEMO 5');
	fid=fopen('ref9005.mat','r+');		
	if (fid>0),
		fclose(fid);
	end;
	if (fid==-1),
   		disp('.....Generating the test set');
		%
		% 	generates a discrete ellipsoid
		%
		w=50;		%  size of voxel space
		h=50;
		d=50;
		a=8;		%  size of ellipsoid
		b=12;
		c=20;
		translimit=10;	% limits of random shifts ( rotlimit in degrees )
		rotlimit=10;
		x=max([a b c]);
		N=ceil(2*pi*x);	%  to make sure that ellipse does not have holes!!
		%
		[X,Y,Z] = sphere(N);
		%
		X=(X*a)+w/2;
		Y=(Y*b)+h/2;
		Z=(Z*c)+d/2;
		%
		nn=(N+1)^2;
		ref=[reshape(X,nn,1) reshape(Y,nn,1) reshape(Z,nn,1)];
		[n,m]=size(X);
		% 
		X1=reshape(X,n*m,1);
		Y1=reshape(Y,n*m,1);
		Z1=reshape(Z,n*m,1);
		%
		P=[X1 Y1 Z1];
		%
		%	rotate and translate randomly
		%
		trans=(rand(1,3)-0.5)*2*translimit;			%  voxels
		rot=(rand(1,3)-0.5)*2*rotlimit;  			%  degrees
		rot=rot*pi/180;
		R=getXYZ(rot(1),rot(2),rot(3)); R=R';
		%
		C=mean(P);
		P=getnewP(R,trans,P',C);
		%
		%  	save reference and test ellipsoid coordinates and transformation
		%
		test=P';
		save ref9005 ref;
   		save test9005 test;
    	%save desiredparam9005 trans rot;
	end;
	%
	m3D(9005,9005,50,50,50,0,100,1,0,1,0,0.000001,0.000001,200,2,2,0,20);
	pb3 =uicontrol(gcf,...
		'Style','push',...
		'Position',[.8 0 .2 .07],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','close all;dmalign2(''demos'')');
	
			
elseif strcmp(action,'start'),
	clf;
	fig = gcf;
	load alignpref;

	txt_dim = uicontrol(fig,...
		'Style','text',...
		'String','Dimensions:',...
		'Position',[.2 .9 .16 .05],...
		'Units','normalized');
	td_2d =  uicontrol(fig,...
		'Style','radio',...
		'String','2D',...
		'Position',[.36 .9 .1 .05],...
		'Units','normalized',...
		'Value',D2,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_2d = ui_handles1(1);'...
			'td_3d = ui_handles1(2);'...
			'txt_d = ui_handles1(5);'...
			'ed_d = ui_handles1(6);'...
			'set(td_2d,''Value'',1);'...
			'set(td_3d,''Value'',0);'...
			'set(txt_d,''Visible'',''off'');'...
			'set(ed_d,''Visible'',''off'');']);
	td_3d =  uicontrol(fig,...
		'Style','radio',...
		'String','3D',...
		'Position',[.46 .9 .1 .05],...
		'Units','normalized',...
		'Value',D3,...
		'CallBack',[...
		    'ui_handles1 =get(gcf,''Userdata'');'...
			'td_2d = ui_handles1(1);'...
			'td_3d = ui_handles1(2);'...
			'txt_d = ui_handles1(5);'...
			'ed_d = ui_handles1(6);'...
			'set(td_3d,''Value'',1);'...
			'set(td_2d,''Value'',0);'...
			'set(txt_d,''Visible'',''on'');'...
			'set(ed_d,''Visible'',''on'');']);
	txt_w = uicontrol(fig,...
		'Style','text',...
		'String','Width:',...
		'Position',[.4 .84 .08 .05],...
		'Units','normalized');
	ed_w = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(width),...
		'Position',[.48 .84 .1 .05],...
		'Units','normalized');
	txt_h = uicontrol(fig,...
		'Style','text',...
		'String','Height:',...
		'Position',[.6 .84 .08 .05],...
		'Units','normalized');
	ed_h = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(height),...
		'Position',[.68 .84 .1 .05],...
		'Units','normalized');
	txt_d = uicontrol(fig,...
		'Style','text',...
		'String','Depth:',...
		'Position',[.8 .84 .08 .05],...
		'Units','normalized');
	ed_d = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(depth),...
		'Position',[.88 .84 .1 .05],...
		'Units','normalized');
	if D2,
		set(txt_d,'Visible','off');
		set(ed_d,'Visible','off');	
	end;
	txt_trans = uicontrol(fig,...
		'Style','text',...
		'String','Transformation:      ',...
		'Position',[.2 .78 .23 .05],...
		'Units','normalized');
	td_rigid =  uicontrol(fig,...
		'Style','radio',...
		'String','Rigid Body',...
		'Position',[.5 .78 .16 .05],...
		'Units','normalized',...
		'Value',RIG,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_rigid = ui_handles1(7);'...
			'td_affine = ui_handles1(8);'...
			'txt_1 = ui_handles1(35);'...
			'ed_1 = ui_handles1(26);'...
			'txt_2 = ui_handles1(36);'...
			'set(td_rigid,''Value'',1);'...
			'set(td_affine,''Value'',0);'...
			'set(txt_1,''Visible'',''on'');'...
			'set(ed_1,''Visible'',''on'');'...
			'set(txt_2,''String'','' Translation threshold (pix/vox):'')']);
	td_affine =  uicontrol(fig,...
		'Style','radio',...
		'String','Affine',...
		'Position',[.7 .78 .16 .05],...
		'Units','normalized',...
		'Value',AFF,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_rigid = ui_handles1(7);'...
			'td_affine = ui_handles1(8);'...
			'txt_1 = ui_handles1(35);'...
			'ed_1 = ui_handles1(26);'...
			'txt_2 = ui_handles1(36);'...
			'set(td_rigid,''Value'',0);'...
			'set(td_affine,''Value'',1);'...
			'set(txt_1,''Visible'',''off'');'...
			'set(ed_1,''Visible'',''off'');'...
			'set(txt_2,''String'',''  Affine threshold :'')']);
	
	txt_rand = uicontrol(fig,...
		'Style','text',...
		'String','Sampling:  ',...
		'Position',[.2 .71 .16 .05],...
		'Units','normalized');
		
	td_rno =  uicontrol(fig,...
		'Style','radio',...
		'String','Use all points',...
		'Position',[.4 .71 .18 .05],...
		'Units','normalized',...
		'Value',randsamp==0,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_rno = ui_handles1(9);'...
			'td_ryes = ui_handles1(10);'...
			'ed_randno = ui_handles1(11);'...
			'txt_randno = ui_handles1(12);'...
			'set(td_rno,''Value'',1);'...
			'set(td_ryes,''Value'',0);'...
			'set(ed_randno,''Visible'',''off'');'...
			'set(txt_randno,''Visible'',''off'');']);
	td_ryes =  uicontrol(fig,...
		'Style','radio',...
		'String','Sample',...
		'Position',[.6 .71 .16 .05],...
		'Units','normalized',...
		'Value',randsamp>0,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_rno = ui_handles1(9);'...
			'td_ryes = ui_handles1(10);'...
			'ed_randno = ui_handles1(11);'...
			'txt_randno = ui_handles1(12);'...
			'set(td_rno,''Value'',0);'...
			'set(td_ryes,''Value'',1);'...
			'set(ed_randno,''Visible'',''on'');'...
			'set(txt_randno,''Visible'',''on'');']);
	ed_randno = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(randsamp),...
		'Position',[.76 .71 .10 .05],...
		'Units','normalized');
	txt_randno = uicontrol(fig,...
		'Style','text',...
		'String','random points',...
		'Position',[.86 .71 .16 .05],...
		'Units','normalized');
	if randsamp==0
		set(ed_randno,'Visible','off');
		set(txt_randno,'Visible','off');
	end;
	
	txt_batch = uicontrol(fig,...
		'Style','text',...
		'String','Alignment:',...
		'Position',[.2 .65 .16 .05],...
		'Units','normalized');
		
	td_single =  uicontrol(fig,...
		'Style','radio',...
		'String','Single',...
		'Position',[.4 .65 .18 .05],...
		'Units','normalized',...
		'Value',batchflag==0,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_single = ui_handles1(13);'...
			'td_batch = ui_handles1(14);'...
			'txt_start = ui_handles1(15);'...
			'ed_start = ui_handles1(16);'...
			'txt_end = ui_handles1(17);'...
			'ed_end = ui_handles1(18);'...
			'txt_12 = ui_handles1(33);'...
			'ed_12 = ui_handles1(34);'...
			'set(td_single,''Value'',1);'...
			'set(td_batch,''Value'',0);'...
			'set(txt_start,''Visible'',''off'');'...
			'set(ed_start,''Visible'',''off'');'...
			'set(txt_end,''Visible'',''off'');'...
			'set(ed_end,''Visible'',''off'');'...
			'set(txt_12,''Visible'',''on'');'...
			'set(ed_12,''Visible'',''on'');']);
			
	td_batch =  uicontrol(fig,...
		'Style','radio',...
		'String','Batch',...
		'Position',[.6 .65 .16 .05],...
		'Units','normalized',...
		'Value',batchflag==1,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_single = ui_handles1(13);'...
			'td_batch = ui_handles1(14);'...
			'txt_start = ui_handles1(15);'...
			'ed_start = ui_handles1(16);'...
			'txt_end = ui_handles1(17);'...
			'ed_end = ui_handles1(18);'...
			'txt_12 = ui_handles1(33);'...
			'ed_12 = ui_handles1(34);'...
			'set(td_single,''Value'',0);'...
			'set(td_batch,''Value'',1);'...
			'set(txt_start,''Visible'',''on'');'...
			'set(ed_start,''Visible'',''on'');'...
			'set(txt_end,''Visible'',''on'');'...
			'set(ed_end,''Visible'',''on'');'...
			'set(txt_12,''Visible'',''off'');'...
			'set(ed_12,''Visible'',''off'');']);
	txt_start = uicontrol(fig,...
		'Style','text',...
		'String','Start at:',...
		'Position',[.8 .65 .1 .05],...
		'Units','normalized');
	ed_start = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(bst),...
		'Position',[.9 .65 .1 .05],...
		'Units','normalized');
	txt_end = uicontrol(fig,...
		'Style','text',...
		'String','End at:',...
		'Position',[.8 .58 .1 .05],...
		'Units','normalized');
	ed_end = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(bend),...
		'Position',[.9 .58 .1 .05],...
		'Units','normalized');
	
	if batchflag==0
		set(txt_start,'Visible','off');
		set(ed_start,'Visible','off');
		set(txt_end,'Visible','off');
		set(ed_end,'Visible','off');
	end;	
	
	txt_dmloc = uicontrol(fig,...
		'Style','text',...
		'String','Distancemap location:     ',...
		'Position',[.2 .52 .26 .05],...
		'Units','normalized');
	td_DMram =  uicontrol(fig,...
		'Style','radio',...
		'String','ON RAM',...
		'Position',[.5 .52 .16 .05],...
		'Units','normalized',...
		'Value',dmloc==0,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_DMram = ui_handles1(19);'...
			'td_DMdisk = ui_handles1(20);'...
			'set(td_DMram,''Value'',1);'...
			'set(td_DMdisk,''Value'',0);']);
	td_DMdisk =  uicontrol(fig,...
		'Style','radio',...
		'String','ON DISK',...
		'Position',[.7 .52 .16 .05],...
		'Units','normalized',...
		'Value',dmloc==1,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_DMram = ui_handles1(19);'...
			'td_DMdisk = ui_handles1(20);'...
			'set(td_DMram,''Value'',0);'...
			'set(td_DMdisk,''Value'',1);']);
	
	cb_1 =  uicontrol(fig,...
		'Style','checkbox',...
		'String','Replace outliers.',...
		'Position',[.6 .44 .34 .05],...
		'Units','normalized',...
		'Value',remout);
	cb_2 =  uicontrol(fig,...
		'Style','checkbox',...
		'String','Clear DM after alignment.',...
		'Position',[.6 .39 .34 .05],...
		'Units','normalized',...
		'Value',exitflag);
	cb_3 =  uicontrol(fig,...
		'Style','checkbox',...
		'String','Save parameters and errors.',...
		'Position',[.6 .34 .34 .05],...
		'Units','normalized',...
		'Value',saveparam);
	cb_4 =  uicontrol(fig,...
		'Style','checkbox',...
		'String','   Plot during the simulation.',...
		'Position',[.6 .29 .34 .05],...
		'Units','normalized',...
		'Value',plotflag);
	cb_5 =  uicontrol(fig,...
		'Style','checkbox',...
		'String','   Plot at the end.',...
		'Position',[.6 .24 .34 .05],...
		'Units','normalized',...
		'Value',plotfinal);
		
	txt_1 = uicontrol(fig,...
		'Horizontalalignment','left',...
		'Style','text',...
		'String',' Angular threshold (in radian): ',...
		'Position',[.12 .44 .34 .05],...
		'Units','normalized');
	ed_1 = uicontrol(fig,...
		'Style','edit',...
		'string',num2str(angtresh),...
		'Position',[.46 .44 .1 .05],...
		'Units','normalized');
		
	txt_2 = uicontrol(fig,...
		'Horizontalalignment','left',...
		'Style','text',...
		'String',' Translation threshold (pix/vox):',...
		'Position',[.12 .39 .34 .05],...
		'Units','normalized');
	ed_2 = uicontrol(fig,...
		'Style','edit',...
		'string',num2str(transtresh),...
		'Position',[.46 .39 .1 .05],...
		'Units','normalized');
	
	if AFF,
		set(txt_1,'Visible','off');
		set(ed_1,'Visible','off');
		set(txt_2,'String','  Affine threshold :');
	end;		
				
	txt_3 = uicontrol(fig,...
		'Horizontalalignment','left',...
		'Style','text',...
		'String',' Maximum iteration no :',...
		'Position',[.12 .34 .34 .05],...
		'Units','normalized');
	ed_3 = uicontrol(fig,...
		'Style','edit',...
		'string',num2str(maxiter),...
		'Position',[.46 .34 .1 .05],...
		'Units','normalized');
		
	txt_5 = uicontrol(fig,...
		'Horizontalalignment','left',...
		'Style','text',...
		'String',' Parameter 1 (alpha):    ',...
		'Position',[.12 .29 .34 .05],...
		'Units','normalized');
	ed_5 = uicontrol(fig,...
		'Style','edit',...
		'string',num2str(alfa),...
		'Position',[.46 .29 .1 .05],...
		'Units','normalized');
		
	txt_6 = uicontrol(fig,...
		'Horizontalalignment','left',...
		'Style','text',...
		'String',' Parameter 2 (vv):',...
		'Position',[.12 .24 .34 .05],...
		'Units','normalized');
	ed_6 = uicontrol(fig,...
		'Style','edit',...
		'string',num2str(vv),...
		'Position',[.46 .24 .1 .05],...
		'Units','normalized');
		
	txt_9 = uicontrol(fig,...
		'Horizontalalignment','left',...
		'Style','text',...
		'String',' Maximum case 3 loop search :',...
		'Position',[.12 .19 .34 .05],...
		'Units','normalized');
	ed_9 = uicontrol(fig,...
		'Style','edit',...
		'string',num2str(maxiterloop),...
		'Position',[.46 .19 .1 .05],...
		'Units','normalized');


	txt_11 = uicontrol(fig,...
		'Style','text',...
		'String','Refno:',...
		'Position',[.25 .59 .08 .05],...
		'Units','normalized');
	ed_11 = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(refno),...
		'Position',[.33 .59 .1 .05],...
		'Units','normalized');		
	txt_12 = uicontrol(fig,...
		'Style','text',...
		'String','Testno:',...
		'Position',[.45 .59 .08 .05],...
		'Units','normalized');
	ed_12 = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(testno),...
		'Position',[.53 .59 .1 .05],...
		'Units','normalized');
				
   	if batchflag,
		set(txt_12,'Visible','off');
		set(ed_12,'Visible','off');
	end;
		
	pb3 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .05 .2 .1],...
		'Units','normalized',...
		'String','SAVE AS DEFAULT',...
		'CallBack','dmalign2(''saveoptions'')');

	pb4 =uicontrol(fig,...
		'Style','push',...
		'Position',[.7 .05 .2 .1],...
		'Units','normalized',...
		'String','ALIGN NOW',...
		'CallBack','dmalign2(''align'')');

	pb5 =uicontrol(fig,...
		'Style','push',...
		'Position',[.45 .05 .2 .1],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','dmalign2(''initialize'')');
				
    dumy(1:10)=[td_2d td_3d ed_w ed_h txt_d ed_d td_rigid td_affine td_rno td_ryes ];
				%1	   2	 3	   4   5	  6		7		8		 9      10			
    dumy(11:18)= [ed_randno txt_randno td_single  td_batch txt_start ed_start txt_end ed_end ];
					%11  		12			13			14        15         16		17		18		
    dumy(19:36)= [td_DMram td_DMdisk cb_1 cb_2 cb_3 cb_4 cb_5 ed_1 ed_2 ed_3 ed_5 ed_6 ed_9 ed_11 txt_12 ed_12 txt_1 txt_2];
				%  19      10        21    22   23   24   25	26	27	28	  29  30	31	32     33     34   35      36
	set(fig,'Userdata',dumy);
	
elseif strcmp(action,'saveoptions'),
	ui_handles1 = get(gcf,'Userdata');
	td_2d = ui_handles1(1);
	td_3d = ui_handles1(2);
	ed_w = ui_handles1(3);
	ed_h = ui_handles1(4);
	ed_d = ui_handles1(6);
	td_rigid = ui_handles1(7);
	td_affine = ui_handles1(8);
	td_rno = ui_handles1(9);
	td_ryes = ui_handles1(10);
	ed_randno = ui_handles1(11);
	td_single  = ui_handles1(13);
	td_batch = ui_handles1(14);
	ed_start = ui_handles1(16);
	ed_end = ui_handles1(18);
	td_DMram = ui_handles1(19);
	td_DMdisk = ui_handles1(20);
	cb_1 = ui_handles1(21);
	cb_2 = ui_handles1(22);
	cb_3 = ui_handles1(23);
	cb_4 = ui_handles1(24);
	cb_5 = ui_handles1(25);
	ed_1 = ui_handles1(26);
	ed_2 = ui_handles1(27);
	ed_3 = ui_handles1(28);
	ed_5 = ui_handles1(29);
	ed_6 = ui_handles1(30);
	ed_9 = ui_handles1(31);
	ed_11 = ui_handles1(32);
	ed_12 = ui_handles1(34);
		
	D2 = get(td_2d,'Value');
	D3 = get(td_3d,'Value');
	width = str2num(get(ed_w,'String'));
	height = str2num(get(ed_h,'String'));
	depth = str2num(get(ed_d,'String'));
	RIG = get(td_rigid,'Value');
	AFF = get(td_affine,'Value');
 	if get(td_rno,'Value');
 		randsamp = 0;
	else,
		randsamp = str2num(get(ed_randno,'String'));
	end;
	batchflag = get(td_batch,'Value');
	bst = str2num(get(ed_start,'String'));
	bend = str2num(get(ed_end,'String'));
	if get(td_DMram,'Value');
		dmloc = 0;
	else,
		dmloc = 1;
	end;
	remout = get(cb_1,'Value');
	exitflag = get(cb_2,'Value');
	saveparam = get(cb_3,'Value');
	plotflag = get(cb_4,'Value');
	plotfinal = get(cb_5,'Value');
	
	angtresh = str2num(get(ed_1,'String'));
	transtresh = str2num(get(ed_2,'String'));
	maxiter = str2num(get(ed_3,'String'));
	alfa = str2num(get(ed_5,'String'));
	vv  = str2num(get(ed_6,'String'));
	maxiterloop = str2num(get(ed_9,'String'));
	
	refno = str2num(get(ed_11,'String'));
	testno = str2num(get(ed_12,'String'));
		
	save alignpref D2 D3 width height depth RIG AFF randsamp remout...
	     dmloc exitflag saveparam plotflag plotfinal batchflag bst bend...
		 angtresh transtresh maxiter alfa vv refno testno maxiterloop;
	

elseif strcmp(action,'grayscale'),
	clf;
	fig = gcf;
	load alignpref;

	txt_intro1 = uicontrol(fig,...
		'Style','text',...
		'String',' This part converts a grayscale raw test image (or image set)',...
		'Position',[.2 .85 .7 .05],...
		'Units','normalized');
	
	txt_intro2 = uicontrol(fig,...
		'Style','text',...
		'String',' using previously calculated transformation parameters ',...
		'Position',[.2 .8 .7 .05],...
		'Units','normalized');
		
	txt_dim = uicontrol(fig,...
		'Style','text',...
		'String','Dimensions:',...
		'Horizontalalignment','left',...
		'Position',[.2 .7 .16 .05],...
		'Units','normalized');
		
	td_2d =  uicontrol(fig,...
		'Style','radio',...
		'String','2D',...
		'Position',[.4 .7 .1 .05],...
		'Units','normalized',...
		'Value',D2,...
		'CallBack',[...
		    'ui_handles1 = get(gcf,''Userdata'');'...
			'td_2d = ui_handles1(1);'...
			'td_3d = ui_handles1(2);'...
			'txt_d = ui_handles1(5);'...
			'ed_d = ui_handles1(6);'...
			'set(td_2d,''Value'',1);'...
			'set(td_3d,''Value'',0);'...
			'set(txt_d,''Visible'',''off'');'...
			'set(ed_d,''Visible'',''off'');']);
			
	td_3d =  uicontrol(fig,...
		'Style','radio',...
		'String','3D',...
		'Position',[.5 .7 .1 .05],...
		'Units','normalized',...
		'Value',D3,...
		'CallBack',[...
		    'ui_handles1 =get(gcf,''Userdata'');'...
			'td_2d = ui_handles1(1);'...
			'td_3d = ui_handles1(2);'...
			'txt_d = ui_handles1(5);'...
			'ed_d = ui_handles1(6);'...
			'set(td_3d,''Value'',1);'...
			'set(td_2d,''Value'',0);'...
			'set(txt_d,''Visible'',''on'');'...
			'set(ed_d,''Visible'',''on'');']);
	txt_w = uicontrol(fig,...
		'Style','text',...
		'String','Width:',...
		'Position',[.4 .6 .08 .05],...
		'Units','normalized');
	ed_w = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(width),...
		'Position',[.48 .6 .06 .05],...
		'Units','normalized');
	txt_h = uicontrol(fig,...
		'Style','text',...
		'String','Height:',...
		'Position',[.6 .6 .08 .05],...
		'Units','normalized');
	ed_h = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(height),...
		'Position',[.68 .6 .06 .05],...
		'Units','normalized');
	txt_d = uicontrol(fig,...
		'Style','text',...
		'String','Depth:',...
		'Position',[.8 .6 .08 .05],...
		'Units','normalized');
	ed_d = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(depth),...
		'Position',[.88 .6 .06 .05],...
		'Units','normalized');
	if D2,
		set(txt_d,'Visible','off');
		set(ed_d,'Visible','off');	
	end;
	
	txt_start = uicontrol(fig,...
		'Style','text',...
		'String','Start at:',...
		'Horizontalalignment','left',...
		'Position',[.2 .5 .2 .05],...
		'Units','normalized');
		
	ed_start = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(bst),...
		'Position',[.4 .5 .1 .05],...
		'Units','normalized');
		
		
	txt_end = uicontrol(fig,...
		'Style','text',...
		'String','End at:',...
		'Position',[.6 .5 .2 .05],...
		'Units','normalized');
	
	ed_end = uicontrol(fig,...
		'Style','edit',...
		'string',int2str(bst),...
		'Position',[.8 .5 .1 .05],...
		'Units','normalized');
		
	
	txt_fname = uicontrol(fig,...
		'Style','text',...
		'String','Filename:',...
		'Horizontalalignment','left',...
		'Position',[.2 .4 .2 .05],...
		'Units','normalized');
		
	ed_fname = uicontrol(fig,...
		'Style','edit',...
		'string','name',...
		'Position',[.4 .4 .2 .05],...
		'Units','normalized');
	
	pb4 =uicontrol(fig,...
		'Style','push',...
		'Position',[.2 .05 .3 .1],...
		'Units','normalized',...
		'String','CORRECT GRAYSCALE',...
		'CallBack','dmalign2(''correct'')');

	
	pb5 =uicontrol(fig,...
		'Style','push',...
		'Position',[.6 .05 .3 .1],...
		'Units','normalized',...
		'String','RETURN',...
		'CallBack','dmalign2(''initialize'')');
					
    	dumy(1:10)=    [	td_2d 	td_3d 	ed_w 	ed_h 	txt_d 	ed_d  	ed_start pb4 	ed_fname ed_end ];
			%1	2	 3	4   	5	6	7	 8	9	  10

	set(fig,'Userdata',dumy);

elseif strcmp(action,'correct'),
	ui_handles1 = get(gcf,'Userdata');
	td_2d = ui_handles1(1);
	td_3d = ui_handles1(2);
	ed_w = ui_handles1(3);
	ed_h = ui_handles1(4);
	ed_d = ui_handles1(6);
	ed_start = ui_handles1(7);
	ed_fname = ui_handles1(9);
	ed_end = ui_handles1(10);
	
	D2 = get(td_2d,'Value');
	D3 = get(td_3d,'Value');
	width = str2num(get(ed_w,'String'));
	height = str2num(get(ed_h,'String'));
	depth = str2num(get(ed_d,'String'));
	testst = str2num(get(ed_start,'String'));
	tested = str2num(get(ed_end,'String'));
	
	fname = get(ed_fname,'String');
	
	if D2,
	   	for testno=testst:tested,
			gsinterp2(width,height,testno,fname);
		end;
	end;
	if D3,
	   	for testno=testst:tested,	
			gsinterp3(width,height,depth,testno,fname);
		end;
	end;
	
elseif strcmp(action,'align'),
	ui_handles1 = get(gcf,'Userdata');
	td_2d = ui_handles1(1);
	td_3d = ui_handles1(2);
	ed_w = ui_handles1(3);
	ed_h = ui_handles1(4);
	ed_d = ui_handles1(6);
	td_rigid = ui_handles1(7);
	td_affine = ui_handles1(8);
	td_rno = ui_handles1(9);
	td_ryes = ui_handles1(10);
	ed_randno = ui_handles1(11);
	td_single  = ui_handles1(13);
	td_batch = ui_handles1(14);
	ed_start = ui_handles1(16);
	ed_end = ui_handles1(18);
	td_DMram = ui_handles1(19);
	td_DMdisk = ui_handles1(20);
	cb_1 = ui_handles1(21);
	cb_2 = ui_handles1(22);
	cb_3 = ui_handles1(23);
	cb_4 = ui_handles1(24);
	cb_5 = ui_handles1(25);
	ed_1 = ui_handles1(26);
	ed_2 = ui_handles1(27);
	ed_3 = ui_handles1(28);
	ed_5 = ui_handles1(29);
	ed_6 = ui_handles1(30);
	ed_9 = ui_handles1(31);
	ed_11 = ui_handles1(32);
	ed_12 = ui_handles1(34);
		
	D2 = get(td_2d,'Value');
	D3 = get(td_3d,'Value');
	width = str2num(get(ed_w,'String'));
	height = str2num(get(ed_h,'String'));
	depth = str2num(get(ed_d,'String'));
	RIG = get(td_rigid,'Value');
	AFF = get(td_affine,'Value');
 	if get(td_rno,'Value');
 		randsamp = 0;
	else,
		randsamp = str2num(get(ed_randno,'String'));
	end;
	batchflag = get(td_batch,'Value');
	bst = str2num(get(ed_start,'String'));
	bend = str2num(get(ed_end,'String'));
	if get(td_DMram,'Value');
		dmloc = 0;
	else,
		dmloc = 1;
	end;
	remout = get(cb_1,'Value');
	exitflag = get(cb_2,'Value');
	saveparam = get(cb_3,'Value');
	plotflag = get(cb_4,'Value');
	plotfinal = get(cb_5,'Value');
	
	angtresh = str2num(get(ed_1,'String'));
	transtresh = str2num(get(ed_2,'String'));
	maxiter = str2num(get(ed_3,'String'));
	alfa = str2num(get(ed_5,'String'));
	vv  = str2num(get(ed_6,'String'));
	maxiterloop = str2num(get(ed_9,'String'));
	
	refno = str2num(get(ed_11,'String'));
	testno = str2num(get(ed_12,'String'));
	%
	%
	%	
	if D2 & RIG,
		if ~batchflag
			m2D(testno,refno,width,height,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,angtresh,transtresh,maxiter,alfa,vv,dmloc,maxiterloop);
		end;
		if batchflag,
			for ii=bst:bend,
				m2D(ii,refno,width,height,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,angtresh,transtresh,maxiter,alfa,vv,dmloc,maxiterloop);
			end;
		end;
	end;
	if D2 & AFF,
		if ~batchflag,
			m2Daff(testno,refno,width,height,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,transtresh,maxiter,alfa,vv,dmloc,maxiterloop);
		end;
		if batchflag,
			for ii=bst:bend;
				m2Daff(ii,refno,width,height,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,transtresh,maxiter,alfa,vv,dmloc,maxiterloop);
			end;
		end;		
	end;
	if D3 & RIG,
		if plotflag,
			disp('Plotting during 3D alignmnent is not recommended.')
			plotflag = 0;
		end;
		if ~batchflag
			m3D(testno,refno,width,height,depth,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,angtresh,transtresh,maxiter,alfa,vv,dmloc,maxiterloop)
		end;
		if batchflag,
			for ii=bst:bend;
				m3D(ii,refno,width,height,depth,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,angtresh,transtresh,maxiter,alfa,vv,dmloc,maxiterloop)
			end;
		end;		
	end;
	if D3 & AFF,
		if plotflag,
			disp('Plotting during 3D alignmnent is not recommended.')
			plotflag = 0;
		end;
		if ~batchflag
			m3Daff(testno,refno,width,height,depth,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,transtresh,maxiter,alfa,vv,dmloc,maxiterloop);
		end;
		if batchflag,
			for ii=bst:bend;
				m3Daff(ii,refno,width,height,depth,saveparam,randsamp,exitflag,plotflag,plotfinal,remout,transtresh,maxiter,alfa,vv,dmloc,maxiterloop);
			end;
		end;		
	end;
end;

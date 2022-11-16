function varargout = guisufnec(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guisufnec_OpeningFcn, ...
                   'gui_OutputFcn',  @guisufnec_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function guisufnec_OpeningFcn(hObject, eventdata, handles, varargin)

axes(handles.adj_mat)
imagesc(rand(5))

set(handles.plot1_button,'Enable','off')
set(handles.plot2_button,'Enable','off')
set(handles.show_net,'Enable','off')
set(handles.SNb1,'Enable','off')
set(handles.slider1,'Enable','off')
set(handles.slider2,'Enable','off')
set(handles.SNb2,'Enable','off')
set(handles.slider3,'Enable','off')
set(handles.slider4,'Enable','off')

handles.output = hObject;
guidata(hObject, handles);

function varargout = guisufnec_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function gen_net_Callback(hObject, eventdata, handles)

handles.X = gen_net;

set(handles.plot1_button,'Enable','on')
set(handles.plot2_button,'Enable','on')
set(handles.show_net,'Enable','on')
% set(handles.SNb1,'Enable','on')
% set(handles.SNb2,'Enable','on')

set(handles.SNb1,'Value',1)
set(handles.SNb2,'Value',1)
set(handles.lb1,'String',{'Select plot 1'})
set(handles.lb2,'String',{'Select plot 2'})

handles.SN1=1;
handles.SN2=1;

guidata(hObject, handles);

function editnet_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_net_Callback(hObject, eventdata, handles)

for n=1:5
    handles.X(:,n) = handles.X(:,n) + 10*n;
end
figure(1)
plot(handles.X)

function plot1_button_Callback(hObject, eventdata, handles)

[c,f]=ginput(1);
f=round(f);
c=round(c);
n1=f; 
n2=c;

set(handles.SNb1,'Enable','on')
set(handles.slider1,'Enable','on')
set(handles.slider2,'Enable','on')

set(handles.SNb1,'Value',1)
handles.SN1=1;

a=handles.X(:,f);
b=handles.X(:,c);

Npts = 1000;

% scale the signal:
a = a /max(a); a=a-min(a);
b = b/max(b); b=b-min(b);

tha=linspace(min(a),max(a),50);
thb=linspace(min(b),max(b),50);

% Binarization part:  simple thresholding
atmp=ones(Npts,1);
btmp=ones(Npts,1);

atmp(a<tha(1))=0;
btmp(b<thb(1))=0;
    
axes(handles.events1_1)        
stem(atmp,'g'); hold on,
plot(a); 
line([0 Npts],[tha(1) tha(1)]); hold off

axes(handles.events2_1)
stem(btmp,'g'); hold on,
plot(b,'r');
line([0 Npts],[thb(1) thb(1)]); hold off

Nab1=ncsty(atmp,btmp);
Sab1=sfcy(atmp,btmp);

for i=1:50
    for j=1:50
    
    % Binarization part:  simple thresholding
    atmp=ones(Npts,1);
    btmp=ones(Npts,1);
    
    atmp(a<tha(i))=0;
    btmp(b<thb(j))=0;
        
    t=[sum(atmp&btmp),sum(atmp&~btmp);sum(~atmp&btmp),sum(~atmp&~btmp)];
    p_y_given_x=t(1,1)/sum(t(:,1));
    p_y_given_not_x=t(1,2)/sum(t(:,2));
    PNS=p_y_given_x-p_y_given_not_x;
   
    %Necessity
    Nab(j,i)=PNS/p_y_given_x;
    if isnan(Nab(j,i)) | isinf(Nab(j,i)) | Nab(j,i)<0
        Nab(j,i)=0;
    end
    %Sufficiency
    Sab(j,i)=PNS/(1-p_y_given_not_x);
    if isnan(Sab(j,i)) | isinf(Sab(j,i)) | Sab(j,i)<0 
        Sab(j,i)=0;
    end
    
    end
end

[X,Y]=meshgrid(tha,thb);

if handles.SN1
    axes(handles.sufnec1)
    surf(X,Y,Nab);
    hold on
    line([tha(1) tha(end)],[thb(1) thb(1)],[1.2 1.2]);
    line([tha(1) tha(1)],[thb(1) thb(end)],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(n1)])
    ylabel(['threshold node ', num2str(n2)])
    zlabel('Necessity')
    axis([min(a) max(a) min(b) max(b) 0 1.2])
    colorbar
else
    axes(handles.sufnec1)
    surf(X,Y,Sab);
    hold on
    line([tha(1) tha(end)],[thb(1) thb(1)],[1.2 1.2]);
    line([tha(1) tha(1)],[thb(1) thb(end)],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(n1)])
    ylabel(['threshold node ', num2str(n2)])
    zlabel('Sufficiency')
    axis([min(a) max(a) min(b) max(b) 0 1.2])
    colorbar
end

handles.n1g1=n1;
handles.n2g1=n2;
handles.atmp1=atmp;
handles.btmp1=btmp;
handles.tha1=min(a);
handles.thb1=min(b);

handles.plot1_3d.X=X;
handles.plot1_3d.Y=Y;
handles.plot1_3d.N=Nab;
handles.plot1_3d.S=Sab;
handles.plot1_3d.tha=[tha(1) tha(end)];
handles.plot1_3d.thb=[thb(1) thb(end)];
handles.plot1_3d.a=[min(a) max(a)];
handles.plot1_3d.b=[min(b) max(b)];

set(handles.slider1,'value',0)
set(handles.slider2,'value',0)

cont1={['Nodes ',num2str(handles.n1g1),'-',num2str(handles.n2g1)],...
    ['N: ',num2str(floor(Nab1*10000)/10000),' / S: ',num2str(floor(Sab1*10000)/10000)],...
    ['Total N = ',num2str(sum(handles.plot1_3d.N(:)))],...
    ['Total S = ',num2str(sum(handles.plot1_3d.S(:)))],...
    ['theshold node ',num2str(handles.n1g1),' = ',num2str(floor(handles.tha1*10000)/10000)],...
    ['theshold node ',num2str(handles.n2g1),' = ',num2str(floor(handles.thb1*10000)/10000)]};
set(handles.lb1,'String',cont1)

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)

guidata(hObject, handles);

function slider1_Callback(hObject, eventdata, handles)

axes(handles.sufnec1)
[az,el] = view;

Npts=1000;
a=handles.X(:,handles.n1g1);
b=handles.X(:,handles.n2g1);

% scale the signal:
a=a/max(a); 
a=a-min(a);

th=get(hObject,'Value')*(max(a)+eps);

% Binarization part: simple thresholding
atmp=ones(Npts,1);
atmp(a<th) = 0;
    
axes(handles.events1_1)        
stem(atmp,'g'); 
hold on,
plot(a); 
line([0 Npts],[th th]); 
hold off
      
handles.tha1=th;
handles.atmp1=atmp;

handles.Nab1=ncsty(atmp,handles.btmp1);
handles.Sab1=sfcy(atmp,handles.btmp1);

cont1={['Nodes ',num2str(handles.n1g1),'-',num2str(handles.n2g1)],...
    ['N: ',num2str(floor(handles.Nab1*10000)/10000),' / S: ',num2str(floor(handles.Sab1*10000)/10000)],...
    ['Total N = ',num2str(sum(handles.plot1_3d.N(:)))],...
    ['Total S = ',num2str(sum(handles.plot1_3d.S(:)))],...
    ['theshold node ',num2str(handles.n1g1),' = ',num2str(floor(handles.tha1*10000)/10000)],...
    ['theshold node ',num2str(handles.n2g1),' = ',num2str(floor(handles.thb1*10000)/10000)]};
set(handles.lb1,'String',cont1)

if handles.SN1
    axes(handles.sufnec1)
    surf(handles.plot1_3d.X,handles.plot1_3d.Y,handles.plot1_3d.N);
    hold on
    line([th th],[handles.plot1_3d.thb(1) handles.plot1_3d.thb(2)],[1.2 1.2]);
    line([handles.plot1_3d.tha(1) handles.plot1_3d.tha(2)],[handles.thb1 handles.thb1],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g1)])
    ylabel(['threshold node ', num2str(handles.n2g1)])
    zlabel('Necessity')
    axis([handles.plot1_3d.a(1) handles.plot1_3d.a(2) handles.plot1_3d.b(1)...
        handles.plot1_3d.b(2) 0 1.2])
    colorbar
else
    axes(handles.sufnec1)
    surf(handles.plot1_3d.X,handles.plot1_3d.Y,handles.plot1_3d.S);
    hold on
    line([th th],[handles.plot1_3d.thb(1) handles.plot1_3d.thb(2)],[1.2 1.2]);
    line([handles.plot1_3d.tha(1) handles.plot1_3d.tha(2)],[handles.thb1 handles.thb1],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g1)])
    ylabel(['threshold node ', num2str(handles.n2g1)])
    zlabel('Sufficiency')
    axis([handles.plot1_3d.a(1) handles.plot1_3d.a(2) handles.plot1_3d.b(1)...
        handles.plot1_3d.b(2) 0 1.2])
    colorbar
end
    
guidata(hObject, handles);

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)

view([az,el]);

function slider2_Callback(hObject, eventdata, handles)

axes(handles.sufnec1)
[az,el] = view;

Npts=1000;
a=handles.X(:,handles.n1g1);
b=handles.X(:,handles.n2g1);

% scale the signal:
b=b/max(b); 
b=b-min(b);

th=get(hObject,'Value')*(max(b)+eps);

% Binarization part: simple thresholding
btmp = ones(Npts,1);
btmp(b<th) = 0;
    
axes(handles.events2_1)        
stem(btmp,'g'); 
hold on,
plot(b); 
line([0 Npts],[th th]); 
hold off
      
handles.thb1=th;
handles.btmp1=btmp;

handles.Nab1=ncsty(handles.atmp1,btmp);
handles.Sab1=sfcy(handles.atmp1,btmp);

cont1={['Nodes ',num2str(handles.n1g1),'-',num2str(handles.n2g1)],...
    ['N: ',num2str(floor(handles.Nab1*10000)/10000),' / S: ',num2str(floor(handles.Sab1*100)/100)],...
    ['Total N = ',num2str(sum(handles.plot1_3d.N(:)))],...
    ['Total S = ',num2str(sum(handles.plot1_3d.S(:)))],...
    ['theshold node ',num2str(handles.n1g1),' = ',num2str(floor(handles.tha1*10000)/10000)],...
    ['theshold node ',num2str(handles.n2g1),' = ',num2str(floor(handles.thb1*10000)/10000)]};
set(handles.lb1,'String',cont1)

if handles.SN1
    axes(handles.sufnec1)
    surf(handles.plot1_3d.X,handles.plot1_3d.Y,handles.plot1_3d.N);
    hold on
    line([handles.tha1 handles.tha1],[handles.plot1_3d.thb(1) handles.plot1_3d.thb(2)],[1.2 1.2]);
    line([handles.plot1_3d.tha(1) handles.plot1_3d.tha(2)],[th th],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g1)])
    ylabel(['threshold node ', num2str(handles.n2g1)])
    zlabel('Necessity')
    axis([handles.plot1_3d.a(1) handles.plot1_3d.a(2) handles.plot1_3d.b(1)...
        handles.plot1_3d.b(2) 0 1.2])
    colorbar
else
    axes(handles.sufnec1)
    surf(handles.plot1_3d.X,handles.plot1_3d.Y,handles.plot1_3d.S);
    hold on
    line([handles.tha1 handles.tha1],[handles.plot1_3d.thb(1) handles.plot1_3d.thb(2)],[1.2 1.2]);
    line([handles.plot1_3d.tha(1) handles.plot1_3d.tha(2)],[th th],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g1)])
    ylabel(['threshold node ', num2str(handles.n2g1)])
    zlabel('Sufficiency')
    axis([handles.plot1_3d.a(1) handles.plot1_3d.a(2) handles.plot1_3d.b(1)...
        handles.plot1_3d.b(2) 0 1.2])
    colorbar
end

guidata(hObject, handles);

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)
view([az,el]);

function plot2_button_Callback(hObject, eventdata, handles)

[c,f]=ginput(1);
f=round(f);
c=round(c);
n1=f; 
n2=c;

set(handles.SNb2,'Enable','on')
set(handles.slider3,'Enable','on')
set(handles.slider4,'Enable','on')

set(handles.SNb2,'Value',1)
handles.SN2=1;

a=handles.X(:,f);
b=handles.X(:,c);

Npts = 1000;

% scale the signal:
a = a /max(a); a=a-min(a);
b = b/max(b); b=b-min(b);

tha=linspace(min(a),max(a),50);
thb=linspace(min(b),max(b),50);

% Binarization part:  simple thresholding
atmp=ones(Npts,1);
btmp=ones(Npts,1);

atmp(a<tha(1))=0;
btmp(b<thb(1))=0;
    
axes(handles.events1_2)        
stem(atmp,'g'); hold on,
plot(a); 
line([0 Npts],[tha(1) tha(1)]); hold off

axes(handles.events2_2)
stem(btmp,'g'); hold on,
plot(b,'r');
line([0 Npts],[thb(1) thb(1)]); hold off

Nab1=ncsty(atmp,btmp);
Sab1=sfcy(atmp,btmp);

for i=1:50
    for j=1:50
    
    % Binarization part:  simple thresholding
    atmp=ones(Npts,1);
    btmp=ones(Npts,1);
    
    atmp(a<tha(i))=0;
    btmp(b<thb(j))=0;    
    
    t=[sum(atmp&btmp),sum(atmp&~btmp);sum(~atmp&btmp),sum(~atmp&~btmp)];
    p_y_given_x=t(1,1)/sum(t(:,1));
    p_y_given_not_x=t(1,2)/sum(t(:,2));
    PNS=p_y_given_x-p_y_given_not_x;
   
    %Necessity
    Nab(j,i)=PNS/p_y_given_x;
    if isnan(Nab(j,i)) | isinf(Nab(j,i)) | Nab(j,i)<0
        Nab(j,i)=0;
    end
    %Sufficiency
    Sab(j,i)=PNS/(1-p_y_given_not_x);
    if isnan(Sab(j,i)) | isinf(Sab(j,i)) | Sab(j,i)<0 
        Sab(j,i)=0;
    end

    
    end
end

[X,Y]=meshgrid(tha,thb);

if handles.SN2
    axes(handles.sufnec2)
    surf(X,Y,Nab);
    hold on
    line([tha(1) tha(end)],[thb(1) thb(1)],[1.2 1.2]);
    line([tha(1) tha(1)],[thb(1) thb(end)],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(n1)])
    ylabel(['threshold node ', num2str(n2)])
    zlabel('Necessity')
    axis([min(a) max(a) min(b) max(b) 0 1.2])
    colorbar
else
    axes(handles.sufnec2)
    surf(X,Y,Sab);
    hold on
    line([tha(1) tha(end)],[thb(1) thb(1)],[1.2 1.2]);
    line([tha(1) tha(1)],[thb(1) thb(end)],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(n1)])
    ylabel(['threshold node ', num2str(n2)])
    zlabel('Sufficiency')
    axis([min(a) max(a) min(b) max(b) 0 1.2])
    colorbar
end

handles.n1g2=n1;
handles.n2g2=n2;
handles.atmp2=atmp;
handles.btmp2=btmp;
handles.tha2=min(a);
handles.thb2=min(b);

handles.plot2_3d.X=X;
handles.plot2_3d.Y=Y;
handles.plot2_3d.N=Nab;
handles.plot2_3d.S=Sab;
handles.plot2_3d.tha=[tha(1) tha(end)];
handles.plot2_3d.thb=[thb(1) thb(end)];
handles.plot2_3d.a=[min(a) max(a)];
handles.plot2_3d.b=[min(b) max(b)];

set(handles.slider3,'value',0)
set(handles.slider4,'value',0)

cont1={['Nodes ',num2str(handles.n1g2),'-',num2str(handles.n2g2)],...
    ['N: ',num2str(floor(Nab1*10000)/10000),' / S: ',num2str(floor(Sab1*10000)/10000)],...
    ['Total N = ',num2str(sum(handles.plot2_3d.N(:)))],...
    ['Total S = ',num2str(sum(handles.plot2_3d.S(:)))],...
    ['theshold node ',num2str(handles.n1g2),' = ',num2str(floor(handles.tha2*10000)/10000)],...
    ['theshold node ',num2str(handles.n2g2),' = ',num2str(floor(handles.thb2*10000)/10000)]};
set(handles.lb2,'String',cont1)

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)

guidata(hObject, handles);

function slider3_Callback(hObject, eventdata, handles)

axes(handles.sufnec2)
[az,el] = view;

Npts=1000;
a=handles.X(:,handles.n1g2);
b=handles.X(:,handles.n2g2);

% scale the signal:
a=a/max(a); 
a=a-min(a);

th=get(hObject,'Value')*(max(a)+eps);

% Binarization part: simple thresholding
atmp=ones(Npts,1);
atmp(a<th) = 0;
    
axes(handles.events1_2)        
stem(atmp,'g'); 
hold on,
plot(a); 
line([0 Npts],[th th]); 
hold off
      
handles.tha2=th;
handles.atmp2=atmp;

handles.Nab2=ncsty(atmp,handles.btmp2);
handles.Sab2=sfcy(atmp,handles.btmp2);

cont1={['Nodes ',num2str(handles.n1g2),'-',num2str(handles.n2g2)],...
    ['N: ',num2str(floor(handles.Nab2*10000)/10000),' / S: ',num2str(floor(handles.Sab2*10000)/10000)],...
    ['Total N = ',num2str(sum(handles.plot2_3d.N(:)))],...
    ['Total S = ',num2str(sum(handles.plot2_3d.S(:)))],...
    ['theshold node ',num2str(handles.n1g2),' = ',num2str(floor(handles.tha2*10000)/10000)],...
    ['theshold node ',num2str(handles.n2g2),' = ',num2str(floor(handles.thb2*10000)/10000)]};
set(handles.lb2,'String',cont1)

if handles.SN2
    axes(handles.sufnec2)
    surf(handles.plot2_3d.X,handles.plot2_3d.Y,handles.plot2_3d.N);
    hold on
    line([th th],[handles.plot2_3d.thb(1) handles.plot2_3d.thb(2)],[1.2 1.2]);
    line([handles.plot2_3d.tha(1) handles.plot2_3d.tha(2)],[handles.thb2 handles.thb2],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g2)])
    ylabel(['threshold node ', num2str(handles.n2g2)])
    zlabel('Necessity')
    axis([handles.plot2_3d.a(1) handles.plot2_3d.a(2) handles.plot2_3d.b(1)...
        handles.plot2_3d.b(2) 0 1.2])
    colorbar
else
    axes(handles.sufnec2)
    surf(handles.plot2_3d.X,handles.plot2_3d.Y,handles.plot2_3d.S);
    hold on
    line([th th],[handles.plot2_3d.thb(1) handles.plot2_3d.thb(2)],[1.2 1.2]);
    line([handles.plot2_3d.tha(1) handles.plot2_3d.tha(2)],[handles.thb2 handles.thb2],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g2)])
    ylabel(['threshold node ', num2str(handles.n2g2)])
    zlabel('Sufficiency')
    axis([handles.plot2_3d.a(1) handles.plot2_3d.a(2) handles.plot2_3d.b(1)...
        handles.plot2_3d.b(2) 0 1.2])
    colorbar
end
    
guidata(hObject, handles);

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)

view([az,el]);

function slider4_Callback(hObject, eventdata, handles)

axes(handles.sufnec2)
[az,el] = view;

Npts=1000;
a=handles.X(:,handles.n1g2);
b=handles.X(:,handles.n2g2);

% scale the signal:
b=b/max(b); 
b=b-min(b);

th=get(hObject,'Value')*(max(b)+eps);

% Binarization part: simple thresholding
btmp = ones(Npts,1);
btmp(b<th) = 0;
    
axes(handles.events2_2)        
stem(btmp,'g'); 
hold on,
plot(b); 
line([0 Npts],[th th]); 
hold off
      
handles.thb2=th;
handles.btmp2=btmp;

handles.Nab2=ncsty(handles.atmp2,btmp);
handles.Sab2=sfcy(handles.atmp2,btmp);

cont1={['Nodes ',num2str(handles.n1g2),'-',num2str(handles.n2g2)],...
    ['N: ',num2str(floor(handles.Nab2*10000)/10000),' / S: ',num2str(floor(handles.Sab2*100)/100)],...
    ['Total N = ',num2str(sum(handles.plot2_3d.N(:)))],...
    ['Total S = ',num2str(sum(handles.plot2_3d.S(:)))],...
    ['theshold node ',num2str(handles.n1g2),' = ',num2str(floor(handles.tha2*10000)/10000)],...
    ['theshold node ',num2str(handles.n2g2),' = ',num2str(floor(handles.thb2*10000)/10000)]};
set(handles.lb2,'String',cont1)

if handles.SN2
    axes(handles.sufnec2)
    surf(handles.plot2_3d.X,handles.plot2_3d.Y,handles.plot2_3d.N);
    hold on
    line([handles.tha2 handles.tha2],[handles.plot2_3d.thb(1) handles.plot2_3d.thb(2)],[1.2 1.2]);
    line([handles.plot2_3d.tha(1) handles.plot2_3d.tha(2)],[th th],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g2)])
    ylabel(['threshold node ', num2str(handles.n2g2)])
    zlabel('Necessity')
    axis([handles.plot2_3d.a(1) handles.plot2_3d.a(2) handles.plot2_3d.b(1)...
        handles.plot2_3d.b(2) 0 1.2])
    colorbar
else
    axes(handles.sufnec2)
    surf(handles.plot2_3d.X,handles.plot2_3d.Y,handles.plot2_3d.S);
    hold on
    line([handles.tha2 handles.tha2],[handles.plot2_3d.thb(1) handles.plot2_3d.thb(2)],[1.2 1.2]);
    line([handles.plot2_3d.tha(1) handles.plot2_3d.tha(2)],[th th],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g2)])
    ylabel(['threshold node ', num2str(handles.n2g2)])
    zlabel('Sufficiency')
    axis([handles.plot2_3d.a(1) handles.plot2_3d.a(2) handles.plot2_3d.b(1)...
        handles.plot2_3d.b(2) 0 1.2])
    colorbar
end

guidata(hObject, handles);

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)
view([az,el]);

function SNb1_Callback(hObject, eventdata, handles)

handles.SN1=get(handles.SNb1,'Value');

axes(handles.sufnec1)
[az,el] = view;

if handles.SN1
    axes(handles.sufnec1)
    surf(handles.plot1_3d.X,handles.plot1_3d.Y,handles.plot1_3d.N);
    hold on
    line([handles.tha1 handles.tha1],[handles.plot1_3d.thb(1) handles.plot1_3d.thb(2)],[1.2 1.2]);
    line([handles.plot1_3d.tha(1) handles.plot1_3d.tha(2)],[handles.thb1 handles.thb1],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g1)])
    ylabel(['threshold node ', num2str(handles.n2g1)])
    zlabel('Necessity')
    axis([handles.plot1_3d.a(1) handles.plot1_3d.a(2) handles.plot1_3d.b(1)...
        handles.plot1_3d.b(2) 0 1.2])
    colorbar    
else
    axes(handles.sufnec1)
    surf(handles.plot1_3d.X,handles.plot1_3d.Y,handles.plot1_3d.S);
    hold on
    line([handles.tha1 handles.tha1],[handles.plot1_3d.thb(1) handles.plot1_3d.thb(2)],[1.2 1.2]);
    line([handles.plot1_3d.tha(1) handles.plot1_3d.tha(2)],[handles.thb1 handles.thb1],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g1)])
    ylabel(['threshold node ', num2str(handles.n2g1)])
    zlabel('Sufficiency')
    axis([handles.plot1_3d.a(1) handles.plot1_3d.a(2) handles.plot1_3d.b(1)...
        handles.plot1_3d.b(2) 0 1.2])
    colorbar    
end

guidata(hObject, handles);

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)
view([az,el]);

function SNb2_Callback(hObject, eventdata, handles)

handles.SN2=get(handles.SNb2,'Value');

axes(handles.sufnec2)
[az,el] = view;

if handles.SN2
    axes(handles.sufnec2)
    surf(handles.plot2_3d.X,handles.plot2_3d.Y,handles.plot2_3d.N);
    hold on
    line([handles.tha2 handles.tha2],[handles.plot2_3d.thb(1) handles.plot2_3d.thb(2)],[1.2 1.2]);
    line([handles.plot2_3d.tha(1) handles.plot2_3d.tha(2)],[handles.thb2 handles.thb2],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g2)])
    ylabel(['threshold node ', num2str(handles.n2g2)])
    zlabel('Necessity')
    axis([handles.plot2_3d.a(1) handles.plot2_3d.a(2) handles.plot2_3d.b(1)...
        handles.plot2_3d.b(2) 0 1.2])
    colorbar    
else
    axes(handles.sufnec2)
    surf(handles.plot2_3d.X,handles.plot2_3d.Y,handles.plot2_3d.S);
    hold on
    line([handles.tha2 handles.tha2],[handles.plot2_3d.thb(1) handles.plot2_3d.thb(2)],[1.2 1.2]);
    line([handles.plot2_3d.tha(1) handles.plot2_3d.tha(2)],[handles.thb2 handles.thb2],[1.2 1.2]);
    hold off
    xlabel(['threshold node ', num2str(handles.n1g2)])
    ylabel(['threshold node ', num2str(handles.n2g2)])
    zlabel('Sufficiency')
    axis([handles.plot2_3d.a(1) handles.plot2_3d.a(2) handles.plot2_3d.b(1)...
        handles.plot2_3d.b(2) 0 1.2])
    colorbar    
end

guidata(hObject, handles);

rotate_axes(handles.events1_1,handles.events1_2,...
    handles.events2_1,handles.events2_2,handles.adj_mat)
view([az,el]);

function rotate_axes(events1_1,events1_2,events2_1,events2_2,adj_mat)

setAllowAxesRotate(rotate3d,events1_1,false);
setAllowAxesRotate(rotate3d,events1_2,false);
setAllowAxesRotate(rotate3d,events2_1,false);
setAllowAxesRotate(rotate3d,events2_2,false);
setAllowAxesRotate(rotate3d,adj_mat,false);
set(rotate3d,'Enable','on');

function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider4_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function lb1_Callback(hObject, eventdata, handles)
function lb2_Callback(hObject, eventdata, handles)
function lb1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lb2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
tic;
global dt cf uston ustoff mods

%default values:
%{
mods.kcap = 2;
mods.kgl = 10;
mods.karf = 50;
mods.kel = 50;
mods.kgna = 10;
mods.kgk = 10;
%}

%function [t,v] = HHrk4(pulsei)
%tspan = [0 200];
t0 = 0;
v0 = [-64.99971605,0.05293425,0.59611197,0.31768108];


% Call the ODE solver ode15s.
%[t,v] = ode15s(@hh,tspan,v0);

%tend = 150;
%dt = 1e-2;

tend = 25;
dt = 1e-4;  % this means that each dt step is 0.1 usec.

halfdt = 0.5*dt;

%cf = 500;  % carrier frequency in kHz
uston = 5;
ustoff = 25;
% note that time is in msec and frequency is in kHz

numsteps=tend/dt;
nout=1;
t=t0;
v=v0;
tsave=t0;
vsave=v0;
for i=1:numsteps
    
    k1 = dt*hheqs(t, v, USflag);
    k2 = dt*hheqs(t + halfdt, v+0.5*k1, USflag);
    k3 = dt*hheqs(t + halfdt, v+0.5*k2, USflag);
    k4 = dt*hheqs(t + dt, v+k3, USflag);
    
    vp1 = v + (1/6)*(k1+2*k2+2*k3+k4);
    
    v = vp1;
    t = t+dt;
    
    if (rem(i,nout)==0)
        tsave=[tsave;t];
        vsave=[vsave;v];
    end
end

%figure
subplot(211)
plot(tsave,vsave(:,1))


%% Now plot the time course of the variable under scrutiny
c=1;
gna=120;
gk=36;
gl=0.3;
ena=50;
ek=-77;
el=-54.4;
iapp=0;

tt=linspace(t0,tend,numsteps);
%cf = 200;  % carrier frequency in Hz
%uston = 5;
%ustoff = 25;
uspw = heaviside(tt-uston).*sin(2*pi*cf.*(tt-uston));

switch USflag
    
    case 1
        title('US modulates capacitance')
        subplot(212)
        
        %kcap = 2;   % weighting of uspw
        
        c = c * (1 + kcap*uspw);
        c(find(c<0))=0;
        plot(tt,c)
        
        %             print -dpng cap
        %             save cap tsave vsave
        
    case 2
        title('US modulates Leak Conductance')
        %kgl = 10;
        
        gl = gl*(1 + kgl*uspw);
        gl(find(gl<0))=0;
        subplot(212)
        plot(tt,gl)
        
        %             print -dpng gl
        %             save gl tsave vsave
        
    case 3
        title('US introduces current through ARF')
        
        %karf = 50;
        
        arf = karf*uspw;
        subplot(212)
        plot(tt,arf)
        
        %             print -dpng iarf
        %             save iarf tsave vsave
        
    case 4
        title('US modulates Leak CUrrent Reversal Potential')
        %  vary leak current reversal potential
        
        %kel = 50;  % varies by +kel and -kel mV
        
        el = el + kel*uspw;
        subplot(212)
        plot(tt,el)
        
        %             print -dpng er
        %             save er tsave vsave
        
    case 5
        title('US modulates Sodium Conductance')
        
        %kgna = 10;
        
        gna = gna*(1 + kgna*uspw);
        gna(find(gna<0))=0;
        subplot(212)
        plot(tt,gna)
        
        %             print -dpng gna
        %             save gna tsave vsave
        
    case 6
        title('US modulates Potassium Conductance')
        
        %kgk = 10;
        
        gk = gk*(1 + kgk*uspw);
        gk(find(gk<0))=0;
        subplot(212)
        plot(tt,gk)
        
        %             print -dpng gk
        %             save gk tsave vsave
        
end
toc





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the ODE function as nested function,
% using the parameter pulsei.
%     function dvdt = hh(t,v)
%     c=1;
%     gna=120;
%     gk=36;
%     gl=0.3;
%     ena=50;
%     ek=-77;
%     el=-54.4;
%     ton=25;
%     toff=200;
%     dvdt(1) = (-gna*(v(2)^3)*v(3)*(v(1)-ena)-gk*(v(4)^4)*(v(1)-ek)...
%             -gl*(v(1)-el) + pulsei*heaviside(t-ton)*heaviside(toff-t))/c;
%     dvdt(2) = alpham(v(1))*(1-v(2))-betam(v(1))*v(2);
%     dvdt(3) = alphah(v(1))*(1-v(3))-betah(v(1))*v(3);
%     dvdt(4) = alphan(v(1))*(1-v(4))-betan(v(1))*v(4);
%     dvdt=dvdt';
%
%         function am=alpham(x)
%             am=-0.1*(x+40)/(exp(-(x+40)/10)-1);
%         end
%
%         function bm=betam(x)
%             bm=4*exp(-(x+65)/18);
%         end
%
%         function ah=alphah(x)
%             ah=0.07*exp(-(x+65)/20);
%         end
%
%         function bh=betah(x)
%             bh=1.0/(exp(-(x+35)/10)+1);
%         end
%
%         function an=alphan(x)
%             an=-0.01*(x+55)/(exp(-(x+55)/10)-1);
%         end
%
%         function bn=betan(x)
%             bn=0.125*exp(-(x+65)/80);
%         end
%     end


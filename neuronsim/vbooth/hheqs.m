function derivs = hheqs(t,v, USflag, uspw, mods)
% function derivs = hheqs(t,v, USflag, uspw, mods)
%
% this function computes the d/dt term of the HH equations.  
% The HH equations are being affected by an ultrasound pressure wavefor, 
% which modulates one of the constants in the model 
% It's intended to work with a 4th order Runge Kutta solver
%
% inputs: 
% t - time
% v - voltage at time t
% USflag - choice of what parameter gets modulated by the ultrasound waveform
% uspw - the ultrasound pressure waveform at time t
% mods - the amplitude of the uspw modulation on said parameter
%

% Default values of the constants
c=1;
gna=120;
gk=36;
gl=0.3;
ena=50;
ek=-77;
el=-54.4;
iapp=0;



% This is where the US pulse modulates one 
% (or more) of the parameters of the HH model 
switch USflag
    case 1
        % capacitance
        kcap = 2;   % weighting of uspw
        
        kcap = mods.kcap;   % weighting of uspw
        
        c = c * (1 + kcap*uspw);
        if c <= 0.1
            c = 0.1;
        end
        
    case 2
        % vary leak conductance - varies between gl*(1-kgl) and gl*(1+kgl), with values cut
        % off at gl=0
        
        kgl = 10;
        
        %%% AP happens between 5.5 and 
        kgl = mods.kgl;
        
        gl = gl*(1 + kgl*uspw);
        if gl < 0
            gl = 0;
        end
        
    case 3
        % add in ARF-induced current that is proportional to uspw
        karf = 50;
        
        karf = mods.karf;
        
        arf = karf*uspw;
        iapp = arf;
        
   
    case 4
        %  vary leak current reversal potential
        kel = 50;  % varies by +kel and -kel mV
        
        kel = mods.kel;
        
        el = el + kel*uspw;
    
    case 5
        % vary Na+ conductance - varies between gl*(1-kgl) and gl*(1+kgl), with values cut
        % off at gl=0
        kgna = 10;
        
        kgna = mods.kgna;
        
        gna = gna*(1 + kgna*uspw);
        if gna < 0
            gna = 0;
        end
        
    case 6
        % vary K+ conductance - varies between gk*(1-kgk) and gk*(1+kgk), with values cut
        % off at gk=0
        kgk = 10;
        
        kgk = mods.kgk;
        
        gk = gk*(1 + kgk*uspw);
        if gk < 0
            gk = 0;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables in vector v: (1) = V, (2) = m, (3) = h, (4) = n
% Right hand side of HH equations in vector dvdt:

V = v(1) ;
m = v(2) ;
h = v(3) ;
n = v(4) ;

% dvdt(1) = (-gna*(m^3)*h*(V-ena) - gk*(n^4)*(V-ek) ...
%     -gl*(V-el) + iapp) / c;
% %            -gl*(V-el) + pulsei*heaviside(t-ton)*heaviside(toff-t))/UScap(t);
% dvdt(2) = alpham(V)*(1-m) - betam(V)*m;   
% dvdt(3) = alphah(V)*(1-h) - betah(V)*h;
% dvdt(4) = alphan(V)*(1-n) - betan(V)*n;

noise = 5*square(t*(1/2.5)*2*pi, 30);

dvdt = (-gna*(m^3)*h*(V-ena) - gk*(n^4)*(V-ek) ...
    -gl*(V-el) + iapp +     noise) / c;
%            -gl*(V-el) + pulsei*heaviside(t-ton)*heaviside(toff-t))/UScap(t);
dmdt = alpham(V)*(1-m) - betam(V)*m;   
dhdt = alphah(V)*(1-h) - betah(V)*h;
dndt = alphan(V)*(1-n) - betan(V)*n;

% now package the derivatives into one vector:
derivs = [dvdt dmdt dhdt dndt];

%dvdt=dvdt';

    function am = alpham(x)
        am = -0.1*(x+40)/(exp(-(x+40)/10)-1);
    end

    function bm = betam(x)
        bm = 4*exp(-(x+65)/18);
    end

    function ah = alphah(x)
        ah = 0.07*exp(-(x+65)/20);
    end

    function bh = betah(x)
        bh = 1.0/(exp(-(x+35)/10)+1);
    end

    function an = alphan(x)
        an = -0.01*(x+55)/(exp(-(x+55)/10)-1);
    end

    function bn = betan(x)
        bn = 0.125*exp(-(x+65)/80);
    end
end

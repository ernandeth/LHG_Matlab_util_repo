function dvdt = hheqs(t,v)
c=1;
gna=120;
gk=36;
gl=0.3;
ena=50;
ek=-77;
el=-54.4;
iapp=0;
%ton=2.5;
%toff=900;
%pulsei=12;

% ultrasound pressure wave  - varies between -1 and 1 - turns on at
% t=uston with value 0
cf = 20000;  % carrier frequency in Hz
uston = 10;
ustoff = 100;
uspw = heaviside(t-uston)*sin(2*pi*cf*(t-uston)/1000);
%uspw = heaviside(t-uston)*heaviside(ustoff-t);

%     % capacitance
%     kcap = 2;   % weighting of uspw
%     UScap = c * (1 + kcap*uspw);
%     if UScap <= 0.1
%         UScap = 0.1;
%     end

% vary leak conductance - varies between gl*(1-kgl) and gl*(1+kgl), with values cut
% off at gl=0
% kgl = 10;
% USgl = gl*(1 + kgl*uspw);
% if USgl < 0
%     USgl = 0;
% end

% % add in ARF-induced current that is proportional to uspw
%     karf = 50;
%     arf = karf*uspw;

% %  vary leak current reversal potential
%     kel = 50;  % varies by +kel and -kel mV
%     USel = el + kel*uspw;

% vary Na+ conductance - varies between gl*(1-kgl) and gl*(1+kgl), with values cut
% off at gl=0
   kgna = 10;
   USgna = gna*(1 + kgna*uspw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables in vector v: (1) = V, (2) = m, (3) = h, (4) = n
% Right hand side of HH equations in vector dvdt:

V = v(1) ;
m = v(2) ;
h = v(3) ;
n = v(4) ;

dvdt(1) = (-gna*(m^3)*h*(V-ena) - gk*(n^4)*(V-ek) ...
    -USgl*(V-el) + iapp) / c;
%            -gl*(V-el) + pulsei*heaviside(t-ton)*heaviside(toff-t))/UScap(t);
dvdt(2) = alpham(V)*(1-m) - betam(V)*m;
dvdt(3) = alphah(V)*(1-h) - betah(V)*h;
dvdt(4) = alphan(V)*(1-n) - betan(V)*n;

v(1)= V ;
v(2)= m ;
v(3)= h ;
v(4)= n ;

%dvdt=dvdt';

    function am=alpham(x)
        am=-0.1*(x+40)/(exp(-(x+40)/10)-1);
    end

    function bm=betam(x)
        bm=4*exp(-(x+65)/18);
    end

    function ah=alphah(x)
        ah=0.07*exp(-(x+65)/20);
    end

    function bh=betah(x)
        bh=1.0/(exp(-(x+35)/10)+1);
    end

    function an=alphan(x)
        an=-0.01*(x+55)/(exp(-(x+55)/10)-1);
    end

    function bn=betan(x)
        bn=0.125*exp(-(x+65)/80);
    end
end

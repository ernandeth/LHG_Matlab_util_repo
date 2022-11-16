%close all
clear all
global FOV Nvox

coilconfig=2

FOV=0.24;
Nvox=128;

% this flag determines whether to scale the field so that the "target"
% always gets the same amount of E-field
targetScaling=0

% Use this flag to choose the simulations from Luis Gomez on a brain:
BrainSims=0

load bestParms_temp.mat

zview= round((0.195-0.0113)*Nvox/FOV);  % 1.13 cm below surface
yview = round(Nvox/2);
xview = round(Nvox/2);

unit_vol=(FOV/Nvox)^3;
unit_area = (FOV/Nvox)^2;

make_masks

% some hard coded parms:
current =1e3/1e-4; % this is actually dI/dt in Amps/sec.
xgap = 0;
R1 = 0.03;  % m
aperture = pi/6;
theta = 0;
zpos = 0.085;



% Calculate E field from coil alone:
%[E1 Ex1 Ey1 Ez1] =  make_fig8(current, R1, xgap, zpos, theta);
%save defaultFig8_128.mat E1 Ex1 Ey1 Ez1



%for solN=1:8
for solN=3

    %     switch coilconfig
    %         case 1
    %             %  square shield
    %             load bestParms_square.mat
    %             load defaultFig8_128.mat
    %
    %             sfactorBottom = parms(1);
    %             Rbottom = parms(2);
    %             Zbottom = parms(3);
    %
    %             Nvox=128;
    %             %[E2 Ex2 Ey2 Ez2] = make_square_shield(current, Rbottom, Zbottom);
    %             %save square_shield_128.mat E2 Ex2 Ey2 Ez2
    %             load square_shield_128.mat
    %
    %
    %
    %         case 2
    % figure 8 shield
    %load bestParms_W080
    %load bestParms_W200

    
    namestr = 'sphere';
    parms = allhst{solN}.parms(end,:);
    W = allhst{solN}.W

    sfactorBottom = parms(1) ;
    Rbottom = parms(2);
    Zbottom = parms(3);
    xgap = -(Rbottom-R1);


    if BrainSims
        namestr='head';
        %%% brain sims from Luis G
        %
        % load set3noshield.mat
        load unshieldedyz.mat
        
        Ex_1 = Emaskx; Ey_1 = Emasky; Ez_1 = Emaskz;
        E1 = Emaskmag2;

        load set4shieldedyz.mat
        E2=Emaskmag1;
        Ex2 = Emaskx; Ey2 = Emasky; Ez2 = Emaskz; 
        Ex1 = Ex_1; Ey1 = Ey_1; Ez1 = Ez_1;
        W=3;
        
        headmask=1; cylmask=1;   % FOR BRAIN - LUIS G. MAPS ONLY
        E=E2;

        h.xdim=128; h.ydim=128; h.zdim=128; h.xsize=.00214; h.ysize=.00214; h.zsize=.00107;
        unit_vol = h.zsize*(h.xsize)^2 ; % THIS IS THE BRAIN CASE
        unit_area = (h.xsize)^2;   % this is the BRAIN case

        FOV = h.zdim*h.zsize;
        zview= round((FOV-0.0128)*Nvox/FOV);  % 1.28 cm below surface
    else
        h=[];
        load defaultFig8_128.mat

        %%%%%

        %[E2 Ex2 Ey2 Ez2] =  make_fig8(current, Rbottom, xgap, Zbottom, theta);
        %eval(['save fig8_shield_128_W' num2str(W) '.mat E2 Ex2 Ey2 Ez2'])
        load (['fig8_shield_128_W' num2str(W) '.mat'])

        % compute total fields
        %
        Ex = sfactorBottom*Ex2  + Ex1;
        Ey = sfactorBottom*Ey2  + Ey1;
        Ez = sfactorBottom*Ez2  + Ez1;

        E = sqrt(Ex.^2 +Ey.^2 + Ez.^2);
    end

    %%%%%%%%%%%%%%%%%%%
    % activation threshold
    % Sthreshold = E1(xview, yview, zview);  % 2 cm below surface
    % look only at the head - both figure eight and combination
    Ehead = E .* headmask;
    E1head = E1 .* headmask;

    % look at the surface of the SPHERE !! we don'nt have a good mask for
    % the head
    Esurf = E .* headmasksurf;
    E1surf = E1 .* headmasksurf;
    
    % find the maximum Efield - ie. at surface
    E1headmax = max(E1head(:));
    Eheadmax = max(Ehead(:));
    
    if BrainSims
        E1headmax = E1head(xview,yview,128);
        Eheadmax = Ehead(xview, yview,128);
    end
    
    E1head = E1head / E1headmax;
    Ehead = Ehead / Eheadmax;

    S1threshold = 1/2;
    Sthreshold = 1/2;
        


    if targetScaling
        % this means that we force the same penetration in both cases
        % we force the same point point to have the same field in both
        % cases
        target1 = (E1head(xview:xview, yview:yview, zview:zview));
        target = (Ehead(xview:xview, yview:yview, zview:zview));
        
        sfactor = sum(target1(:))/sum(target(:));
        
        Ehead = Ehead * sfactor;
    else
        % in this case we force the scalp field to be the same in both
        % cases - i.e., = 1
        E1surf = E1surf / E1headmax;
        Esurf = Esurf / Eheadmax;

    end



    % find the supra-threshold voxels (above the maximum at the surface)

    Esupra = Ehead;
    Esupra(Ehead < Sthreshold ) = 0;

    E1supra = E1head;
    E1supra(E1head < S1threshold ) = 0;

    act_vol=length(find(Esupra(:))) * unit_vol
    act_vol1=length(find(E1supra(:))) * unit_vol

    % find the SURFACE supra-threshold voxels
    Esurfsupra = zeros(size(E));
    inds = find(Esupra);
    [x y z] = ind2sub([Nvox, Nvox, Nvox], inds);
    for n=1:length(x)
        tmp = Esupra(x(n), y(n), :);  % a vertical column through the pixel 
        zmax = max(find(tmp));
        Esurfsupra( x(n), y(n), zmax) = Esupra(x(n), y(n), zmax);
    end
    surfN = length(find(Esurfsupra));
    act_surf = surfN * unit_area
        
    E1surfsupra = zeros(size(E1));
    inds = find(E1supra);
    [x y z] = ind2sub([Nvox, Nvox, Nvox], inds);
    for n=1:length(x)
        tmp = E1supra(x(n), y(n), :);
        zmax = max(find(tmp));
        E1surfsupra( x(n), y(n), zmax) = E1supra(x(n), y(n), zmax);
    end
    surfN1 = length(find(E1surfsupra));
    act_surf1 = surfN1* unit_area

    % find the penetration
    inds = find(Esupra);
    [x y z] = ind2sub([Nvox, Nvox, Nvox], inds);
    lowestZ =  min(z)*FOV/Nvox;
    pen = 0.12 + 0.075 - min(z)*FOV/Nvox

    inds1 = find(E1supra);
    [x y z] = ind2sub([Nvox, Nvox, Nvox], inds1);
    lowestZ1 = min(z)*FOV/Nvox;
    pen1 = 0.12 + 0.075 - min(z)*FOV/Nvox

    % find the penetration
    if BrainSims
        inds = find(Esupra(xview, yview,:));
        pen = FOV - min(inds)*FOV/Nvox

        inds1 = find(E1supra(xview, yview,:));
        pen1 = FOV - min(inds1)*FOV/Nvox
    end

    penetration_change = (pen-pen1);

    %
    % switch coilconfig
    %     case 1
    %         save E_squareshield E Ex Ey Ez E1 Ex1 Ey1 Ez1
    %     case 2
    %
    %         save E_fig8shield E Ex Ey Ez E1 Ex1 Ey1 Ez1
    % end

    Escaled = Ehead;% /Eheadmax;
    E1scaled = E1head;%/E1headmax;

    LE = 20*log10(Escaled);
    LE1 = 20*log10(E1scaled);



    %LE = 20*log10(E/Eheadmax);
    %LE1 = 20*log10(E1/E1headmax);

    sfigure(5)
    subplot (311)
    plot(linspace(-FOV/2,FOV/2, 128), squeeze(LE1(:,yview, zview))), title('x axis')
    hold on
    plot(linspace(-FOV/2,FOV/2, 128),squeeze(LE(:,yview, zview)),'r'), title('x axis')
    ylabel('20log(E)')
    axis tight

    hold off

    subplot (312)
    plot(linspace(-FOV/2,FOV/2, 128),squeeze(LE1(xview,:,zview))), title('y axis')
    hold on
    plot(linspace(-FOV/2,FOV/2, 128),squeeze(LE(xview,:,zview)),'r'), title('y axis')
    hold off
    ylabel('20log(E)')
    axis tight

    subplot (313)
    plot(linspace(-FOV/2,FOV/2, 128),squeeze(LE1(xview, yview,:))), title('z axis')
    hold on
    plot(linspace(-FOV/2,FOV/2, 128),squeeze(LE(xview, yview,:)),'r'), title('z axis')
    xlabel('distance (m.)')
    legend('no shield','shield 1','Location', 'NorthWest')
    legend BOXOFF
    hold off
    ylabel('20log(E)')
    axis tight


%     LE = 20*log10(Esurfsupra);
%     LE1 = 20*log10(E1surfsupra);

%     LE = 20*log10(E/Eheadmax);
%     LE1 = 20*log10(E1/E1headmax);

    cmax = 20*log10(1);
    cmin = 20*log10(1e-2);

    sfigure(3);
    ov02(h,  LE1 ,...
        xview, yview, zview,0, cmin, cmax);
    colormap gray
    set(gcf, 'Name',['No shield']);
    eval(['print -djpeg ' namestr '_noshield'])
    eval(['print -deps ' namestr '_noshield'])

    sfigure(4);
    ov02(h,  (LE) ,...
        xview, yview, zview,0, cmin, cmax);
    set(gcf, 'Name',['Optimized with W = ' num2str(W)]);
    colormap gray
    eval(['print -djpeg ' namestr '_shield_' num2str(targetScaling) '_W=' num2str(W) ])
    eval(['print -depsc ' namestr '_shield_' num2str(targetScaling) '_W=' num2str(W) ])

    [   Esupra(xview, yview, zview)   E1supra(xview, yview, zview)]

    all_maxE(solN) = max(Esupra(:));
    all_pen(solN)= pen;
    all_act_surf(solN) = act_surf;
    all_act_vol(solN) = act_vol;
    all_W(solN) = W;
    drawnow;
end
startFields = [1e2*pen1        1e4*act_surf1       1e6*act_vol1		max(E1supra(:)) ]
optFields = [(1e2*all_pen)'    1e4*all_act_surf'   1e6*all_act_vol'	all_maxE' ]  % expressed as cm

percents = -(startFields(end,:) -optFields(end,:) )./startFields

if ~BrainSims

    OptimEndParms
    names = '       num     W    SharpChange  PenChange  shFactor  Radius     Zpos     Penet.    ActSurf  ActVol  MaxField'
    allResults =  [optparms optFields]
    allResults = [zeros(1,7) startFields; allResults]
    save OptimEndParms.txt allResults -ASCII

end
    

return
% put some test code lines here:
LE = 20*log10(Esupra);
LE1 = 20*log10(E1supra);

LE = 20*log10(Esurfsupra);
LE1 = 20*log10(E1surfsupra);

%     LE = 20*log10(E/Eheadmax);
%     LE1 = 20*log10(E1/E1headmax);

cmax = 20*log10(1);
cmin = 20*log10(1e-1);
for zview=90:120
    sfigure(3);
    ov02(h,  LE1 ,...
        xview, yview, zview,0, cmin, cmax);

    set(gcf, 'Name',['No shield']);
    %eval('print -djpeg head_noshield')


    sfigure(4);
    ov02(h,  (LE) ,...
        xview, yview, zview,0, cmin, cmax);
    set(gcf, 'Name',['Optimized with W = ' num2str(W)]);


    drawnow
end


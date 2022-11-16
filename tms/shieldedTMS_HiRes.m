        
load bestParms_temp.mat
for c=[3 4 5 8]
 	W = allhst{c}.W
	parms = allhst{c}.parms(end,:);
	R1 = 0.03;
        sfactorBottom = parms(1) ;
        Rshield = parms(2);
        Zshield = parms(3);
        FOV=0.24;
        Nvox=128;

        xgap = -(Rshield-R1);

        [E2 Ex2 Ey2 Ez2] =  make_fig8(current, Rshield, xgap, Zshield, theta);
	str = ['save fig8_shield_128_W' num2str(W) '.mat  E2 Ex2 Ey2 Ez2']
	eval(str)
	str = ['save bestParms_W' num2str(W) '.mat  parms']
	eval(str)
end


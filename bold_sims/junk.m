TI=[.1 0.4 0.8 1.2 1.5 2 2.5 3 3.5 5 7 9];
T1art=2.8
dt=0



for n=1:size(TI,2)
    FAIR_Name=sprintf('FAIR%03d',n)
	outName=sprintf('Flow%03d',n)
	ez_fair(FAIR_Name, TI(n), T1art,dt,outName);
end
 

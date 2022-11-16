% batch file to go through a bunch of subjects and see if the blind retroicor 
% increases the Zscores
% uses spmJr

subjs = [ ...]

for s=1:length(subjs)
	theshold = 2.2
	isASL=0
	
	wDir= ''
	cd(wDir)
	dummy = dir('*.img')
	datafile = dumm(1).name;
	
	load(SPM.mat)
	DesMat = SPM.xX.X
	contrast = [1 0 0 -1 0 0 0 ]
	 
	spmJr(datafile,DesMat,contrast)
	!mv Zmap_0001.img Zbefore.img
	!mv Zmap_0001.hdr Zbefore.hdr
	
	blind_retroicor(datafile, threshold, isASL)
	
	spmJr('residuals', DesMat, contrast)
	!mv Zmap_0001.img Zafter.img
	!mv Zmap_0001.hdr Zafter.hdr
	
		
	Zdiff = lightbox('Zafter')- lightbox('Zbefore');
	lightbox(Zdiff)

end

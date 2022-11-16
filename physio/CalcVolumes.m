function [volume, P] = CalcVolumes

P = spm_get(inf, '*.img','Select the images to calculate volume');

num = size(P,1);
volume = zeros(num,1);

V = spm_vol(P);

[mypath,myname,myext] = fileparts(deblank(P(1,:)));
volume_stats = fullfile(mypath,'volume_stats.txt');
fid = fopen(volume_stats,'wt');

for i=1:num,
    Y = spm_read_vols(V(i));
    volume(i) = length(find(Y>eps))*V(i).mat(1,1)*V(i).mat(2,2)*V(i).mat(3,3);
    fprintf(fid, '"%s"\t%d\n',P(i,:),volume(i));
end    

fclose(fid);
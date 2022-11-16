function vcut = vsi_Vcutoff(Mz_c, Mz_l)
%% calculate vcut based on laminar flow
max_v = 10; % max v
vstep = 0.01;
v = -max_v:vstep:max_v;
nv = length(v);
vmid = (nv+1)/2;
min_v = vstep;
v_lam = reshape(squeeze(min_v: vstep: max_v),[],1); %(-(50+step): step1: 50);

% calculate ASL signal
asl = Mz_c - Mz_l; % you can directly provide the ASL signal from Bloch simulation
asl_positive = asl(vmid+1:end);
P = [];
Mzasl = [];
for i = 1:length(v_lam)
    P = zeros(size(v_lam));
    Vmax = v_lam(i);
    nV = Vmax/vstep;
    prob = 1/nV;
    P(1:round(nV)) = prob;
    Mzv = asl_positive(:).*P;
    Mzasl(i) = sum(Mzv(:));

end
%% show v-profile vs. mean V
figure
plot(v_lam/2,Mzasl); % v/2 for converting max velocity to 
                     % mean velocity with laminar flow
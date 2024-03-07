function convert_VSI2prep_pulses(pulselen, controltype)
% convert VS labeling pulses to the new format for asl3dflex sequence
if floor(pulselen/2) ~= (pulselen/2)
    pulselen = pulselen+1;
end

rho_str=['myVSI_' num2str(pulselen) '.rho.txt'];
th_str=['myVSI_' num2str(pulselen) '.theta.txt'];
grad_str=['myVSI_' num2str(pulselen) '.grad.txt'];

rho = load(rho_str);
theta = load(th_str);
grad = load(grad_str);

rho = [rho rho];
theta = [theta theta];
if controltype==0
	grad = [grad zeros(size(grad))];
else
	grad = [grad abs(grad)];
end

mkdir(num2str(pulselen))
rho_str=['./' num2str(pulselen) '/rho.txt'];
th_str=['./' num2str(pulselen) '/theta.txt'];
grad_str=['./' num2str(pulselen) '/grad.txt'];

save(rho_str,"rho", '-ascii');
save(th_str,"theta", '-ascii');
save(grad_str,"grad", '-ascii');

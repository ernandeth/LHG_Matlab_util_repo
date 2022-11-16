tau1=17;
tau2=42;
allh_err = linspace(-2,2,3);

for h_err = allh_err  % iterate over errors in the HRF model

	H = HRF_mat(tau1+h_err, tau2, T);
	plot(H(:,1),'k')
	hold on
end

axis([0 30 -1 1.5])
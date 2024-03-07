function [loss,gradients] = modelLoss(netParms,Nlocs,Sdat,Tdat)
% define the loss function for the model: 

% compute the Loss: MSE of the calculated targets to the true targets
T_est = model(netParms,Nlocs, Sdat);
alpha = 1;
loss = alpha*mse(T_est,Tdat,DataFormat="BC");%  + 1e-1*sqrt(sum(netParms.mult10.Weights(:).^2)) + (1-alpha)*sum(abs(T_est(:)-Tdat(:)));
% compute the network weight gradients
gradients = dlgradient(loss, netParms);

plot(Tdat(:), T_est(:),'.')
ylabel('Estimate')
xlabel('Truth')
drawnow
end


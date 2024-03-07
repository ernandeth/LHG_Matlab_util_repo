function Y = model(netParms,X, source)

% This function defines the multilayer perceptron network (model)
%
Z0 = X;

Z1 = Z0 * netParms.mult1.Weights;
Z1 = norm_center(Z1);
Z1 = myReLU(Z1);  % myReLU() can't use gpuArray so I make my own reLU


Z2 =  Z1 * netParms.mult2.Weights;
Z2 = norm_center(Z2);
Z2 = myReLU(Z2) ;

Z3 =  Z2 * netParms.mult3.Weights;
Z3 = norm_center(Z3);
Z3 = myReLU(Z3) ;

Z4 =  Z3 * netParms.mult4.Weights;
Z4 = norm_center(Z4);
Z4 = myReLU(Z4);

Z5 =  Z4 * netParms.mult5.Weights;
Z5 = norm_center(Z5);
Z5 = myReLU(Z5);

Z6 =  Z5 * netParms.mult6.Weights;
Z6 = norm_center(Z6);
Z6 = myReLU(Z6);

Z7 =  Z6 * netParms.mult7.Weights;
Z7 = norm_center(Z7);
Z7 = myReLU(Z7);


Z8 =  Z7 * netParms.mult8.Weights;
Z8 = norm_center(Z8);
Z8 = myReLU(Z8);

Z9 =  Z8 * netParms.mult9.Weights;
Z9 = norm_center(Z9);
Z9 = myReLU(Z9);

Z10 = Z9 * netParms.mult10.Weights;


% Now we'll mutiply the output of this layer (grappa weights)
% by the neighbor signals to obtain the target signals.
% Those target signals will be the output of the network.
% however ....
% the MLP model and the loss function don't take complex numbers
% so we need to rebuild the complex numbers first, and then break them back
% down into real and imaginary
%
% Gweights: Npatches  x  (Nnbrs*Ncoils x Ncoils) - complex
% sources: Npatches  x  (Nnbrs*Ncoils x 1)      - complex
% tagets : Npatches  x  Ncoils                  - complex
%
% target_signal = Gweights * source_signal

%  Figure out dimensions
Nnbrs = size(X,2)/3;
Ncoils = size(source,2)/Nnbrs/2;
Npatches = size(source,1);

G_weights = complex(Z10(:,1:end/2) , Z10(:,end/2+1:end));
G_weights = reshape(G_weights, Npatches,[] , Ncoils);


% make the source signals complex numbers
source_tmp = complex(source(:,1:end/2) , source(:,end/2+1:end));
% now compute the targets using those G_weights
tmp = dlarray(zeros(Npatches, Ncoils));
for n=1:Npatches
    % dimensions: 
    % [1 x Ncoils]  =  [ 1 x Ncoils*Nnbrs ] * [Ncoils*Nnbrs  x Ncoils ]
    tmp(n,:) = source_tmp(n,:) * squeeze(G_weights(n,:,:));
end

% convert the output to a dlarray of real and imaginary numbers
Y = dlarray([real(tmp) imag(tmp)]);

clear G_weights tmp source_tmp
%}
end

function out = myReLU(in)
out=in;
out(in<0)=0;
end

function out = norm_center(in)
m=mean(in,1);
s=std(in,[],1);
out = (in-m)./s;

end
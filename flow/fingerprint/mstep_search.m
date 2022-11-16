load dictionary_big
% 
% for n=1:length(parms)
%     mtis0(n) = parms(n).mtis0;
%     f(n) = parms(n).f;
%     cbva(n) = parms(n).cbva;
%     transit(n) = parms(n).transit;
%     kfor(n) = parms(n).kfor;
%     r1tis(n) = parms(n).r1tis;
%     beta(n) = parms(n).beta;
%     Disp(n) = parms(n).Disp;
%     %Ptime(n) = parms(n).Ptime;
% end
% 

fprintf('\nSearching for T1, beta and kfor with fixed flow parms ...  ')
%step 1:  fix f, transit time, and cbva
f_0 = 40/6000;
cbva_0 = 0.02;
transit_0 = 1.2;

mask = zeros(length(parms), 1);
tmask = mask;
vmask = mask;
fmask = mask;

fmask( abs(f - f_0) < f_0/10) = 1;      

tmask( abs(transit - transit_0 ) < transit_0/10) = 1;      

vmask( abs(cbva - cbva_0) < cbva_0/10) = 1;      

mask = fmask .* vmask .* tmask;

inds = find(mask);

% % the reduced parameter space:
r_parms = parms(inds);
% 
% % the reduced dictionary
r_dict = dict(inds,:);

% Now search for the best r1tis
[best r ] = best_match(r_dict, timecourse);
best_parms = parms(best);

%% step 2: fix r1tis and flip angle 
% search for best transit time (allow all flows)
fprintf('\nSearching   fixed flow parms with fixed T1, beta and kfor with ...  ')

mask = zeros(length(parms), 1);
r1mask = mask;
betamask = mask;
kformask = mask;


r1tis_0 = best_parms.r1tis;
r1mask( abs(r1tis - r1tis_0) < r1tis_0/10 ) = 1;      

beta_0 = best_parms.beta;
betamask( abs(beta - beta_0) < beta_0/10 ) = 1;      

kfor_0 = best_parms.kfor;
kformask( abs(kfor - kfor_0) < kfor_0/10 ) = 1;      

mask = kformask .* r1mask .* betamask ;

inds = find(mask);

% the new reduced parameter space:
r_parms = parms(inds);

% the new reduced dictionary
r_dict = dict(inds,:);

% note that we get more sensitivity to flow and transit time if we
% differentiate the data and the dictionary and get rid of the first point:
r_dict = r_dict(:,2:2:end) - r_dict(:, 1:2:end);
for n=1:size(r_dict,1)
    tmp = r_dict(n,:);
    tmp = tmp - mean(tmp);
    tmp = tmp / norm(tmp);
    r_dict(n,:) = tmp;
end

tmp = timecourse(2:2:end) - timecourse(1:2:end);
tmp = tmp - mean(tmp);
tmp = tmp / norm(tmp);


% Now search for the best flow parms
[best R ] = best_match(r_dict, tmp);
best_parms = r_parms(best);

return


%% step 3: fix r1tis , transit time. Search for best flow and cbv

mask = zeros(length(parms), 1);
r1mask = mask;
tmask = mask;

r1tis_0 = best_parms.r1tis;
transit_0 = best_parms.transit;

r1mask( abs(r1tis - r1tis_0) < r1tis_0/10 ) = 1;      

tmask( abs(transit - transit_0) < transit_0/10 ) = 1;  

mask = r1mask .* tmask ;
inds = find(mask);

% the new reduced parameter space:
r_parms = parms(inds);

% the new reduced dictionary
r_dict = dict(inds,:);

% note that we get more sensitivity to flow and transit time if we
% differentiate the data and the dictionary and get rid of the first point:
%r_dict = r_dict(:,4:2:end) - r_dict(:, 3:2:end);
r_dict = r_dict(:,2:2:end) - r_dict(:, 1:2:end);
for n=1:size(r_dict,1)
    tmp = r_dict(n,:);
    tmp = tmp - mean(tmp);
    tmp = tmp / norm(tmp);
    r_dict(n,:) = tmp;
end

tmp = timecourse(2:2:end) - timecourse(1:2:end);
tmp = tmp - mean(tmp);
tmp = tmp / norm(tmp);
%tmp = timecourse;

% Now search for the best flow parms
[best R ] = best_match(r_dict, tmp);
best_parms = r_parms(best);




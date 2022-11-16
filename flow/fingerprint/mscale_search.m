function [out_parms] = mscale_search(all_timecourses, msk)

global dict1 dict2 parms1 parms2

%
fprintf('\nLoading  dictionaries from file ...  ')
load dictionary_coarse_flow
dict1 = dict;
parms1 = parms;

load dictionary_coarse_relax
dict2 = dict(:,5:end);
parms2 = parms;
%

% allocate space for output parameters
tmp = parms2(1);
tmp.mtis0 = 0;
tmp.f = 0;
tmp.cbva = 0;
tmp.transit = 0;
tmp.kfor=0;
tmp.r1tis=0;
tmp.beta=0;
tmp.Disp = 0;
out_parms = repmat(tmp, size(all_timecourses,2), 1);


Npix = size(all_timecourses,2);

for p=1:Npix
    
    if msk(p)==1
        timecourse = all_timecourses(:,p);
        timecourse2 = timecourse(5:end);
        
        fprintf('\nSearching  T1, beta and kfor with coarse flow dictionary ... %d  ', p)
        [best R] = best_match(dict1', timecourse);
       	best = best(1); % in case there is more than one max
 
        %fprintf('\n Best match in Fine relaxation dictionary: ...');
        best_parms = parms1(best);
        out_parms(p) =best_parms;
        
        %% step 2: fix r1tis and flip angle
        r1tis_0 = best_parms.r1tis;
        beta_0 = best_parms.beta;
        kfor_0 = best_parms.kfor;
        
        fprintf('\nSearching for flow parms with fixed T1, beta and kfor  ...  ')
        

	% COnstruct a mask to restrict the search        
        mtis0 = zeros(length(parms2),1);
        f = zeros(size(mtis0));
        cbva =zeros(size(mtis0));
        transit = zeros(size(mtis0));
        kfor = zeros(size(mtis0));
        r1tis =zeros(size(mtis0));
        beta = zeros(size(mtis0));
        Disp = zeros(size(mtis0));
        
        for n=1:length(parms2)
            mtis0(n) = parms2(n).mtis0;
            f(n) = parms2(n).f;
            cbva(n) = parms2(n).cbva;
            transit(n) = parms2(n).transit;
            kfor(n) = parms2(n).kfor;
            r1tis(n) = parms2(n).r1tis;
            beta(n) = parms2(n).beta;
            Disp(n) = parms2(n).Disp;
            %Ptime(n) = parms(n).Ptime;
        end
        
        mask = zeros(length(parms2), 1);
        r1mask = mask;
        betamask = mask;
        kformask = mask;
        
        % We restrict the new dictionary to those entries that have relaxation
        % parameters within 20% of the value obtained from the T1 dictionary
        range = 0.2;
        r1mask( abs(1./r1tis - 1./r1tis_0) < 1./r1tis_0 .* 1.5*range ) = 1;
        betamask( abs(beta - beta_0) < beta_0 * range ) = 1;
        kformask( abs(kfor - kfor_0) < kfor_0 * range ) = 1;
        mask = kformask .* r1mask .* betamask ;
        
        inds = find(mask);
        
        if isempty(inds),   % go to the next pixel (end this iteration though for loop)
            fprintf('\n no good values here ...')
        else
            % the new reduced parameter space:
            r_parms = parms2(inds);
            
            % the new reduced dictionary
            r_dict = dict2(inds,:);
            
            % note that we get more sensitivity to flow and transit time if we
            % differentiate the data and the dictionary and get rid of the first point:
            %{
            r_dict = r_dict(:,2:2:end) - r_dict(:, 1:2:end);
            for n=1:size(r_dict,1)
            tmp = r_dict(n,:);
            tmp = tmp - mean(tmp);
            tmp = tmp / norm(tmp);
            r_dict(n,:) = tmp;
            end
            
            % or fcourse,  we also need to subtract the timecourse
            timecourse = timecourse(2:2:end) - timecourse(1:2:end);
            %}
            
            % finally, search for the best flow params in the constricted dictionary
            tmp = timecourse2;
            tmp = tmp - mean(tmp);
            timecourse2 = tmp / norm(tmp);
            
            
            [best R] = best_match(r_dict', timecourse2);
            best = best(1); % in case there is more than one max
            out_parms(p) = r_parms(best);
            
            
            % but also stick in the fine resolution T1, kfor, and beta parms here
            out_parms(p).beta = beta_0;
            out_parms(p).r1tis = r1tis_0;
            out_parms(p).kfor = kfor_0;
            
            fprintf('\n Best Combination : ...');
            out_parms(p)
        end
    end
end


return


function [bestparms,bestind,recon_curve] = search_dictionary(obssignals, useFileDictionary, doFigs,timing_parms,r1)

if ~useFileDictionary
    fprintf('\nGenerating Dictionary ... \n')
    [dict, parms] = gen_dictionary_140220(timing_parms,r1);
    dict=abs(dict);
%     dict=abs(dict(:,2:2:end)-dict(:,1:2:end));
    fprintf('\n ... done ')
else
    load dictionary.mat
%         dict=abs(dict);

end

P = size(obssignals,1);

for m=1:P
    
    y=(obssignals(m,:));
    dict_temp=abs(dict).';
    
    %% Inner Product
    meany=mean(y);
    y = y - meany;
    dict_temp = dict_temp-repmat(mean(dict_temp,1),[25 1]);
    dict_temp_norm = sqrt(sum(dict_temp.*conj(dict_temp)));
    inner_product = conj(y)*dict_temp./dict_temp_norm./norm(y);
    [best bestind] = max(abs(inner_product));
    
    bestparms(m) = parms(bestind);
    
    
    %% Correlation
%         for n=1:size(dict,1);
%             cc = corrcoef(obssignals(m,:), dict(n,:));
%             score(n) = abs(cc(1,2));
%         end
%         best = max(score(:));
%     
%         bestind = find(score==best);
%     
%         bestparms(m) = parms(bestind);
    
    %%Normalized best-fit curve
    dictcol=squeeze(dict_temp(:,bestind));
    coef_vector =  pinv(dictcol)*y.';
    recon_curve = (dictcol*coef_vector) + meany;
end


if doFigs
    
end
return

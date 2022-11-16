function bestparms = search_dictionary(obssignals, useFileDictionary, doFigs)

if ~useFileDictionary
    fprintf('\nGenerating Dictionary ... \n')
    [dict, parms] = gen_dictionary;
    fprintf('\n ... done ')
else
    load dictionary.mat
end

%fprintf('\n Comparing the observed signal at each voxel with each dictionary entry ... \n');
score = zeros(size(dict,1), 1);
P = size(obssignals,1);

for m=1:P
    
    %     fprintf('\r percent pixels .... %03f', 100*m/P);
    parfor n=1:size(dict,1);
        cc = corrcoef(obssignals(m,:), dict(n,:));
        score(n) = abs(cc(1,2));
        %{
        figure(66)
        plot(obssignals(m,:))
        hold on,
        plot(dict(n,:));
        hold off
        %}
        
        
    end
    best = max(score(:));
    
    
    if ~isnan(best)
        
        bestind = find(score==best);
        
        if length(bestind)>1
            fprintf('\nMore than one entry fits the maximum value ... put in Nan')
            bestind=nan;
            bestparms(m).mtis0 = nan;
            bestparms(m).f = nan;
            bestparms(m).cbva = nan;
            bestparms(m).transit = nan;
            bestparms(m).kfor = nan;
            bestparms(m).r1tis = nan;
            
        else
            % found the best match in the dictionary
            bestparms(m) = parms(bestind);
            
            if doFigs
                figure(3)
                
                subplot(211)
                plot(obssignals(m,:))
                hold on
                plot(dict(bestind,:),'r')
                legend('Observed Signal', 'Dictionary Match');
                hold off
                title('Subtractions')
                xlabel('Acq. num.')
                ylabel('Signal (a.u.)')
                parms(bestind)
                
                subplot(212)
                plot(score), title('Correlation with dictionary entries')
                hold on
                plot(bestind,score(bestind),'*r')
                hold off
                xlabel('Dictionaly Entry')
                ylabel('Correlation Coefficient')
                
%                 subplot(313)
%                 hist(score,100);
                drawnow
            end
            
        end
        
    else
        bestparms(m) = nan;
        bestparms(m).mtis0 = nan;
        bestparms(m).f = nan;
        bestparms(m).cbva = nan;
        bestparms(m).transit = nan;
        bestparms(m).kfor = nan;
        bestparms(m).r1tis = nan;
        
    end
end


% save bestparmsinds bestparms score
% fprintf('\n .... done');



if doFigs
    
end
return

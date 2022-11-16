function prettyPlots2(ResultsFile)

    load(ResultsFile)

    %some extra code to make the final figures.
   
    set(0,'DefaultFigurePosition',[1 400 1200 400])
    figure

    subplot(131)
    plot(1./nlevels, 100*Spower_SNR'), hold on
    subplot(132)
    plot(1./nlevels, 100*Spower_SNR./repmat(Spower_SNR(6,:),6,1)), hold on
    subplot(133)
    plot(1./nlevels, eff_SNR'), hold on

    subplot(131)
    xlabel('SNR'), ylabel('% Probability of Detection')
    legend('no sub.','pairwise','running','surround','sinc', 'no sub GLS')
    legend('Location', 'SouthEast'), legend boxoff
    leg=legend;
    tit=title('Statistical Power');
    %axis ([0 2 0 100]), 
    fatlines,  dofontsize(10)
    hold off

    subplot(132)
    xlabel('SNR'), ylabel('Power as % of nosub+GLS')
    tit=[tit title('Relative Power')];
    %legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.', 'no sub + GLS','Location', 'SouthEast'), legend boxoff
    fatlines, dofontsize(10)
    hold off

    subplot(133)
    xlabel('SNR'), ylabel('1 / Var(c\beta)')
    tit=[tit title('Efficiency')];
    %legend('no sub.','pairwise sub.','running sub.','surround sub.','sinc sub.', 'no sub + GLS','Location', 'SouthEast'), legend boxoff
    fatlines, dofontsize(10)
    hold off

    set(tit,'fontsize',14)
    set(leg,'FontSize',9)

    set(gcf,'PaperPosition',[.5 1 [1200 400]/1200*7.5])
    str = sprintf('print -depsc2 %s', ResultsFile);
    eval(str)
return
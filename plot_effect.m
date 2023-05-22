function [] = plot_effect(effects, signifclusBoth, mask, vari, surf, folder)
% Adjusted code from SurfStat script SurfStatView.m to fit for effects &
% colourmap etc.

    [tt, maxShift] = data_to_plot(signifclusBoth, mask, effects);

    titleStr = ['Effects ' vari ' LTLE & RTLE pre vs post'];

    plot_data(tt,surf, titleStr, 'white');

    % Change colourmap
    cm=[ones(1,3)*0.7;
        213/255 62/255 79/255;
        linspace(213,244,32)'/255 linspace(62,109,32)'/255 linspace(79,67,32)'/255;
        linspace(244,253,32)'/255 linspace(109,174,32)'/255 linspace(67,97,32)'/255;
        linspace(253,254,31)'/255 linspace(174,224,31)'/255 linspace(97,139,31)'/255;
        linspace(254,255,31)'/255 linspace(224,255,31)'/255 linspace(139,191,31)'/255;
        linspace(255,230,31)'/255 linspace(255,245,31)'/255 linspace(191,152,31)'/255;
        linspace(230,171,31)'/255 linspace(245,221,31)'/255 linspace(152,164,31)'/255;
        linspace(171,102,32)'/255 linspace(221,194,32)'/255 linspace(164,165,32)'/255;
        linspace(102,50,32)'/255 linspace(194,136,32)'/255 linspace(165,189,32)'/255;
        50/255 136/255 189/255;
        ones(1,3)*0.99];
    SurfStatColormap(cm);
    
    cb = SurfStatColLim( [0 255]*maxShift );
    set(cb,'XLim',[0 255]*maxShift);
    h=get(cb,'Children');
    set(h,'XData',[0 255]*maxShift);
    set(cb,'XTick',[0 256/2 255]*maxShift);
    set(cb,'TickLabels',{'     >0.5\newlineneg. effect' '0' '     >0.5\newlinepos. effect'});
    cb.Label.String = 'f-squared';
    cb.Label.Position(2) = 3;
    set(gcf,'Position',[680 200 747 560]);
    sgtitle(['LTLE & RTLE pre vs post ' vari])


    print(gcf,[folder '/effect_' vari '_LTLEandRTLE.png'],'-dpng','-r300');
    saveas(gcf,[folder '/effect_' vari '_LTLEandRTLE.fig']);

end
function [tt, maxVal] = data_to_plot(signifclusBoth, maskTL, vals)
% Adjusted code from SurfStat script SurfStatView.m to work for effects &
% colourmap etc.

    effectShift = vals + 0.5; % So sign. ones start at 0
    maxVal = 1;

    t1=signifclusBoth.*(effectShift*253 + 1);
    t1(isnan(t1)) = 0;

    %%% TODO add fn for all negative/ all positive data
    if max(vals) <= 0 || min(vals) >= 0
        error("All effects are positive/ negative. Plotting incorrectly.")
    end

    t3=(1-signifclusBoth)*255;
    tt=(t1+t3).*maskTL;
end
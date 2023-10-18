function envwband = PEPs_EnvelopeCalc(audIn, srate, time, audionom, fbtime)
    fsteps = 8;
    freqbm = equal_xbm_bands(0,256,fsteps); % Subroutinte to divide the frequency interval into N bands of equal width along the BM. (Chimera toolbox)
    b = quad_filt_bank(freqbm, srate);      % Create a bank of fir complex filters. 
    wavfilt = zeros(size(b,2),length(audIn));
    wavenv = zeros(size(b,2),length(audIn));

    for bcnt = 1:size(b,2)
        plot(1:size(b,1),abs(b(:,bcnt)))
        hold on
    end

    for k = 1:size(b,2)
        wavfilt(k,:) = fftfilt(b(:,k), audIn);        % Apply the filter bank (fir complex filters)
        wavenv(k,:) = smooth(abs(wavfilt(k,:)),1000,'lowess');   % The envelope is the abs of complex signal as real & imaginary parts in quadrature. 

    end

    figure('Color',[1 1 1], 'Position', [0 1 1600 1000]);
    for pcnt = 1:size(wavfilt,1)
        subplot(fsteps,1,pcnt)
        plot(time(1:end),wavfilt(pcnt,:))
        title(['Cut-off frequency as BM width: ',num2str(ceil(freqbm(pcnt))),'-',num2str(ceil(freqbm(pcnt+1))), 'Hz']);
    end

    figure('Color',[1 1 1], 'Position', [0 1 1600 1000]);
    for ecnt = 1:size(wavenv,1)
        subplot(fsteps,1,ecnt)
        plot(time(1:end), wavenv(ecnt,:),'k')   %apply smoothing for visualisation
        title(['Narrow-band Envelope: Cut-off frequency as BM width: ',num2str(ceil(freqbm(ecnt))),'-',num2str(ceil(freqbm(ecnt+1))), 'Hz']);
    end
    
    %% To get wide band envelope, sum the narrow-band envelopes.

    envwband = sum(wavenv,1);

    figure('Color',[1 1 1],'Position',[0 1 1800 300])
    plot(time(1:end), envwband,'r')
    title('Wide band envelope (0Hz - 256Hz)')
    set(gca,'YGrid','on', 'XGrid','on')
    
%     fbtime = fbtime.*.01
  
    figure
    plot(time, envwband)
    hold on
    stem(time,fbtime,'r')
    scrollplot(gca,'WindowSizeX',5, 'MinX',0)
    x = gca;
    x.Title.String = audionom;

end
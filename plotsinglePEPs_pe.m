function plotsinglePEPs_pe (hdl,~)

D = get(hdl,'UserDat');
ecnt = D{1,1};
Cond1 = D{1,2}';
Cond2 = D{1,3}';
colrz = D{1,4};
time = D{1,5};
legendnoms = D{1,6};
eois = D{1,7}; 
topermut = D{1,8};
fbstat = D{1,9};

if size(eois,1)==1
    eois =eois';
end

assignin('base','Cond1',Cond1)
assignin('base','Cond2',Cond2)

f1 = figure;
set(f1,'Color',[1 1 1]);

axs = zeros(1,size(topermut,2));
linz = zeros(1,size(topermut,2));
           
[axs1,ploth1] = plot_ConfInt2(Cond1',time,nanmean(Cond1,1)',colrz(1,:));   %Call of plotConfInt2() function to calculate and plot the 95% CI (sem)
hold on
axs(1,1) = axs1;
linz(1,1) = ploth1;

[axs2,ploth2]= plot_ConfInt2(Cond2',time,nanmean(Cond2,1)',colrz(2,:));
hold on
axs(1,2) = axs2;
linz(1,2) = ploth2;

% Add duration annotations
ha1 = annotation('arrow');
ha1.Parent = f1.CurrentAxes;  % associate annotation with current axes
% now you can use data units
ha1.X = [0 fbstat(1,1)*1000];
ha1.Y = [-1.5 -1.5];
ha1.Color = 'blue';
ha1.LineWidth = 2;

ha2 = annotation('arrow');
ha2.Parent = f1.CurrentAxes;  % associate annotation with current axes
% now you can use data units
ha2.X = [0 fbstat(1,2)*1000];
ha2.Y = [-1.3 -1.3];
ha2.Color = 'red';
ha2.LineWidth = 2;

% Add the standard deviation information of dur1
std1_upper = fbstat(1,1)+fbstat(2,1);
if (std1_upper*1000)>time(end)
    std1_upper = time(end)/1000;
end
std1_lower = fbstat(1,1)-fbstat(2,1);
ha1b = annotation('doublearrow');
ha1b.Parent = f1.CurrentAxes;
ha1b.X = [std1_lower*1000 std1_upper*1000];
ha1b.Y = [-1.5 -1.5];
ha1b.Color = [0.3 0.3 0.3];
ha1b.LineStyle = '--';


% Add the standard deviation information of dur2
std2_upper = fbstat(1,2)+fbstat(2,2);
if (std2_upper*1000)>time(end)
    std2_upper = time(end)/1000;
end
std2_lower = fbstat(1,2)-fbstat(2,2);
ha2b = annotation('doublearrow');
ha2b.Parent = f1.CurrentAxes;
ha2b.X = [std2_lower*1000 std2_upper*1000];
ha2b.Y = [-1.3 -1.3];
ha2b.Color = [0.3 0.3 0.3];
ha2b.LineStyle = '--';

set(gca,'YDir','reverse','XAxisLocation','origin','YAxisLocation','origin','Box','off',...
        'YGrid','off','XGrid','off') % Set the current axis properties
Titre = title(eois{ecnt,1});
Titre.FontSize = 16;
xlabel(gca,'Time (ms)','FontSize',12)
ylabel(gca,'Potential (\muV)','FontSize', 12);

[tvals_curr,pvalue_corr,permout] = plot_perm(topermut,time,8,[time(1) time(2)],'yes');    
time_sig = time(~isnan(permout));
if ~isempty(time_sig)
    lowlims = time_sig(1);
    hilims = time_sig(end);
    ylimits = get(gca,'YLim');
    hold(gca,'on');
    shadcols={[1 1 0.2],[0.7 1 1],[1 0.8 0.8],[0.8 0.8 1],[0.8 1 0.8],[10.8 1],[1 1 0.7],[0.7 1 1],[1 0.8 0.8],[0.8 0.8 1],[0.8 1 0.8],[10.8 1]};
    for icnt=1:length(lowlims)
        yintval =ones(1,length(lowlims(icnt):hilims(icnt)))*ylimits(1);
        ha =fill([lowlims(icnt) lowlims(icnt):hilims(icnt) hilims(icnt)],[ylimits(2) yintval ylimits(2)],shadcols{1,icnt},'FaceAlpha',0.3,'EdgeColor','none');
        assignin('base','ha',ha);
        hold on
    end
end


legend(linz,legendnoms,'FontSize', 14);


end 

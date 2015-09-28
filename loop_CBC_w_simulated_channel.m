
list=transpose(1:10:200);
ResSim=[];
figure('Position',[100 300 400 400])

for n=1:length(list);
    
    [Ca, RHO, coloc]=CBC_w_simulated_channel(0,list(n));
    
    ResSim(n,1)=list(n);
    ResSim(n,2)=mean(Ca(~isnan(Ca)));
    ResSim(n,3)=var(Ca(~isnan(Ca)));
    ResSim(n,4)=median(Ca(~isnan(Ca)));
    ResSim(n,5)=coloc;
    
    ksdensity(Ca); hold on;
    
    
    clear Ca Rho coloc
end

figure('Position',[500 300 700 700])

subplot(2,2,1)
scatter(ResSim(:,1),ResSim(:,2));
% title('Mean C_A');
xlabel('loc density ratio (Ch1/Ch2)','FontSize',10);
ylabel('Mean C_A','FontSize',10);
box on;

subplot(2,2,2)
scatter(ResSim(:,1),ResSim(:,3));
% title('Var C_A');
xlabel('loc density ratio (Ch1/Ch2)','FontSize',10);
ylabel('variance C_A','FontSize',10);
box on;

subplot(2,2,3)
scatter(ResSim(:,1),ResSim(:,4));
% title('Median C_A');
xlabel('loc density ratio (Ch1/Ch2)','FontSize',10);
ylabel('Median C_A','FontSize',10);
box on;

subplot(2,2,4)
scatter(ResSim(:,1),ResSim(:,5));
% title('% colocalized');
xlabel('loc density ratio (Ch1/Ch2)','FontSize',10);
ylabel('% colocalized (C_A>0)','FontSize',10);
box on;
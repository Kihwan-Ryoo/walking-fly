clear all; clc;
files = dir('data\*.EDR');
    % 210811_HCS_001 ~ 004  xxxxx
    % 210811_YW_001 ~ 002   &   210812_YW_003 ~ 007     xxxxx
    % 까지 총 11개는 xvel, yvel 값이 2배 크게 측정되었다.     xxxxx
    halflist = ['HCS_Verror_001.EDR','HCS_Verror_002.EDR','HCS_Verror_003.EDR','HCS_Verror_004.EDR'];
    
    f1 = figure('Name','HCS_Verror_heading','NumberTitle','off');   
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f11 = figure('Name','n_HCS_Verror_heading','NumberTitle','off');   
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f2 = figure('Name','HCS_heading','NumberTitle','off');         
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f22 = figure('Name','n_HCS_heading','NumberTitle','off');  
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f3 = figure('Name','JH_SS54295_heading','NumberTitle','off');   
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f33 = figure('Name','n_JH_SS54295_heading','NumberTitle','off');  
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f4 = figure('Name','KH_SS02718_heading','NumberTitle','off');   
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f44 = figure('Name','n_KH_SS02718_heading','NumberTitle','off');  
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f5 = figure('Name','SY_SS00096_heading','NumberTitle','off');   
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f55 = figure('Name','n_SY_SS00096_heading','NumberTitle','off');  
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f6 = figure('Name','YW_SS00078_heading','NumberTitle','off');   
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    f66 = figure('Name','n_YW_SS00078_heading','NumberTitle','off');  
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    
    % how many flies ~?
    n_HCS_Verror = 0;
    n_HCS = 0;
    n_JH_SS54295 = 0;
    n_KH_SS02718 = 0;
    n_SY_SS00096 = 0;
    n_YW_SS00078 = 0;
    for n_file = 1:length(files)
        if contains(files(n_file).name,'Verror')
            n_HCS_Verror = n_HCS_Verror + 1;
        elseif contains(files(n_file).name,'HCS')
            n_HCS = n_HCS + 1;
        elseif contains(files(n_file).name,'JH')
            n_JH_SS54295 = n_JH_SS54295 + 1;
        elseif contains(files(n_file).name,'KH')
            n_KH_SS02718 = n_KH_SS02718 + 1;
        elseif contains(files(n_file).name,'SY')
            n_SY_SS00096 = n_SY_SS00096 + 1;
        elseif contains(files(n_file).name,'YW')
            n_YW_SS00078 = n_YW_SS00078 + 1;
        end
    end

    m1 = 0; m2 = 0; m3 = 0; m4 = 0; m5 = 0; m6 = 0;
    M1 = zeros(n_HCS_Verror,round(10^3/(2*pi)) + 1); M2 = zeros(n_HCS,round(10^3/(2*pi)) + 1); M3 = zeros(n_JH_SS54295,round(10^3/(2*pi)) + 1);
    M4 = zeros(n_KH_SS02718,round(10^3/(2*pi)) + 1); M5 = zeros(n_SY_SS00096,round(10^3/(2*pi)) + 1); M6 = zeros(n_YW_SS00078,round(10^3/(2*pi)) + 1);
    W1 = zeros(n_HCS_Verror,1); W2 = zeros(n_HCS,1); W3 = zeros(n_JH_SS54295,1);
    W4 = zeros(n_KH_SS02718,1); W5 = zeros(n_SY_SS00096,1); W6 = zeros(n_YW_SS00078,1);
    
for file_idx = 1:length(files)
    [data,h] = import_edr(strcat('data\',files(file_idx).name));
    t = data(:,1);
    xpos = data(:,2); heading = 10 - xpos;
    xvel = data(:,4); yvel = data(:,5);
    if contains(halflist,files(file_idx).name) == 1
        xvel = xvel/2; yvel = yvel/2;
    end
    pat_id = data(:,6);
    patid = round(pat_id*10);
    
    %pattern이 켜져있을 때만의 정보를 가져온다.
    t = t(patid>0); t = t - t(1);
    xpos = xpos(patid>0);   Xpos = xpos; xpos_raw = xpos;
    xvel = xvel(patid>0);
    yvel = yvel(patid>0);
    heading = heading(patid>0);
    
    % cut slightly data for good looking
    if (size(t,1)>1200000)
        t = t(1:1200000);
        xpos = xpos(1:1200000); Xpos = Xpos(1:1200000); xpos_raw = xpos_raw(1:1200000);
        xvel = xvel(1:1200000); yvel = yvel(1:1200000);
        heading = heading(1:1200000);
    end
    
    %get information only for walking movements (remove stop points)
    windowSize = 100*2;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    Xvel = filter(b,a,xvel);
    Yvel = filter(b,a,yvel);
    std = sqrt(Xvel.^2 + Yvel.^2);
    num = 0.11;
    T = t.*(std>num);
    xpos = xpos(std > num);
    BAR = T > 0;
    
    pct = (sum(BAR)/size(BAR,1))*100; pct = round(pct,2);
    
    % Let's get info only for continuous walking
    % sampling interval of WinEDR which we used for recording is 1ms
        % then sampling rate of WinEDR is 1kHz
    WS = 1000*5;
    nwindow = 1 + size(BAR,1) - WS;
    contiSTD = 0.8;
    BARstd = zeros(size(BAR));
    for i = 1:nwindow
        if sum(BAR(i:i+WS-1)) > WS * contiSTD
            BARstd(i) = 1;
        else
            BARstd(i) = 0;
        end
    end
    
    BARspec = zeros(size(BAR));
    for j = 1:nwindow
        if BARstd(j) == 1
            BARspec(j:j+WS-1) = 1;
        end
    end
    
    pctspec = (sum(BARspec)/size(BARspec,1))*100; pctspec = round(pctspec,2);
    
    Xpos = Xpos(BARspec>0);
    
    % Change range of xpos & Xpos from 0 ~ 10 to 0 ~ 2*pi
    xpos = xpos*2*pi/10;
    Xpos = Xpos*2*pi/10;
    
    w = ones(size(Xpos));
    n = sum(w);
    r = sum(w.*exp(1i*Xpos));
    
    [s,s0] = circ_std(Xpos,w,(pi/10));
    tt = 1.96 * s0 / sqrt(n);
    
    mu = angle(r); %%% *****
%     if mu < 0
%         mu = 2*pi + mu;
%     end
    MU = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
    rho = MU/s0/(5*(10^3));
    
    ul = mu + tt;
    ll = mu - tt;
    UL = linspace(ll,ul,100);
    RHO = ones(100,1)*max(rho);
    
    % figure를 그리고 png로 저장한다.
    figure('Name',files(file_idx).name,'NumberTitle','off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'Color','w');
    
    subplot(6,1,1);
    plot(xpos_raw,'r');
    hold on; b_raw = bar(BARspec*10,'c');
    title('heading','Interpreter','none');
    
    subplot(6,1,2);
    plot(xvel); hold on; plot(Xvel,'r');
    xlabel('t (au)');ylabel('V (au)');
    set(gca,'Box','off','TickDir','out');
    title('x velocity','Interpreter','none');
    
    subplot(6,1,3);
    plot(yvel); hold on; plot(Yvel,'r');
    xlabel('t (au)');ylabel('V (au)');
    set(gca,'Box','off','TickDir','out');
    title('y velocity','Interpreter','none');
    
    subplot(6,1,4);
    bar(BAR);
    title(pct);
    
    subplot(615);
    bar(BARspec);
    title(pctspec);
    
    subplot(6,1,6);
    plot(std);
    hold on;
    plot(ones(size(std))*num,'m');
    
    cd result\;
    print(files(file_idx).name(1:(end-4)),'-dpng');
    cd ..\;
    
    BE = linspace(0,(2*pi),21);
    if contains(files(file_idx).name,'Verror')
        figure(f1); nexttile;
        polarhistogram(xpos,BE,'FaceAlpha',.1,'Normalization','probability');
        set(gca,'ThetaZeroLocation','top');
        title(strcat(files(file_idx).name(1:(end-4)), ' (', string(pctspec), '%) '),'Interpreter','none');
        hold on;
        polarhistogram(Xpos,BE,'FaceColor','g','Normalization','probability');
        hold on;
        polarplot(MU,rho,'LineWidth',2,'Color','r');
        hold on;
        polarplot(UL,RHO,'LineWidth',2,'Color','r');
        
        figure(f11);
        m1 = m1 + 1;
        W1(m1) = pctspec;
        M1(m1) = mu;
        M1(m1,:) = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
        R1 = M1(m1,:)/s0/(5*(10^3));
        polarplot(M1(m1,:),R1,'LineWidth',2,'Color','k'); hold on;
        set(gca,'ThetaZeroLocation','top');
        if M1(n_HCS_Verror,:) ~= 0
            n = n_HCS_Verror;
            r = sum(W1.*exp(1i*M1(:,1)));
            
            [s,s0] = circ_std(M1(:,1),W1);
            tt = 1.96 * s0 / sqrt(n);
            
            mu = angle(r);
%             if mu < 0
%                 mu = 2*pi + mu;
%             end
            MU = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
            rho = MU/s0/(5*(10^3));
    
            ul = mu + tt;
            ll = mu - tt;
            UL = linspace(ll,ul,100);
            RHO = ones(100,1)*max(rho);
            
            polarplot(MU,rho,'LineWidth',2,'Color','r');
            hold on; 
            polarplot(UL,RHO,'LineWidth',2,'Color','r');
            rlim([0 0.41]);
        end
    elseif contains(files(file_idx).name,'HCS')
        figure(f2); nexttile;
        polarhistogram(xpos,BE,'FaceAlpha',.1,'Normalization','probability');
        set(gca,'ThetaZeroLocation','top');     % set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
        title(strcat(files(file_idx).name(1:(end-4)), ' (', string(pctspec), '%) '),'Interpreter','none');
        hold on;
        polarhistogram(Xpos,BE,'FaceColor','g','Normalization','probability');
        hold on;
        polarplot(MU,rho,'LineWidth',2,'Color','r');
        hold on;
        polarplot(UL,RHO,'LineWidth',2,'Color','r');
        
        figure(f22);
        m2 = m2 + 1;
        W2(m2) = pctspec;
        M2(m2) = mu;
        M2(m2,:) = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
        R2 = M2(m2,:)/s0/(5*(10^3));
        polarplot(M2(m2,:),R2,'LineWidth',2,'Color','k'); hold on;
        set(gca,'ThetaZeroLocation','top');
        if M2(n_HCS,:) ~= 0
            n = n_HCS;
            r = sum(W2.*exp(1i*M2(:,1)));
            
            [s,s0] = circ_std(M2(:,1),W2);
            tt = 1.96 * s0 / sqrt(n);
            
            mu = angle(r);
%             if mu < 0
%                 mu = 2*pi + mu;
%             end
            MU = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
            rho = MU/s0/(5*(10^3));
    
            ul = mu + tt;
            ll = mu - tt;
            UL = linspace(ll,ul,100);
            RHO = ones(100,1)*max(rho);
            
            polarplot(MU,rho,'LineWidth',2,'Color','r');
            hold on; 
            polarplot(UL,RHO,'LineWidth',2,'Color','r');
            rlim([0 0.41]);
        end
    elseif contains(files(file_idx).name,'JH')
        figure(f3); nexttile;
        polarhistogram(xpos,BE,'FaceAlpha',.1,'Normalization','probability');
        set(gca,'ThetaZeroLocation','top');
        title(strcat(files(file_idx).name(1:(end-4)), ' (', string(pctspec), '%) '),'Interpreter','none');
        hold on;
        polarhistogram(Xpos,BE,'FaceColor','g','Normalization','probability');
        hold on;
        polarplot(MU,rho,'LineWidth',2,'Color','r');
        hold on;
        polarplot(UL,RHO,'LineWidth',2,'Color','r');
        
        figure(f33); 
        m3 = m3 + 1;
        W3(m3) = pctspec;
        M3(m3) = mu;
        M3(m3,:) = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
        R3 = M3(m3,:)/s0/(5*(10^3));
        polarplot(M3(m3,:),R3,'LineWidth',2,'Color','k'); hold on;
        set(gca,'ThetaZeroLocation','top');
        if M3(n_JH_SS54295,:) ~= 0
            n = n_JH_SS54295;
            r = sum(W3.*exp(1i*M3(:,1)));
            
            [s,s0] = circ_std(M3(:,1),W3);
            tt = 1.96 * s0 / sqrt(n);
            
            mu = angle(r);
%             if mu < 0
%                 mu = 2*pi + mu;
%             end
            MU = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
            rho = MU/s0/(5*(10^3));
    
            ul = mu + tt;
            ll = mu - tt;
            UL = linspace(ll,ul,100);
            RHO = ones(100,1)*max(rho);
            
            polarplot(MU,rho,'LineWidth',2,'Color','r');
            hold on; 
            polarplot(UL,RHO,'LineWidth',2,'Color','r');
            rlim([0 0.41]);
        end
    elseif contains(files(file_idx).name,'KH')
        figure(f4); nexttile;
        polarhistogram(xpos,BE,'FaceAlpha',.1,'Normalization','probability');
        set(gca,'ThetaZeroLocation','top');
        title(strcat(files(file_idx).name(1:(end-4)), ' (', string(pctspec), '%) '),'Interpreter','none');
        hold on;
        polarhistogram(Xpos,BE,'FaceColor','g','Normalization','probability');
        hold on;
        polarplot(MU,rho,'LineWidth',2,'Color','r');
        hold on;
        polarplot(UL,RHO,'LineWidth',2,'Color','r');
        
        figure(f44);
        m4 = m4 + 1;
        W4(m4) = pctspec;
        M4(m4) = mu;
        M4(m4,:) = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
        R4 = M4(m4,:)/s0/(5*(10^3));
        polarplot(M4(m4,:),R4,'LineWidth',2,'Color','k'); hold on;
        set(gca,'ThetaZeroLocation','top');
        if M4(n_KH_SS02718,:) ~= 0
            n = n_KH_SS02718;
            r = sum(W4.*exp(1i*M4(:,1)));
            
            [s,s0] = circ_std(M4(:,1),W4);
            tt = 1.96 * s0 / sqrt(n);
            
            mu = angle(r);
%             if mu < 0
%                 mu = 2*pi + mu;
%             end
            MU = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
            rho = MU/s0/(5*(10^3));
    
            ul = mu + tt;
            ll = mu - tt;
            UL = linspace(ll,ul,100);
            RHO = ones(100,1)*max(rho);
            
            polarplot(MU,rho,'LineWidth',2,'Color','r');
            hold on; 
            polarplot(UL,RHO,'LineWidth',2,'Color','r');
            rlim([0 0.41]);
        end
    elseif contains(files(file_idx).name,'SY')
        figure(f5); nexttile;
        polarhistogram(xpos,BE,'FaceAlpha',.1,'Normalization','probability');
        set(gca,'ThetaZeroLocation','top');
        title(strcat(files(file_idx).name(1:(end-4)), ' (', string(pctspec), '%) '),'Interpreter','none');
        hold on;
        polarhistogram(Xpos,BE,'FaceColor','g','Normalization','probability');
        hold on;
        polarplot(MU,rho,'LineWidth',2,'Color','r');
        hold on;
        polarplot(UL,RHO,'LineWidth',2,'Color','r');
        
        figure(f55);
        m5 = m5 + 1;
        W5(m5) = pctspec;
        M5(m5) = mu;
        M5(m5,:) = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
        R5 = M5(m5,:)/s0/(5*(10^3));
        polarplot(M5(m5,:),R5,'LineWidth',2,'Color','k'); hold on;
        set(gca,'ThetaZeroLocation','top');
        if M5(n_SY_SS00096,:) ~= 0
            n = n_SY_SS00096;
            r = sum(W5.*exp(1i*M5(:,1)));
            
            [s,s0] = circ_std(M5(:,1),W5);
            tt = 1.96 * s0 / sqrt(n);
            
            mu = angle(r);
%             if mu < 0
%                 mu = 2*pi + mu;
%             end
            MU = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
            rho = MU/s0/(5*(10^3));
    
            ul = mu + tt;
            ll = mu - tt;
            UL = linspace(ll,ul,100);
            RHO = ones(100,1)*max(rho);
            
            polarplot(MU,rho,'LineWidth',2,'Color','r');
            hold on; 
            polarplot(UL,RHO,'LineWidth',2,'Color','r');
            rlim([0 0.41]);
        end
    elseif contains(files(file_idx).name,'YW')
        figure(f6); nexttile;
        polarhistogram(xpos,BE,'FaceAlpha',.1,'Normalization','probability');
        set(gca,'ThetaZeroLocation','top');
        title(strcat(files(file_idx).name(1:(end-4)), ' (', string(pctspec), '%) '),'Interpreter','none');
        hold on;
        polarhistogram(Xpos,BE,'FaceColor','g','Normalization','probability');
        hold on;
        polarplot(MU,rho,'LineWidth',2,'Color','r');
        hold on;
        polarplot(UL,RHO,'LineWidth',2,'Color','r');
        
        figure(f66); 
        m6 = m6 + 1;
        W6(m6) = pctspec;
        M6(m6) = mu;
        M6(m6,:) = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
        R6 = M6(m6,:)/s0/(5*(10^3));
        polarplot(M6(m6,:),R6,'LineWidth',2,'Color','k'); hold on;
        set(gca,'ThetaZeroLocation','top');
        if M6(n_YW_SS00078,:) ~= 0
            n = n_YW_SS00078;
            r = sum(W6.*exp(1i*M6(:,1)));
            
            [s,s0] = circ_std(M6(:,1),W6);
            tt = 1.96 * s0 / sqrt(n);
            
            mu = angle(r);
%             if mu < 0
%                 mu = 2*pi + mu;
%             end
            MU = linspace(mu,mu + round(10^3/(2*pi))*(2*pi),round(10^3/(2*pi)) + 1);
            rho = MU/s0/(5*(10^3));
    
            ul = mu + tt;
            ll = mu - tt;
            UL = linspace(ll,ul,100);
            RHO = ones(100,1)*max(rho);
            
            polarplot(MU,rho,'LineWidth',2,'Color','r');
            hold on;
            polarplot(UL,RHO,'LineWidth',2,'Color','r');
            rlim([0 0.41]);
        end
    end
end
cd result\;
print(f1,f1.Name,'-dpng'); print(f2,f2.Name,'-dpng'); print(f3,f3.Name,'-dpng'); 
print(f4,f4.Name,'-dpng'); print(f5,f5.Name,'-dpng'); print(f6,f6.Name,'-dpng');
print(f11,f11.Name,'-dpng'); print(f22,f22.Name,'-dpng'); print(f33,f33.Name,'-dpng'); 
print(f44,f44.Name,'-dpng'); print(f55,f55.Name,'-dpng'); print(f66,f66.Name,'-dpng');
cd ..\;
close all;

function  Plot_r2(data,dh_vals)
% for a single value of J this plot r varying with L and t and
% overlays the theoretical values for localised and delocalised systems

% get vals
t_vals=unique(data.t);
%t_vals=0.1:0.2:1.5;
if nargin==1
    dh_vals=unique(data.dh);
end
h_vals=unique(data.h);
L_vals=unique(data.L);


% get a plot
fig=figure(1);
clf;
hold on;

% legend data
leg={};

% magic numbers
r_goe=4-2*sqrt(3);
r_poi=2*log(2)-1;
% markers
mkr={'^','v','>','<','o','+','*','x','s','d','p','h'};
% colours
clr={};
for ih=1:length(dh_vals)
    clr{end+1}=hsv2rgb([ih/length(dh_vals),1,0.85]);
end

% add the theoretical values
plot([0,max(t_vals)],[r_goe,r_goe],'k');
plot([0,max(t_vals)],[r_poi,r_poi],'k');
xlim([0,max(t_vals)]);
ylim([0,0.6]);
leg{1}='GOE (0.53...)';
leg{2}='Poisson (0.39...)';

% for each value of L
for iL=1:length(L_vals)
    % get L
    L=L_vals(iL);
    % get indices
    IL=(data.L==L);
    % get marker
    imkr=mkr{iL};
    
    % for each value of t
    for ih=1:length(dh_vals)
        % get dh
        dh=dh_vals(ih);
        % get indices
        Idh=abs(data.dh-dh)<1E-13;
        % get marker
        iclr=clr{ih};
        
        % get index for relevant data
        I=Idh&IL;
        % get the data
        it=data.t(I);
        ir=data.r(I);
        % sort it
        [~,J]=sort(it);
        %disp([t,length(idh)])
        % plot the data
        plot(it(J),ir(J),[':',imkr],...
            'Color',iclr);
        % legend data
        leg{end+1}=['L = ',num2str(L),' dh = ',num2str(dh)];
        
        
    end
end

% add legend
legend(leg,'Location','eastoutside');
% 
set(fig, 'Position', [20 20 1200 600])

xlabel('t');
ylabel('\langler\rangle');
title('Plot of t versus level spacing ratio \langler\rangle');

end


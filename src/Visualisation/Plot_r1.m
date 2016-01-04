function  Plot_r1(data,t_vals)
% for a single value of J this plot r varying with L and dh and
% overlays the theoretical values for localised and delocalised systems

% get vals
dh_vals=unique(data.dh);
%t_vals=0.1:0.2:1.5;
if nargin==1
    t_vals=unique(data.t);
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
for it=1:length(t_vals)
    clr{end+1}=hsv2rgb([it/length(t_vals),1,0.85]);
end

% add the theoretical values
plot([0,max(dh_vals)],[r_goe,r_goe],'k');
plot([0,max(dh_vals)],[r_poi,r_poi],'k');
xlim([0,max(dh_vals)]);
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
    for it=1:length(t_vals)
        % get t
        t=t_vals(it);
        % get indices
        It=abs(data.t-t)<1E-13;
        % get marker
        iclr=clr{it};
        
        % get index for relevant data
        I=It&IL;
        % get the data
        idh=data.dh(I);
        ir=data.r(I);
        % sort it
        [~,J]=sort(idh);
        %disp([t,length(idh)])
        % plot the data
        plot(idh(J),ir(J),[':',imkr],...
            'Color',iclr);
        % legend data
        leg{end+1}=['L = ',num2str(L),' t = ',num2str(t)];
        
        
    end
end

% add legend
legend(leg,'Location','eastoutside');
% 
set(fig, 'Position', [20 20 1200 600]);

xlabel('dh');
ylabel('\langler\rangle');
title('Plot of \deltah versus level spacing ratio \langler\rangle');

end


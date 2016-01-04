function  Plot_x1(data,t_vals)
% for a single value of J this plot 1/xi varying with L and dh

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

% markers
mkr={'^','v','>','<','o','+','*','x','s','d','p','h'};
% colours
clr={};
for it=1:length(t_vals)
    clr{end+1}=hsv2rgb([it/length(t_vals),1,0.85]);
end

% add the theoretical values
xlim([0,max(dh_vals)]);
ylim([0,ceil(max(data.b))]);

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
        ib=data.b(I);
        % sort it
        [~,J]=sort(idh);
        
        % plot the data
        plot(idh(J),ib(J),[':',imkr],...
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
ylabel('1/\xi');
title('Plot of \deltah versus 1/\xi');

end


function  Plot_x2(data,dh_vals)
% for a single value of J this plot 1/xi varying with L and t

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

% markers
mkr={'^','v','>','<','o','+','*','x','s','d','p','h'};
% colours
clr={};
for ih=1:length(dh_vals)
    clr{end+1}=hsv2rgb([ih/length(dh_vals),1,0.85]);
end

% add the theoretical values
xlim([0,max(t_vals)]);
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
    for ih=1:length(dh_vals)
        % get t
        dh=dh_vals(ih);
        % get indices
        Idh=abs(data.dh-dh)<1E-13;
        % get marker
        iclr=clr{ih};
        
        % get index for relevant data
        I=Idh&IL;
        % get the data
        it=data.t(I);
        ib=data.b(I);
        % sort it
        [~,J]=sort(it);
        
        % plot the data
        plot(it(J),ib(J),[':',imkr],...
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
ylabel('1/\xi');
title('Plot of t versus 1/\xi');


end


run_name='run_1';
pool_size=1;


params={};
J_vals=0;
h_vals=1;
t_vals=0.1:0.1:1.5;
dh_vals=linspace(0.1,1.5,15);
L_vals=[6,8,10];


for J=J_vals
    for h=h_vals
        for t=t_vals
            for dh=dh_vals
                for L=L_vals
                    p.J=J;
                    p.h=h;
                    p.t=t;
                    p.dh=dh;
                    p.L=L;
                    if L==6
                        N=200;
                    elseif L==8
                        N=100;
                    elseif L==10;
                        N=50;
                    end
                    params{end+1}=p;
                end
            end
        end
    end
end
length(params)
include_all;
data_gen_east_model_stats(params,run_name,pool_size);
            
run_name='run_1';
pool_size=30;


params={};
J_vals=0;
h_vals=1;
t_vals=0.1:0.1:1.5;
dh_vals=linspace(0.1,1.5,15);
L_vals=[6,8,10,12];


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
                    p.left_flip=true;
                    p.right_flip=false;
                    p.boundary_cyclic=true;
                    switch L
                        case 6
                            N=200;
                        case 8
                            N=100;
                        case 10
                            N=50;
                        case 12
                            N=25;
                        otherwise
                            N=round(200/2^(L-6));
                    end

                    params{end+1}=p;
                end
            end
        end
    end
end
length(params)
include_all;
parlocal_assert(pool_size);
data_gen_east_model_stats(params,run_name,pool_size);
            
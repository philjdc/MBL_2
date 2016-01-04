run_name='run_3';
pool_size=30;


params={};
J_vals=0;
h_vals=1;
t_vals=0.1:0.1:1.5;
dh_vals=unique([linspace(0.1,1.5,20),logspace(log10(1E-3),log10(1.5),20)]);
L_vals=[14];


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
                    p.N=round(500/2^(L-6));
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
            
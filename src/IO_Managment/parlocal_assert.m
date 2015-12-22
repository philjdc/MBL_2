function []=parlocal_assert(pool_size)
% makes sure local cluster can take the pool size

cluster_prof=parcluster('local');
if (cluster_prof.NumWorkers<pool_size)
    cluster_prof.NumWorkers=pool_size;
    saveProfile(cluster_prof);
end

end
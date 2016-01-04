function agg_data = data_merge(data1,data2)

% initialise agg_data
agg_data=[];
% check same field values
assert(isempty(setxor(fieldnames(data1),fieldnames(data2))));
% get field names
f=fieldnames(data1);
% merge data
for n=1:length(f)
    x1=getfield(data1,f{n});
    x2=getfield(data2,f{n});
    
    agg_data=setfield(agg_data,f{n},[x1,x2]);
end



end


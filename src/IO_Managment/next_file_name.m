function fname = next_file_name(fname_pattern)
% given a string fname_pattern which is a full pathname containg the 
% charachter *, this replaces * with integers returning the fname, the 
% unused filename with the smallest integer in the place of *.

% check template
if (length(strfind(fname_pattern,'*'))~=1)
    error('Error: "next_file_name.m" argument "fname_template" must contain exactly one instance of char "*"');
end

% search
n=1;
while true
    fname=strrep(fname_pattern,'*',num2str(n));
    % check valid
    if ~isempty(regexp(fname, '[*?"<>|]', 'once'))
        error('Error %s is not a valid filename',fname);
    end
    if ~exist(fname,'file')
        break;
    end
    n=n+1;
end

end


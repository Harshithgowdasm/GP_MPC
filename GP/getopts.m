function opts=getopts(opts,newOpts)
    if isempty(newOpts); return; end
    ch1=fieldnames(newOpts);
    ch2=fieldnames(opts);
    for i=1:length(ch1)
        for j=1:length(ch2)
            if strcmp(ch1{i},ch2{j})
                if isstruct(opts.(ch2{j})) && isstruct(newOpts.(ch1{i}))
                    opts.(ch2{j})=getopts(opts.(ch2{j}),newOpts.(ch1{i}));
                else
                    opts.(ch2{j})=newOpts.(ch1{i});
                end
                break;
            end
            if j == length(ch2)
                opts.(ch1{i})=newOpts.(ch1{i});
            end
        end
    end
end


function yaxis_val = raincloud_yaxis(data)
yaxis_val = zeros(size(mean(data,2)));
datum1 = round(mean(data,2),2);
mt_len = unique(datum1);
for i=1:size(mt_len,1)
    found = find(mt_len(i) == datum1);
    if ~isempty(found)
        if size(found,1) == 1
            yaxis_val(found) = 0;
        elseif size(found,1)-1 == 1
            yaxis_val(found) = [0,1].*0.3;
        elseif mod(size(found,1)-1,2) == 0
            yaxis_val(found) = [-((size(found,1)-1)/2):((size(found,1)-1)/2)].*0.3;
        else
            yaxis_val(found) = [-(size(found,1))/2+1:((size(found,1))/2)].*0.3;
        end
    end
end
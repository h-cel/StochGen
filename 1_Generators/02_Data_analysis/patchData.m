function [xpatch, ypatch] = patchData(data,w)
% Nog nodig?

row = size(data,1);
col = size(data,2);
dep = size(data,3);

cnt = 0;
for i = 1:row
    for j = 1:col
        for k = 1:dep
            cnt = cnt + 1;
            if nargin == 1
                [f{cnt},x{cnt}] = ksdensity(data{i,j,k});
            else
                [f{cnt},x{cnt}] = ksdensity(data{i,j,k},'Width',w);
            end
        end
    end
end

f1 = max(cell2mat(f'));
f2 = min(cell2mat(f'));
x1 = max(cell2mat(x'));
x2 = min(cell2mat(x'));

xpatch = [x1 fliplr(x2)];
ypatch = [f1 fliplr(f2)];

end


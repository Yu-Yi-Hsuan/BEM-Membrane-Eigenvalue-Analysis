%% ==========================================================
% 繪製圓形 Im(U) 前十個 singular values 的 Bar Chart
%
% 讀取檔案：
%   ImU_First10Sigma.TXT
%
% 橫軸：1~10，代表 sigma_1 ~ sigma_10
% 縱軸：對應的 singular value
%% ==========================================================

clear; clc; close all;

%% ===== 使用者設定 =====
inputFile = 'ImU_First10Sigma.TXT';
a = 1.0;   % 圓半徑（僅用於圖標題）

%% ===== 讀取檔案 =====
fid = fopen(inputFile,'r');

if fid == -1
    error('找不到檔案：%s',inputFile);
end

sigma = nan(10,1);
kValue = NaN;

while ~feof(fid)

    lineText = strtrim(fgetl(fid));

    if isempty(lineText)
        continue;
    end

    % 讀取 K
    if startsWith(lineText,'K =')
        temp = sscanf(lineText,'K = %f');
        if ~isempty(temp)
            kValue = temp(1);
        end
    end

    % 讀取前十個 sigma
    temp = sscanf(lineText,'%d %f');

    if numel(temp) == 2
        id = temp(1);

        if id >= 1 && id <= 10
            sigma(id) = temp(2);
        end
    end
end

fclose(fid);

if any(isnan(sigma))
    error('無法完整讀取前十個 sigma 值，請檢查 %s 格式。',inputFile);
end

%% ===== 顯示數值 =====
fprintf('K = %.10f\n',kValue);

for i = 1:10
    fprintf('sigma%d = %.15e\n',i,sigma(i));
end

%% ===== Bar Chart =====
figure;

bar(1:10,sigma,0.65);

% 因為十個 sigma 的量級可能差很多，使用 log scale 較容易比較
set(gca,'YScale','log');

xlabel('Singular value number');
ylabel('\sigma');

if isfinite(kValue)
    title(sprintf('Circle, a = %g, Im(U), k = %.4f',a,kValue));
else
    title(sprintf('Circle, a = %g, First Ten Singular Values of Im(U)',a));
end

xticks(1:10);
xticklabels(string(1:10));

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'LineWidth',1.2);

grid off;
box on;

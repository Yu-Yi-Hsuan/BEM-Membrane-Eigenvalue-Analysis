clc; clear; close all;

%% ===== 使用者設定 =====
a = 1.0;                % 圓半徑（僅用於圖標題）
valleyWindow = 0.10;    % 谷值左右比較範圍
dropRatio = 1.5;        % 最小下降倍率

%% ===== 讀取 Im(U) 資料 =====
fileName = 'Im_U.DAT';

if ~isfile(fileName)
    error('找不到檔案：%s',fileName);
end

data = load(fileName);

% 至少需要兩欄：k、sigma_min
if size(data,2) < 2
    error('%s 必須包含 k 與最小奇異值兩欄。',fileName);
end

k = data(:,1);
sigma = abs(data(:,2));

%% ===== 移除無效資料 =====
valid = isfinite(k) & k > 0 & ...
        isfinite(sigma) & sigma > 0;

k = k(valid);
sigma = sigma(valid);

%% ===== 依 k 排序 =====
[k,order] = sort(k);
sigma = sigma(order);

%% ===== 尋找並列出谷值 =====
idx = findDeepValleys(k,sigma,valleyWindow,dropRatio);

fprintf('\nIm(U) 谷值 k：\n');
disp(k(idx));

%% ===== 繪製 Im(U) =====
figure;

hImU = semilogy(k,sigma,'-', ...
    'LineWidth',2.0, ...
    'DisplayName','Im(U)');

xlabel('k');
ylabel('\sigma_{min}');

legend(hImU,'Location','best');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'LineWidth',1.2);

axis tight;
box on;
grid off;

disp('請調整圖形位置與大小，完成後回 Command Window 按 Enter 儲存...');
pause;
exportgraphics(gcf,'Fig_Circle_ImU.tif', ...
 'Resolution',2400);
%% ==========================================================
% 尋找具有明顯下降的局部谷值
%% ==========================================================
function idx = findDeepValleys(k,sigma,valleyWindow,dropRatio)

    candidate = find(islocalmin(sigma));
    idx = false(size(sigma));

    for n = 1:numel(candidate)
        i = candidate(n);

        left  = k >= k(i)-valleyWindow & k < k(i);
        right = k > k(i) & k <= k(i)+valleyWindow;

        if ~any(left) || ~any(right)
            continue;
        end

        sideLevel = min(max(sigma(left)),max(sigma(right)));

        if sideLevel/max(sigma(i),eps) >= dropRatio
            idx(i) = true;
        end
    end
end

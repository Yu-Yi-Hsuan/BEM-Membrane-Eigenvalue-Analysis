%% ==========================================================
% 圓形域：比較 Re(U)、Im(U)、Complex U 的最小奇異值
% MATLAB R2024a
%% ==========================================================

clc; clear; close all;

%% ===== 使用者設定 =====
valleyWindow = 0.15;    % 谷值左右比較範圍
dropRatio    = 1.5;     % 最小下降倍率

%% ===== 讀取資料 =====
fileNames = {
    'Re_U.DAT'
    'Im_U.DAT'
    'Complex_U.DAT'
};

labels = {
    'Re(U)'
    'Im(U)'
    'Complex(U)'
};

% 檢查檔案是否存在
for i = 1:numel(fileNames)
    if ~isfile(fileNames{i})
        error('找不到檔案：%s',fileNames{i});
    end
end

data = cellfun(@load,fileNames,'UniformOutput',false);

% 每個檔案至少須包含 k 與 sigma_min 兩欄
if any(cellfun(@(A) size(A,2) < 2,data))
    error('Re_U.DAT、Im_U.DAT、Complex_U.DAT 都必須包含 k 與最小奇異值兩欄。');
end

% 以 Complex_U.DAT 的 k 為基準
k = data{3}(:,1);
sigma = cellfun(@(A) abs(A(:,2)),data,'UniformOutput',false);

% 檢查資料筆數及 k 值是否一致
if any(cellfun(@(A) size(A,1) ~= numel(k),data))
    error('三個 DAT 檔的資料筆數不一致。');
end

for i = 1:numel(data)
    if any(abs(data{i}(:,1)-k) > 1.0e-10)
        error('三個 DAT 檔的 k 值不一致。');
    end
end

% 移除無效資料
valid = isfinite(k) & k > 0;
for i = 1:numel(sigma)
    valid = valid & isfinite(sigma{i}) & sigma{i} > 0;
end

k = k(valid);
for i = 1:numel(sigma)
    sigma{i} = sigma{i}(valid);
end

% 依 k 排序
[k,order] = sort(k);
for i = 1:numel(sigma)
    sigma{i} = sigma{i}(order);
end

%% ===== 尋找並列出谷值 =====
idx = cell(numel(sigma),1);

for i = 1:numel(sigma)
    idx{i} = findDeepValleys(k,sigma{i},valleyWindow,dropRatio);

    fprintf('\n%s 谷值 k：\n',labels{i});
    disp(k(idx{i}));
end

%% ===== 繪製比較圖 =====
figure;

% 先畫 Complex U，再畫 Im(U)，最後畫 Re(U)，
% 讓 Re(U) 在曲線重疊時顯示於最上層。
hComplex = semilogy(k,sigma{3},'-','LineWidth',2.5);
hold on;
hIm = semilogy(k,sigma{2},'-.','LineWidth',1.5);
hRe = semilogy(k,sigma{1},'--','LineWidth',1.8);

formatPlot();

% 明確指定圖例控制物件，避免圖例線型與名稱錯位。
legend([hRe,hIm,hComplex],labels,'Location','best');

%% ==========================================================
% 尋找具有明顯下降的局部谷值
%% ==========================================================
function idx = findDeepValleys(k,sigma,valleyWindow,dropRatio)

    candidate = find(islocalmin(sigma));
    idx = false(size(sigma));

    for n = 1:numel(candidate)
        i = candidate(n);

        left = k >= k(i)-valleyWindow & k < k(i);
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

%% ==========================================================
% 統一圖形格式
%% ==========================================================
function formatPlot()

    xlabel('k');
    ylabel('\sigma_{min}');

    set(gca, ...
        'FontName','Times New Roman', ...
        'FontSize',12, ...
        'LineWidth',1.2);

    axis tight;
    box on;
    grid off;
    
end

disp('請調整圖形位置與大小，完成後回 Command Window 按 Enter 儲存...');
pause;

exportgraphics(gcf,'Fig_Circle.tif', ...
'Resolution',2400);

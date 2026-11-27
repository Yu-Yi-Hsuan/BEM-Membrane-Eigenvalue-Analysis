clc; clear; close all;

%% ===== 使用者設定 =====
N = 7;                 % (N,k)-gon 的 N
kgon = 3;              % (N,k)-gon 的 k
valleyWindow = 0.20;   % 谷值左右比較範圍
dropRatio = 1.5;       % 最小下降倍率

%% ===== 讀取資料 =====
fileNames = {
    'Complex_U.DAT'
    'Complex_L.DAT'
    'ReU_ReL.DAT'
    'Complex_U_L.DAT'
};

data = cellfun(@load,fileNames,'UniformOutput',false);

% 每個檔案皆須包含 k 與 sigma_min
if any(cellfun(@(A) size(A,2) < 2,data))
    error('四個 DAT 檔都必須包含 k 與最小奇異值兩欄。');
end

k = data{4}(:,1);
sigma = cellfun(@(A) abs(A(:,2)),data,'UniformOutput',false);

% 檢查資料筆數
if any(cellfun(@(A) size(A,1) ~= numel(k),data))
    error('四個 DAT 檔的資料筆數不一致。');
end

% 移除無效資料
valid = isfinite(k) & k > 0;

for i = 1:4
    valid = valid & isfinite(sigma{i}) & sigma{i} > 0;
end

k = k(valid);

for i = 1:4
    sigma{i} = sigma{i}(valid);
end

% 依 k 排序
[k,order] = sort(k);

for i = 1:4
    sigma{i} = sigma{i}(order);
end

%% ===== 尋找谷值 =====
labels = {
    'Complex U'
    'Complex L'
    '[Re(U);Re(L)]'
    'Complex [U;L]'
};

idx = cell(4,1);

for i = 1:4
    idx{i} = findDeepValleys( ...
        k,sigma{i},valleyWindow,dropRatio);

    fprintf('\n%s 谷值 k：\n',labels{i});
    disp(k(idx{i}));
end

%% ===== 圖 1：U、L 與複數增廣矩陣 =====
figure;

semilogy(k,sigma{1},'--','LineWidth',1.5);
hold on;
semilogy(k,sigma{2},'-.','LineWidth',1.5);
semilogy(k,sigma{4},'-','LineWidth',2);

formatPlot(N,kgon);
legend(labels([1,2,4]),'Location','best');

%% ===== 圖 2：實部與複數增廣矩陣 =====
figure;

semilogy(k,sigma{3},'--','LineWidth',1.5);
hold on;
semilogy(k,sigma{4},'-','LineWidth',2);

formatPlot(N,kgon);
legend(labels([3,4]),'Location','best');

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
function formatPlot(N,kgon)

    xlabel('k');
    ylabel('\sigma_{min}');
    title(sprintf('(%d,%d)-gon',N,kgon));

    set(gca, ...
        'FontName','Times New Roman', ...
        'FontSize',12, ...
        'LineWidth',1.2);

    axis tight;
    box on;
    grid off;
end

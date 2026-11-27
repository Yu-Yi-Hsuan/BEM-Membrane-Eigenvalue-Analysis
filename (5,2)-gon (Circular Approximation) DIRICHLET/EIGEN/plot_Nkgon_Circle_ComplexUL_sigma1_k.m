clc; clear; close all;

%% ===== 使用者設定 =====
% 雙聯通幾何說明
% 這些文字只用於圖標題，不影響 Complex_U_L.DAT 的讀取與繪圖。
%
% 目前幾何：
%   外邊界：(5,2)-gon
%   內邊界：circle

boundaryInfo = {
    '(5,2)-gon';    % 外邊界
    'circle'        % 內圓
};

valleyWindow = 0.15;   % 谷值左右比較範圍
dropRatio    = 1.5;    % 最小下降倍率

%% ===== 讀取 Complex [U;L] 資料 =====
inputFile = 'Complex_U_L.DAT';

if ~isfile(inputFile)
    error('找不到檔案：%s',inputFile);
end

CUL = load(inputFile);

if size(CUL,2) < 2
    error('%s 必須包含 k 與最小奇異值兩欄。',inputFile);
end

k        = CUL(:,1);
sigmaCUL = abs(CUL(:,2));

%% ===== 移除無效資料並排序 =====
valid = isfinite(k) & k > 0 & ...
        isfinite(sigmaCUL) & sigmaCUL > 0;

k        = k(valid);
sigmaCUL = sigmaCUL(valid);

[k,order] = sort(k);
sigmaCUL  = sigmaCUL(order);

%% ===== 尋找谷值 =====
idxCUL = findDeepValleys( ...
    k,sigmaCUL,valleyWindow,dropRatio);

fprintf('\nComplex [U;L] 的真根候選 k：\n');
disp(k(idxCUL));

%% ===== 繪圖 =====
figure;

hCUL = semilogy(k,sigmaCUL,'-', ...
    'LineWidth',2.0, ...
    'DisplayName','Complex [U;L]');

xlabel('k');
ylabel('\sigma_{min}');

% 自動產生多連通幾何標題
titleText = strjoin(boundaryInfo,' / ');
title(['Multiply-connected domain: ',titleText]);

legend(hCUL,'Location','best');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'LineWidth',1.2);

axis tight;
box on;
grid off;

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

        sideLevel = min( ...
            max(sigma(left)), ...
            max(sigma(right)));

        if sideLevel/max(sigma(i),eps) >= dropRatio
            idx(i) = true;
        end
    end
end

clc; clear; close all;

%% ===== 使用者設定 =====
N = 7;                 % (N,k)-gon 的 N
kgon = 3;              % (N,k)-gon 的 k
valleyWindow = 0.15;   % 谷值左右比較範圍
dropRatio = 1.5;       % 最小下降倍率

%% ===== 讀取資料 =====
CU = load('Complex_U.DAT');

if size(CU,2) < 2
    error('Complex_U.DAT 必須包含 k 與最小奇異值兩欄。');
end

k = CU(:,1);
sigmaCU = abs(CU(:,2));

%% ===== 移除無效資料並排序 =====
valid = isfinite(k) & k > 0 & isfinite(sigmaCU) & sigmaCU > 0;
k = k(valid);
sigmaCU = sigmaCU(valid);

[k,order] = sort(k);
sigmaCU = sigmaCU(order);

%% ===== 尋找谷值 =====
idxCU = findDeepValleys(k,sigmaCU,valleyWindow,dropRatio);

fprintf('\nComplex U 的真根候選 k：\n');
disp(k(idxCU));

%% ===== 繪圖 =====
figure;
semilogy(k,sigmaCU,'-','LineWidth',2.0);

formatPlot(N,kgon);
legend('Complex(U)','Location','best');

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

exportgraphics(gcf,'Fig.tif', ...
'Resolution',2400);
%% ==========================================================
% 環形域四種 SVD 搜根結果比較
%
% 圖 1：Complex U、Complex L、Complex [U;L]
% 圖 2：[Re(U);Re(L)]、Complex [U;L]
%
% 只列出具有明顯下降的局部谷值
%% ==========================================================

clear; clc; close all;

%% ===== 使用者設定 =====
valleyWindow = 0.15;   % 左右比較範圍
dropRatio    = 1.5;    % 至少下降倍率

%% ===== 讀取資料 =====
CU = load('Complex_U.DAT');

k = CU(:,1);
sigmaCU = abs(CU(:,2));

%% ===== 尋找明顯谷值 =====
idxCU = findDeepValleys(k,sigmaCU,valleyWindow,dropRatio);

fprintf('\nComplex [U] 的真根候選 k：\n');
disp(k(idxCU));

%% ===== 繪圖 =====
figure;
semilogy(k,sigmaCU,'-','LineWidth',2.0);

xlabel('k');
ylabel('\sigma_{min}');
legend('Complex [U]','Location','best');

set(gca,'FontName','Times New Roman',...
          'FontSize',12,...
          'LineWidth',1.2);

axis tight;
box on;
grid off;

%% ==========================================================
% 保留具有明顯下降的局部谷值
%
% valleyWindow：左右比較範圍
% dropRatio   ：最低點至少下降倍率
%% ==========================================================
function idx = findDeepValleys(k,sigma,valleyWindow,dropRatio)

candidate = find(islocalmin(sigma));
idx = false(size(sigma));

for n = 1:length(candidate)

    i = candidate(n);

    left  = k>=k(i)-valleyWindow & k<k(i);
    right = k>k(i) & k<=k(i)+valleyWindow;

    if ~any(left) || ~any(right)
        continue;
    end

    peak = min(max(sigma(left)),max(sigma(right)));

    if peak/max(sigma(i),eps) >= dropRatio
        idx(i)=true;
    end

end

end
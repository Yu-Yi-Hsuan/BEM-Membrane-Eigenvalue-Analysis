%% ==========================================================
% 比較四種 SVD 搜根結果
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
U    = load('Complex_U.DAT');
L    = load('Complex_L.DAT');
ReUL = load('ReU_ReL.DAT');
CUL  = load('Complex_U_L.DAT');

k = CUL(:,1);

sigmaU    = abs(U(:,2));
sigmaL    = abs(L(:,2));
sigmaReUL = abs(ReUL(:,2));
sigmaCUL  = abs(CUL(:,2));

%% ===== 尋找明顯谷值 =====
idxU    = findDeepValleys(k,sigmaU,valleyWindow,dropRatio);
idxL    = findDeepValleys(k,sigmaL,valleyWindow,dropRatio);
idxReUL = findDeepValleys(k,sigmaReUL,valleyWindow,dropRatio);
idxCUL  = findDeepValleys(k,sigmaCUL,valleyWindow,dropRatio);

fprintf('\nComplex U 的谷值 k：\n');
disp(k(idxU));

fprintf('Complex L 的谷值 k：\n');
disp(k(idxL));

fprintf('[Re(U);Re(L)] 的谷值 k：\n');
disp(k(idxReUL));

fprintf('Complex [U;L] 的真根候選 k：\n');
disp(k(idxCUL));

%% ===== 圖 1：單一矩陣與複數增廣矩陣 =====
figure;

semilogy(k,sigmaU,'--','LineWidth',1.5);
hold on;
semilogy(k,sigmaL,'-.','LineWidth',1.5);
semilogy(k,sigmaCUL,'-','LineWidth',2.0);

xlabel('k');
ylabel('\sigma_{min}');

legend('Complex U','Complex L','Complex [U;L]', ...
       'Location','best');

set(gca,'FontName','Times New Roman', ...
         'FontSize',12, ...
         'LineWidth',1.2);

axis tight;
box on;
grid off;

%% ===== 圖 2：實部與複數增廣矩陣 =====
figure;

semilogy(k,sigmaReUL,'--','LineWidth',1.5);
hold on;
semilogy(k,sigmaCUL,'-','LineWidth',2.0);

xlabel('k');
ylabel('\sigma_{min}');

legend('[Re(U);Re(L)]','Complex [U;L]', ...
       'Location','best');

set(gca,'FontName','Times New Roman', ...
         'FontSize',12, ...
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
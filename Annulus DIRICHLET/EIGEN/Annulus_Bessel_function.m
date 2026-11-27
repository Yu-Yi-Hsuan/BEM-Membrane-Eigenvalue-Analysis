%% ==========================================================
% 環形 Dirichlet 特徵波數解析解比對程式
%
% 解析方程：
% J_n(k*r1)Y_n(k*r2) - J_n(k*r2)Y_n(k*r1) = 0
%
% 功能：
% 1. 輸入數值波數 k_num
% 2. 搜尋最接近的環形解析波數
% 3. 顯示對應階數 n、根的型態及百分比誤差
%
% MATLAB R2024a
%% ==========================================================

clear; clc;

%% ===== 使用者設定 =====
k_num = 5.52;      % 數值程式得到的波數 k
r1 = 0.43565;           % 內半徑
r2 = 1.0;           % 外半徑

nMax = 20;          % 搜尋最高 Bessel 階數
kMax = 20;          % 搜尋解析根的最大 k
dk = 0.001;         % 掃描間距

%% ===== 參數檢查 =====
if r1 <= 0 || r2 <= r1
    error('需滿足 0 < r1 < r2。');
end

if k_num <= 0
    error('k_num 必須大於 0。');
end

%% ===== 搜尋環形解析根 =====
rootTable = [];

for n = 0:nMax

    f = @(k) besselj(n,k*r1).*bessely(n,k*r2) ...
           - besselj(n,k*r2).*bessely(n,k*r1);

    kScan = dk:dk:kMax;
    fScan = f(kScan);

    for i = 1:length(kScan)-1

        if isfinite(fScan(i)) && isfinite(fScan(i+1)) ...
                && fScan(i)*fScan(i+1) < 0

            kRoot = fzero(f,[kScan(i),kScan(i+1)]);

            if isempty(rootTable) || ...
                    all(abs(rootTable(:,2)-kRoot) > 1e-6)
                rootTable = [rootTable; n, kRoot]; %#ok<AGROW>
            end
        end
    end
end

if isempty(rootTable)
    error('找不到解析根，請增加 kMax 或 nMax。');
end

%% ===== 找出最接近的解析根 =====
[~,idx] = min(abs(rootTable(:,2)-k_num));

nMatch = rootTable(idx,1);
kExact = rootTable(idx,2);

if nMatch == 0
    rootType = '單根';
else
    rootType = '二重根';
end

percentError = abs(k_num-kExact)/abs(kExact)*100;

%% ===== 顯示結果 =====
fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('環形 Dirichlet 特徵波數比對結果\n');
fprintf('----------------------------------------\n');
fprintf('內半徑 r1       = %.6f\n',r1);
fprintf('外半徑 r2       = %.6f\n',r2);
fprintf('輸入數值 k      = %.6f\n',k_num);
fprintf('最接近階數 n    = %d\n',nMatch);
fprintf('解析波數 k      = %.10f\n',kExact);
fprintf('根的型態        = %s\n',rootType);
fprintf('百分比誤差      = %.6f %%\n',percentError);
fprintf('----------------------------------------\n');
%% ========================================================== 
% 自動比對輸入波數 k 與圓形解析解 Bessel 零點
%
% 功能：
% 1. 輸入數值波數 k 與半徑 a
% 2. 自動尋找最接近的 Bessel 階數 n
% 3. 回傳對應解析波數 k
% 4. 計算 J_n(k*a)
% 5. 判斷單根或二重根
% 6. 計算百分比誤差
%% ==========================================================

clear; clc;

%% ====== 使用者輸入 ======
k_num = 5.5201;     % 你算到的波數 k
a = 0.43565;            % 圓半徑

%% ====== 搜尋設定 ======
nmax = 10;          % 搜尋 J0 ~ J10
mmax = 10;          % 每一階最多找前 10 個零點
xmax = 50;          % 搜尋範圍

ka_num = k_num * a;

best_error = inf;

%% ====== 搜尋所有 Bessel 零點 ======
for n = 0:nmax

    x = linspace(1e-6, xmax, 200000);
    J = besselj(n,x);

    roots = [];

    for i = 1:length(x)-1

        if J(i)*J(i+1) < 0

            root = fzero(@(t)besselj(n,t), [x(i),x(i+1)]);
            roots = [roots; root];

            if length(roots) >= mmax
                break;
            end

        end

    end

    %% 找最接近輸入 k 的解析根
    for j = 1:length(roots)

        k_exact = roots(j) / a;

        error_percent = abs(k_num-k_exact) / abs(k_exact) * 100;

        if error_percent < best_error
            best_error = error_percent;
            best_n = n;
            best_k = k_exact;
            J_input = besselj(n,ka_num);
        end

    end

end

%% ====== 判斷重根性 ======
if best_n == 0
    multiplicity_text = '單根';
else
    multiplicity_text = '二重根';
end

%% ====== 顯示結果 ======
fprintf('=====================================\n');
fprintf('輸入波數 k          = %.10f\n',k_num);
fprintf('半徑 a              = %.10f\n\n',a);

fprintf('最接近 Bessel 階數 n = %d\n',best_n);
fprintf('解析波數 k          = %.10f\n',best_k);
fprintf('J_%d(k*a)            = %.15e\n\n',best_n,J_input);

fprintf('重根性              = %s\n',multiplicity_text);
fprintf('百分比誤差          = %.6f %%\n',best_error);
fprintf('=====================================\n');
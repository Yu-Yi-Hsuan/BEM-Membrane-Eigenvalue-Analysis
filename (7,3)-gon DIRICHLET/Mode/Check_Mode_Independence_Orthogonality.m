%% 模態獨立性與正交性檢查
% 比較兩個實數模態的線性獨立性與離散正交性。
% 請將程式與下列檔案放在 MATLAB 目前工作資料夾：
%   xspmoU2.dat、xspmoU2_2.dat：域內座標及模態值 [x, y, u]
%   xspedU2.DAT、xspedU2_2.DAT：邊界元素編號及模態值 [編號, u]
% 檔名含 _2 的資料須對應同一特徵波數的第二個模態。
% 不必重新執行 Fortran，且不會建立額外資料夾。
clear; clc; close all;

%% ===== 1. 使用者設定 =====
workDir = pwd;             % MATLAB 目前工作資料夾
checkBoundary = true;      % 是否一併檢查邊界模態
coordTol = 1e-7;           % 兩份域內資料的座標比對容許誤差
rankTolRel = 1e-8;        % 數值秩判定的相對門檻
orthTol = 1e-3;           % 正規化內積低於此值，視為近似正交
saveReport = false;       % 設為 true：將結果另存文字檔（不建立資料夾）

%% ===== 2. 讀取域內模態 =====
inside1 = findFile(workDir, 'xspmoU2.dat');
inside2 = findFile(workDir, 'xspmoU2_2.dat');
raw1 = readmatrix(inside1);
raw2 = readmatrix(inside2);
validateattributes(raw1,{'numeric'},{'2d','ncols',3,'nonempty','finite'});
validateattributes(raw2,{'numeric'},{'2d','ncols',3,'nonempty','finite'});

if size(raw1,1) ~= size(raw2,1)
    error('兩份域內模態資料筆數不同：%d 與 %d。',size(raw1,1),size(raw2,1));
end

% 依 x、y 座標排序，避免兩份檔案的場點順序不同而造成錯誤比較。
[xy1, p1] = sortrows(raw1(:,1:2),[1 2]);
[xy2, p2] = sortrows(raw2(:,1:2),[1 2]);
maxCoordDiff = max(abs(xy1(:)-xy2(:)));
if maxCoordDiff > coordTol
    error('域內座標不一致：最大差異 %.4g，大於容許誤差 %.4g。',maxCoordDiff,coordTol);
end

u1 = raw1(p1,3);
u2 = raw2(p2,3);
resultInside = checkTwoModes(u1,u2,rankTolRel,orthTol);

fprintf('\n========== 域內模態檢查 ==========\n');
fprintf('檔案：%s / %s\n',inside1,inside2);
fprintf('共同場點數：%d\n',numel(u1));
fprintf('最大座標差異：%.8g\n',maxCoordDiff);
printResult(resultInside);
fprintf('注意：此處為等權重離散內積，不等同於嚴格的域內質量積分。\n');

%% ===== 3. 讀取邊界模態（可選） =====
resultBoundary = [];
if checkBoundary
    try
        bound1 = findFile(workDir,'xspedU2.DAT');
        bound2 = findFile(workDir,'xspedU2_2.DAT');
        b1 = readmatrix(bound1);
        b2 = readmatrix(bound2);
        validateattributes(b1,{'numeric'},{'2d','ncols',2,'nonempty','finite'});
        validateattributes(b2,{'numeric'},{'2d','ncols',2,'nonempty','finite'});

        if size(b1,1) ~= size(b2,1)
            error('兩份邊界模態的元素數不同。');
        end

        % 依元素編號排列，確保兩個模態的元素一一對應。
        [id1,q1] = sort(b1(:,1));
        [id2,q2] = sort(b2(:,1));
        if ~isequal(id1,id2)
            error('兩份邊界模態的元素編號不一致。');
        end

        resultBoundary = checkTwoModes(b1(q1,2),b2(q2,2),rankTolRel,orthTol);
        fprintf('\n========== 邊界模態檢查 ==========\n');
        fprintf('檔案：%s / %s\n',bound1,bound2);
        fprintf('共同邊界元素數：%d\n',numel(id1));
        printResult(resultBoundary);
        fprintf('注意：此處為未加權離散內積，不等同於邊界積分。\n');
    catch ME
        warning('已略過邊界模態檢查：%s',ME.message);
    end
end

%% ===== 4. 輸出文字報告（可選） =====
if saveReport
    reportFile = fullfile(workDir,'Mode_Independence_Report.txt');
    fid = fopen(reportFile,'w');
    if fid < 0
        error('無法建立文字報告檔案。');
    end

    fprintf(fid,'域內模態：場點數 = %d，數值秩 = %d，s2/s1 = %.12g，正規化內積 = %.12g\n', ...
        numel(u1),resultInside.rank,resultInside.sRatio,resultInside.eta);
    if ~isempty(resultBoundary)
        fprintf(fid,'邊界模態：數值秩 = %d，s2/s1 = %.12g，正規化內積 = %.12g\n', ...
            resultBoundary.rank,resultBoundary.sRatio,resultBoundary.eta);
    end
    fclose(fid);
    fprintf('\n文字報告已儲存：%s\n',reportFile);
end

%% ===== 5. 附屬函式 =====
function filepath = findFile(folder,wanted)
    % 忽略副檔名大小寫差異，例如 .dat 和 .DAT 均可讀取。
    listing = dir(folder);
    idx = find(~[listing.isdir] & strcmpi({listing.name},wanted));
    if numel(idx) ~= 1
        error('資料夾 %s 中應有且僅有一份檔案 %s，目前找到 %d 份。',folder,wanted,numel(idx));
    end
    filepath = fullfile(folder,listing(idx).name);
end

function R = checkTwoModes(a,b,rankTolRel,orthTol)
    a = double(a(:));
    b = double(b(:));
    if any(~isfinite(a)) || any(~isfinite(b))
        error('模態資料包含非有限數值。');
    end

    na = norm(a);
    nb = norm(b);
    if na <= eps || nb <= eps
        R = struct('rank',0,'sRatio',NaN,'eta',NaN,'norm1',na,'norm2',nb, ...
                   'independent',false,'orthogonal',false);
        return;
    end

    % 先將兩個模態分別正規化，避免振幅大小影響數值秩的判定。
    X = [a/na, b/nb];
    s = svd(X,'econ');
    R.sRatio = s(2)/s(1);
    R.rank = sum(s > rankTolRel*s(1));
    R.eta = abs(X(:,1)'*X(:,2));
    R.norm1 = na;
    R.norm2 = nb;
    R.independent = (R.rank == 2);
    R.orthogonal = (R.eta < orthTol);
end

function printResult(R)
    fprintf('模態 1 / 模態 2 的範數：%.10g / %.10g\n',R.norm1,R.norm2);
    fprintf('數值秩（兩個正規化模態）：%d\n',R.rank);
    fprintf('奇異值比值 s2/s1：%.12e\n',R.sRatio);
    fprintf('正規化內積絕對值：%.12e\n',R.eta);
    if R.independent
        fprintf('線性獨立性：通過\n');
    else
        fprintf('線性獨立性：未通過或數值精度不足\n');
    end
    if R.orthogonal
        fprintf('離散近似正交性：通過\n');
    else
        fprintf('離散近似正交性：未通過\n');
    end
end

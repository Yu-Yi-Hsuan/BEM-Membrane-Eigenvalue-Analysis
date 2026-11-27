%% ==========================================================
% 正三角形 Dirichlet 特徵波數解析解比對程式
%
% 解析公式：
% k_mn = 4*pi/(3*a)*sqrt(m^2 + m*n + n^2)
%
% MATLAB R2024a
%% ==========================================================

clear; clc;

%% ===== 使用者設定 =====
k_num = 8.7180;          % 數值程式得到的波數 k
inputFile = 'F15.DAT';

kMax = 20;               % 搜尋解析根的最大 k

%% ===== 參數檢查 =====
if k_num <= 0
    error('k_num 必須大於 0。');
end

%% ===== 讀取 F15 邊界座標 =====
fid = fopen(inputFile, 'r');

if fid == -1
    error('無法開啟 %s。', inputFile);
end

% 跳過前兩行控制資料
fgetl(fid);
fgetl(fid);

xb = [];
yb = [];

while ~feof(fid)
    lineText = strtrim(fgetl(fid));

    if isempty(lineText)
        continue;
    end

    % 節點資料以 -1 結束
    if strcmp(lineText, '-1')
        break;
    end

    data = sscanf(lineText, '%f');

    % F15 第 5、6 欄為 x、y 座標
    if numel(data) >= 6
        xb(end+1,1) = data(5); %#ok<AGROW>
        yb(end+1,1) = data(6); %#ok<AGROW>
    end
end

fclose(fid);

if numel(xb) < 3
    error('邊界點不足，請檢查 F15.DAT。');
end

%% ===== 取得正三角形三個頂點 =====
points = [xb, yb];
nPoint = size(points,1);

distanceMatrix = zeros(nPoint);

% 計算所有邊界點之間的距離
for i = 1:nPoint
    for j = i+1:nPoint
        distance = hypot(points(i,1)-points(j,1), ...
                         points(i,2)-points(j,2));

        distanceMatrix(i,j) = distance;
        distanceMatrix(j,i) = distance;
    end
end

% 取距離最遠的兩個頂點
[~,index] = max(distanceMatrix(:));
[i1,i2] = ind2sub(size(distanceMatrix),index);

% 取得第三個頂點
distanceSum = distanceMatrix(:,i1) + distanceMatrix(:,i2);
[~,i3] = max(distanceSum);

vertices = points([i1,i2,i3],:);

%% ===== 計算三角形平均邊長 =====
d12 = norm(vertices(1,:) - vertices(2,:));
d23 = norm(vertices(2,:) - vertices(3,:));
d31 = norm(vertices(3,:) - vertices(1,:));

a = mean([d12,d23,d31]);

%% ===== 搜尋正三角形解析根 =====
mMax = ceil(3*a*kMax/(4*pi)) + 2;

modeTable = [];

% Dirichlet 模態索引：n >= m > 0
for m = 1:mMax
    for n = m:mMax
        ell = m^2 + m*n + n^2;
        kRoot = 4*pi/(3*a)*sqrt(ell);

        if kRoot <= kMax
            % m = n：只有對稱模態，重數貢獻 1
            % m < n：對稱與反對稱模態，重數貢獻 2
            if m == n
                modeCount = 1;
            else
                modeCount = 2;
            end

            modeTable = [modeTable; ...
                m, n, ell, kRoot, modeCount]; %#ok<AGROW>
        end
    end
end

if isempty(modeTable)
    error('找不到解析根，請增加 kMax。');
end

%% ===== 合併相同解析特徵波數 =====
ellUnique = unique(modeTable(:,3));

kRoot = zeros(numel(ellUnique),1);
multiplicity = zeros(numel(ellUnique),1);

for i = 1:numel(ellUnique)
    idxSame = modeTable(:,3) == ellUnique(i);

    kRoot(i) = modeTable(find(idxSame,1),4);
    multiplicity(i) = sum(modeTable(idxSame,5));
end

[kRoot,order] = sort(kRoot);
multiplicity = multiplicity(order);

%% ===== 找出最接近的解析根 =====
[~,idx] = min(abs(kRoot-k_num));

kExact = kRoot(idx);
rootMultiplicity = multiplicity(idx);

if rootMultiplicity == 1
    rootType = '單根';
elseif rootMultiplicity == 2
    rootType = '二重根';
else
    rootType = sprintf('%d 重根',rootMultiplicity);
end

percentError = abs(k_num-kExact)/abs(kExact)*100;

%% ===== 顯示結果 =====
fprintf('\n');
fprintf('----------------------------------------\n');
fprintf('正三角形 Dirichlet 特徵波數比對結果\n');
fprintf('----------------------------------------\n');
fprintf('平均邊長 a      = %.6f\n',a);
fprintf('輸入數值 k      = %.6f\n',k_num);
fprintf('解析波數 k      = %.10f\n',kExact);
fprintf('根的型態        = %s\n',rootType);
fprintf('百分比誤差      = %.6f %%\n',percentError);
fprintf('----------------------------------------\n');

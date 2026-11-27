%% 正多邊形薄膜特徵值分析（FEM）
% 使用三節點線性三角形元素，求解二維 Helmholtz 固定邊界特徵值問題：
%   -∇²u = k²u   （域內）
%       u = 0     （邊界）
% 離散後求解 K*u = lambda*M*u，其中 lambda = k²。
clear; clc; close all;

%% ===== 0. 使用者設定 =====
nSide     = 5;       % 正多邊形邊數，例如 3、4、5、6
R         = 0.8506508084;     % 外接圓半徑（中心到頂點的距離）
h         = 0.01;    % 目標網格尺寸（越小，元素越多）
nModes    = 10;      % 要計算的特徵值數量（重根也分別計算）
angle0    = pi/2;    % 第一個頂點的角度（弧度），只影響旋轉方向
plotModes = 6;      % 要繪製的模態數量

% 基本參數檢查
assert(isscalar(nSide) && nSide>=3 && nSide==floor(nSide), 'nSide 必須是至少 3 的整數。');
assert(isscalar(R) && isfinite(R) && R>0, 'R 必須為正數。');
assert(isscalar(h) && isfinite(h) && h>0, 'h 必須為正數。');
assert(isscalar(nModes) && isfinite(nModes) && nModes>=1 && nModes==floor(nModes), 'nModes 必須為正整數。');
assert(isscalar(angle0) && isfinite(angle0), 'angle0 必須為有限數值。');
assert(isscalar(plotModes) && isfinite(plotModes) && plotModes>=0 && plotModes==floor(plotModes), 'plotModes 必須為非負整數。');


%% ===== 1. 建立正多邊形頂點 =====
theta = angle0 + (0:nSide-1)'*2*pi/nSide;
V = R*[cos(theta),sin(theta)];

%% ===== 2. 建立三角形網格 =====
% 從中心連接各相鄰頂點，形成 nSide 個大三角形；
% 每個大三角形再規則細分為 nDiv^2 個小三角形。
sideLength = 2*R*sin(pi/nSide);  % 正多邊形每條邊的長度
nDiv = max(2,ceil(max(R,sideLength)/h)); % 每個大三角形的分割數

% 預先配置空間，避免在迴圈中反覆擴增陣列。
nRaw = nSide*(nDiv+1)*(nDiv+2)/2;
nElem = nSide*nDiv^2;
Praw = zeros(nRaw,2);
Traw = zeros(nElem,3);
pCount = 0;
tCount = 0;
for s = 1:nSide
    A = V(s,:);
    B = V(mod(s,nSide)+1,:);
    loc = zeros(nDiv+1,nDiv+1);  % 局部節點編號（i+1,j+1）

    % 網格點位置：P=(i*A+j*B)/nDiv，條件 i+j<=nDiv
    for i = 0:nDiv
        for j = 0:(nDiv-i)
            pCount = pCount+1;
            Praw(pCount,:) = (i*A+j*B)/nDiv;
            loc(i+1,j+1) = pCount;
        end
    end

    % 每個格子拆成一或兩個三角形
    for i = 0:nDiv-1
        for j = 0:(nDiv-1-i)
            p1 = loc(i+1,j+1);
            p2 = loc(i+2,j+1);
            p3 = loc(i+1,j+2);
            tCount = tCount+1;
            Traw(tCount,:) = [p1 p2 p3];
            if i+j <= nDiv-2
                p4 = loc(i+2,j+2);
                tCount = tCount+1;
                Traw(tCount,:) = [p2 p4 p3];
            end
        end
    end
end

% 合併共用節點：先將座標取至固定精度，再統一節點編號。
scale = max(R,1);
Pkey = round((Praw/scale)*1e12)/1e12;
[~,first,oldToNew] = unique(Pkey,'rows','stable');
P = Praw(first,:);
T = oldToNew(Traw);

nP = size(P,1);
nT = size(T,1);

% 找出邊界節點；Dirichlet 邊界的位移為零。
boundaryMask = false(nP,1);
tol = 1e-9*max(R,1);
for s = 1:nSide
    A = V(s,:);
    B = V(mod(s,nSide)+1,:);
    AB = B-A;
    t = ((P-A)*AB')/(AB*AB');
    distance = abs((P(:,1)-A(1))*AB(2) - (P(:,2)-A(2))*AB(1))/norm(AB);
    boundaryMask = boundaryMask | (t>=-1e-9 & t<=1+1e-9 & distance<=tol);
end
freeNodes = find(~boundaryMask);
assert(numel(freeNodes)>nModes, '內部節點太少，請減小 h 或減少 nModes。');

%% ===== 3. 組裝 FEM 剛度矩陣與質量矩陣 =====
% 三節點線性三角形元素：
%   Ke(i,j) = ∫ grad(N_i)·grad(N_j) dΩ
%   Me(i,j) = ∫ N_i*N_j dΩ（採一致質量矩陣）
rows = zeros(9*nT,1);
cols = zeros(9*nT,1);
kval = zeros(9*nT,1);
mval = zeros(9*nT,1);
minArea = inf;
for e = 1:nT
    ind = T(e,:);
    x = P(ind,1);
    y = P(ind,2);
    areaE = abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))/2;
    minArea = min(minArea,areaE);
    assert(areaE>1e-14*R^2, '第 %d 個三角形面積過小，請檢查網格設定。',e);

    b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    Ke = (b*b' + c*c')/(4*areaE);
    Me = areaE/12 * [2 1 1; 1 2 1; 1 1 2];

    q = (e-1)*9+(1:9);
    [ii,jj] = ndgrid(ind,ind);
    rows(q) = ii(:);
    cols(q) = jj(:);
    kval(q) = Ke(:);
    mval(q) = Me(:);
end
K = sparse(rows,cols,kval,nP,nP);
M = sparse(rows,cols,mval,nP,nP);

%% ===== 4. 加入固定邊界條件並求解特徵值 =====
% 只保留內部自由節點，求解 Kff*u = lambda*Mff*u。
% 與 BIQM 比較時，應使用 k = sqrt(lambda)，而非 lambda。
Kff = K(freeNodes,freeNodes);
Mff = M(freeNodes,freeNodes);
opts.tol = 1e-10;
opts.maxit = 2000;
% 矩陣直接傳入 eigs，不必指定 issym 或 isreal。
[U,D] = eigs(Kff,Mff,nModes,'smallestabs',opts);
lambda = real(diag(D));
[lambda,order] = sort(lambda);
U = real(U(:,order));
assert(all(lambda>0), '求得非正特徵值，請檢查網格與邊界條件。');
k = sqrt(lambda);

% 把求得的內部節點位移補回完整網格，邊界位移設為零。
modeShape = zeros(nP,nModes);
modeShape(freeNodes,:) = U;
for j = 1:nModes
    amplitude = max(abs(modeShape(:,j)));
    assert(amplitude>0, '第 %d 個模態的振幅為零。',j);
    modeShape(:,j) = modeShape(:,j)/amplitude;
end

%% ===== 5. 顯示特徵值與網格資訊 =====
areaExact = nSide*R^2*sin(2*pi/nSide)/2;
areaMesh = sum(abs((P(T(:,2),1)-P(T(:,1),1)).*(P(T(:,3),2)-P(T(:,1),2)) ...
                 -(P(T(:,3),1)-P(T(:,1),1)).*(P(T(:,2),2)-P(T(:,1),2)))/2);
fprintf('正 %d 邊形：外接圓半徑 R = %.6f\n',nSide,R);
fprintf('網格：%d 個節點、%d 個三角形、%d 個內部自由節點\n',nP,nT,numel(freeNodes));
fprintf('分割數 nDiv = %d，最小元素面積 = %.6g\n',nDiv,minArea);
fprintf('幾何面積 = %.9f，網格面積 = %.9f\n',areaExact,areaMesh);
assert(abs(areaMesh-areaExact)<=1e-8*areaExact, '網格面積與幾何面積不符。');
fprintf('  模態序號          lambda = k^2              k\n');
for j = 1:nModes
    fprintf('    %3d          %15.9f       %15.9f\n',j,lambda(j),k(j));
end

Results = table((1:nModes)',lambda,k,'VariableNames',{'Mode','Lambda','k'});
% 如需另存特徵值表格，取消下一行最前面的 % 即可：
% writetable(Results,'FEM_RegularPolygon_Eigenvalues.csv');

%% ===== 6. 繪製三角形網格 =====
figure('Color','w','Name','FEM 三角形網格');
triplot(T,P(:,1),P(:,2),'k-');
axis equal tight;
xlabel('x'); ylabel('y');
title(sprintf('正 %d 邊形 FEM 網格（%d 個三角形）',nSide,nT));

%% ===== 7. 繪製各階振動模態 =====
nShow = min(plotModes,nModes);
if nShow > 0
    nCol = min(3,nShow);
    nRow = ceil(nShow/nCol);
    figure('Color','w','Name','FEM 固定邊界振動模態');
    for j = 1:nShow
        subplot(nRow,nCol,j);
        trisurf(T,P(:,1),P(:,2),modeShape(:,j),'EdgeColor','none');
        view(2); axis equal tight; shading interp;
        colormap(parula); caxis([-1,1]); colorbar;
        xlabel('x'); ylabel('y');
        title(sprintf('第 %d 模態，k = %.6f',j,k(j)));
    end
end

% 工作區常用變數：
% P：所有節點的二維座標；T：各三角形的三個節點編號
% K、M：整體剛度矩陣與質量矩陣；lambda：特徵值；k：特徵波數
% modeShape：各階模態；Results：特徵值結果表格

%% FEM_StarPolygon_Eigen_P2_LocalRefine.m
% 二維 Helmholtz 方程：星形 (N,k)-gon 的固定邊界特徵值問題
%   -∇²u = lambda*u  （域內）
%       u = 0         （邊界）
% 使用六節點二次三角形元素（P2 FEM），在星形凹角附近局部加密。
% 網格為符合式三角網格：相鄰元素共用邊中點，沒有懸掛節點。
% 這是指定凹角的局部加密，尚非依誤差估計器驅動的自適應 AMR。
% 注意：星形跨步參數 k 與特徵波數 k_eig 不同。
clear; clc; close all;

%% ===== 0. 使用者設定 =====
N         = 5;              % 星形外尖點數，例如 5、7、9
k         = 2;              % 星形跨步參數，例如 (5,2)、(7,2)、(7,3)
R         = 2.2270327288;   % 外尖點到中心距離（與 BIQM 相同）
h         = 0.2;           % 基礎網格尺寸；不必設成 0.001
nModes    = 12;             % 特徵波數數量（包含重根）
angle0    = pi/2;           % 第一個外尖點方向
plotModes = 6;              % 繪圖模態數，0 表示不繪圖
plotMesh  = true;           % 是否繪製 FEM 網格

% 局部加密設定：逐層縮小凹角附近的加密範圍
refineLevels = 4;           % 建議先 2~3 層
refineRadiusRatio = [0.30, 0.15, 0.075, 0.0375]; % 每層半徑/R

% 收斂測試：保持基礎 h 不變，比較局部加密層數（非完全巢狀的 h 比較）
doConvergence = false;
levelsList = [0 1 2 3];     % 由粗到細，0 表示不加密

% 記憶體保護：先限制單組三角形數量再組裝 P2 矩陣
maxElements = 10000000;

%% ===== 1. 檢查參數，建立星形邊界 =====
assert(isscalar(N) && isfinite(N) && N==floor(N) && N>=3, ...
    'N 必須是大於等於 3 的整數。');
assert(isscalar(k) && isfinite(k) && k==floor(k) && k>=1 && k<N/2, ...
    'k 必須是整數，且滿足 1 <= k < N/2。');
assert(gcd(N,k)==1, '(N,k) 必須互質。');
assert(isfinite(R) && isscalar(R) && R>0, 'R 必須為正數。');
assert(isfinite(h) && isscalar(h) && h>0, 'h 必須為正數。');
assert(nModes==floor(nModes) && nModes>=1, 'nModes 必須為正整數。');
assert(plotModes==floor(plotModes) && plotModes>=0, 'plotModes 必須為非負整數。');

% 與原有 Nk_gon_infield 幾何一致：外尖點與內凹點交替排列。
rInner = R*cos(k*pi/N)/cos((k-1)*pi/N);
outerTheta = angle0 + (0:N-1)'*2*pi/N;
innerTheta = outerTheta + pi/N;
V = zeros(2*N,2);
V(1:2:end,:) = R*[cos(outerTheta),sin(outerTheta)];
V(2:2:end,:) = rInner*[cos(innerTheta),sin(innerTheta)];
areaExact = polyarea(V(:,1),V(:,2));

assert(refineLevels>=0 && refineLevels==floor(refineLevels), ...
    'refineLevels 必須為非負整數。');
assert(numel(refineRadiusRatio)>=refineLevels, ...
    'refineRadiusRatio 的長度需大於等於 refineLevels。');
assert(all(refineRadiusRatio>0) && all(isfinite(refineRadiusRatio)), ...
    '加密半徑比例必須為正數。');
if doConvergence
    runLevels = unique(levelsList(:).','sorted');
    assert(all(runLevels>=0) && all(runLevels==floor(runLevels)), ...
        'levelsList 必須為非負整數。');
    assert(max(runLevels)<=numel(refineRadiusRatio), ...
        'levelsList 不得超過 refineRadiusRatio 的長度。');
else
    runLevels = refineLevels;
end

nRuns = numel(runLevels);
allK = nan(nRuns,nModes);
allElements = zeros(nRuns,1);
allDOF = zeros(nRuns,1);
allMinArea = zeros(nRuns,1);

%% ===== 2. 逐組網格進行 P2 FEM 特徵值計算 =====
for runID = 1:nRuns
    nLevel = runLevels(runID);
    fprintf('\n===== 星形 P2 FEM：h = %.4g，局部加密 = %d 層（%d/%d）=====\n', ...
        h,nLevel,runID,nRuns);

    % 初始網格沿星形中心分區；之後只在內凹頂點附近加密。
    edgeLen = sqrt(sum((V([2:end,1],:)-V).^2,2));
    nDiv = max(2,ceil(max(R,max(edgeLen))/h));
    estBase = 2*N*nDiv^2;
    assert(estBase<=maxElements, ...
        '初始網格已超過 maxElements；請增大 h。');
    [P,T] = makeStarMesh(V,nDiv,R);
    concavePoints = V(2:2:end,:);
    for lev = 1:nLevel
        radius = R*refineRadiusRatio(lev);
        [P,T] = refineAroundCorners(P,T,concavePoints,radius,maxElements);
        fprintf('  第 %d 層加密：半徑 = %.5g；三角形 = %d\n', ...
            lev,radius,size(T,1));
    end
    nP = size(P,1);
    nT = size(T,1);
    fprintf('初始 nDiv = %d；局部加密後三角形 = %d\n',nDiv,nT);

    % 把 P1 三角形升階成 P2：每條網格邊加入一個共用中點。
    % 六個局部節點順序：[頂點1,頂點2,頂點3,邊12,邊23,邊31]。
    [P2,T6,isBoundary] = upgradeToP2(P,T);
    freeNodes = find(~isBoundary);
    nDOF = numel(freeNodes);
    assert(nDOF>nModes,'內部自由度不足，請減小 h 或減少 nModes。');

    fprintf('分割數 nDiv = %d；三角形 = %d；P2 總節點 = %d；內部自由度 = %d\n', ...
        nDiv,nT,size(P2,1),nDOF);

    % 六點正權重三角形積分：對四次多項式精確，足以精確積分 P2 質量矩陣。
    % 參考三角形的重心座標 [L1,L2,L3]，權重總和為 1。
    a1 = 0.816847572980459; b1 = 0.091576213509771;
    a2 = 0.108103018168070; b2 = 0.445948490915965;
    Lq = [a1 b1 b1; b1 a1 b1; b1 b1 a1; ...
          a2 b2 b2; b2 a2 b2; b2 b2 a2];
    wq = [repmat(0.109951743655322,3,1); ...
          repmat(0.223381589678011,3,1)];

    % 預先計算 P2 形狀函數及其對重心座標的偏導數。
    Nq = zeros(6,6);
    dNq = zeros(6,6,3);
    for q = 1:6
        L1=Lq(q,1); L2=Lq(q,2); L3=Lq(q,3);
        Nq(:,q) = [L1*(2*L1-1); L2*(2*L2-1); L3*(2*L3-1); ...
                   4*L1*L2; 4*L2*L3; 4*L3*L1];
        dNq(:,q,1) = [4*L1-1;0;0;4*L2;0;4*L3];
        dNq(:,q,2) = [0;4*L2-1;0;4*L1;4*L3;0];
        dNq(:,q,3) = [0;0;4*L3-1;0;4*L2;4*L1];
    end

    % 先建立節點到內部自由度的對照，直接組裝自由度矩陣，節省記憶體。
    freeMap = zeros(size(P2,1),1);
    freeMap(freeNodes) = 1:nDOF;
    capacity = 36*nT;  % 記憶體需求隨元素數增加，勿任意提高 maxElements
    rowID = zeros(capacity,1);
    colID = zeros(capacity,1);
    stiffVal = zeros(capacity,1);
    massVal = zeros(capacity,1);
    ptr = 0;
    minArea = inf;
    areaMesh = 0;

    for e = 1:nT
        ind3 = T(e,:);
        X = P(ind3,:);
        J = [X(2,:)-X(1,:);X(3,:)-X(1,:)];
        twiceArea = abs(det(J));
        areaE = twiceArea/2;
        assert(areaE>1e-14*R^2,'第 %d 個三角形面積太小。',e);
        areaMesh = areaMesh+areaE;
        minArea = min(minArea,areaE);

        % L1、L2、L3 在實際元素中的梯度（與元素方向無關）。
        gL2 = (J\[1;0]).';
        gL3 = (J\[0;1]).';
        gradL = [-gL2-gL3;gL2;gL3];

        Ke = zeros(6,6);
        Me = zeros(6,6);
        for q = 1:6
            gradN = squeeze(dNq(:,q,:))*gradL; % 6x2
            Ke = Ke + areaE*wq(q)*(gradN*gradN.');
            Nv = Nq(:,q);
            Me = Me + areaE*wq(q)*(Nv*Nv.');
        end

        loc = freeMap(T6(e,:));
        active = find(loc>0);  % 直接排除固定邊界的自由度
        m = numel(active);
        [ii,jj] = ndgrid(loc(active),loc(active));
        pos = ptr+(1:m*m);
        rowID(pos) = ii(:);
        colID(pos) = jj(:);
        subK = Ke(active,active);
        subM = Me(active,active);
        stiffVal(pos) = subK(:);
        massVal(pos) = subM(:);
        ptr = ptr+m*m;
    end

    Kff = sparse(rowID(1:ptr),colID(1:ptr),stiffVal(1:ptr),nDOF,nDOF);
    Mff = sparse(rowID(1:ptr),colID(1:ptr),massVal(1:ptr),nDOF,nDOF);
    clear rowID colID stiffVal massVal;

    % 對稱化以消除浮點數累積造成的極小非對稱項。
    Kff = (Kff+Kff.')/2;
    Mff = (Mff+Mff.')/2;

    % shift=0 的 eigs 求解最接近零的廣義特徵值。
    opts.tol = 1e-10;
    opts.maxit = 2000;
    [U,D] = eigs(Kff,Mff,nModes,0,opts);
    lambda = real(diag(D));
    [lambda,order] = sort(lambda);
    U = real(U(:,order));
    assert(all(isfinite(lambda)) && all(lambda>0), ...
        '出現非正特徵值，請檢查網格、矩陣和邊界條件。');
    k_eig = sqrt(lambda);

    allK(runID,:) = k_eig(:).';
    allElements(runID) = nT;
    allMinArea(runID) = minArea;
    allDOF(runID) = nDOF;
    fprintf('幾何面積 = %.9f；網格面積 = %.9f；最小元素面積 = %.7g\n', ...
        areaExact,areaMesh,minArea);
    fprintf('  模態            lambda = k_eig^2                k_eig\n');
    for j=1:nModes
        fprintf('  %3d             %15.9f            %15.9f\n',j,lambda(j),k_eig(j));
    end

    % 保留最後（最細）一組網格的模態圖與工作區變數。
    if runID==nRuns
        Results = table((1:nModes)',lambda,k_eig, ...
            'VariableNames',{'Mode','Lambda','k_eig'});
        if plotMesh
            figure('Color','w','Name','星形 FEM P2 三角形網格');
            triplot(T,P(:,1),P(:,2),'k-');
            axis equal tight; xlabel('x'); ylabel('y');
            title(sprintf('(%d,%d)-gon P2 網格：%d 個三角形',N,k,nT));
        end
        nShow=min(plotModes,nModes);
        if nShow>0
            modeShape=zeros(size(P2,1),nShow);
            modeShape(freeNodes,:) = U(:,1:nShow);
            modeShape = modeShape./max(abs(modeShape),[],1);
            nCol=min(3,nShow);
            nRow=ceil(nShow/nCol);
            figure('Color','w','Name','星形 FEM P2 模態');
            % P2 模態圖不能只用三個角點插值，需使用四分小三角形繪製中點資訊。
            Tplot = [T6(:,[1 4 6]);T6(:,[4 2 5]); ...
                     T6(:,[6 5 3]);T6(:,[4 5 6])];
            for j=1:nShow
                subplot(nRow,nCol,j);
                trisurf(Tplot,P2(:,1),P2(:,2),modeShape(:,j), ...
                    'EdgeColor','none');
                view(2); axis equal tight; shading interp;
                colormap(parula); clim([-1,1]); colorbar;
                xlabel('x'); ylabel('y');
                title(sprintf('第 %d 模態，k = %.6f',j,k_eig(j)));
            end
        end
    end
end

%% ===== 3. 輸出加密層數收斂性比較 =====
if doConvergence
    Convergence = table(runLevels(:),allElements,allDOF,allMinArea, ...
        'VariableNames',{'Levels','Elements','FreeDOF','MinArea'});
    for j=1:nModes
        Convergence.(sprintf('k%d',j)) = allK(:,j);
    end
    disp('===== P2 FEM 局部加密層數收斂性比較 =====');
    disp(Convergence);
    if nRuns>1
        relativeChange = 100*abs(diff(allK,1,1))./abs(allK(2:end,:));
        disp('相鄰加密層數之相對變化（%），每欄對應一個模態：');
        disp(relativeChange);
    end
    % writetable(Convergence,sprintf('FEM_P2_LocalRefine_%d_%d.csv',N,k));
end

%% ===== 本檔案使用的局部函數 =====
function [P,T] = makeStarMesh(V,nDiv,R)
% 使用中心扇形的規則三角形細分，沿共用邊合併重複節點。
    nSector=size(V,1);
    nLocal=(nDiv+1)*(nDiv+2)/2;
    Praw=zeros(nSector*nLocal,2);
    Traw=zeros(nSector*nDiv^2,3);
    pCount=0; tCount=0;
    for s=1:nSector
        A=V(s,:); B=V(mod(s,nSector)+1,:);
        loc=zeros(nDiv+1,nDiv+1);
        for i=0:nDiv
            for j=0:nDiv-i
                pCount=pCount+1;
                Praw(pCount,:)=(i*A+j*B)/nDiv;
                loc(i+1,j+1)=pCount;
            end
        end
        for i=0:nDiv-1
            for j=0:nDiv-1-i
                p1=loc(i+1,j+1);
                p2=loc(i+2,j+1);
                p3=loc(i+1,j+2);
                tCount=tCount+1;
                Traw(tCount,:)=[p1 p2 p3];
                if i+j<=nDiv-2
                    p4=loc(i+2,j+2);
                    tCount=tCount+1;
                    Traw(tCount,:)=[p2 p4 p3];
                end
            end
        end
    end
    Pkey=round((Praw/max(R,1))*1e12)/1e12;
    [~,first,oldToNew]=unique(Pkey,'rows','stable');
    P=Praw(first,:);
    T=oldToNew(Traw);
end

function [P2,T6,isBoundary] = upgradeToP2(P,T)
% 三節點 -> 六節點；相鄰三角形共用同一個邊中點節點。
% 透過出現一次的網格邊自動識別外邊界，包含所有 P2 邊界中點。
    nT=size(T,1);
    allEdges=[T(:,[1 2]);T(:,[2 3]);T(:,[3 1])];
    [uniqueEdges,~,idx]=unique(sort(allEdges,2),'rows');
    nE=size(uniqueEdges,1);
    edgeID=reshape(idx,nT,3);
    counts=accumarray(idx,1,[nE,1]);
    boundaryEdges=uniqueEdges(counts==1,:);
    isBoundary=false(size(P,1)+nE,1);
    isBoundary(unique(boundaryEdges(:)))=true;
    isBoundary(size(P,1)+find(counts==1))=true;
    P2=[P;(P(uniqueEdges(:,1),:)+P(uniqueEdges(:,2),:))/2];
    T6=[T,size(P,1)+edgeID];
end

function [Pnew,Tnew] = refineAroundCorners(P,T,corners,radius,maxElements)
% 在凹角附近進行符合式局部加密。
% 1) 距離凹角小於 radius 的三角形全部標記三條邊。
% 2) 全域建立唯一邊中點，鄰接三角形共用中點。
% 3) 依每個三角形被標記的邊數，採 1->2、2->3、3->4 拆分。
%    因此不產生懸掛節點；幾何邊界仍保持原本直線段。
    nt=size(T,1);
    choose=false(nt,1);
    for c=1:size(corners,1)
        v=corners(c,:);
        d2=(P(:,1)-v(1)).^2+(P(:,2)-v(2)).^2;
        choose=choose | any(reshape(d2(T(:)),nt,3)<radius^2,2);
    end
    if ~any(choose)
        Pnew=P; Tnew=T; return;
    end
    Edges=[T(:,[1 2]);T(:,[2 3]);T(:,[3 1])];
    [E,~,ic]=unique(sort(Edges,2),'rows');
    eID=reshape(ic,nt,3);
    edgeMarked=false(size(E,1),1);
    edgeMarked(eID(choose,:))=true;
    midID=zeros(size(E,1),1);
    selected=find(edgeMarked);
    midID(selected)=size(P,1)+(1:numel(selected));
    Pnew=[P;(P(E(selected,1),:)+P(E(selected,2),:))/2];
    markedLocal=midID(eID);  % nt x 3；0 表示未分割
    nMark=sum(markedLocal>0,2);
    expected=nt+sum(nMark==1)+2*sum(nMark==2)+3*sum(nMark==3);
    assert(expected<=maxElements, ...
        ['局部加密後預估有 %d 個三角形，超過 maxElements=%d。' ...
         '請降低 refineLevels、縮小 refineRadiusRatio 或增大 h。'], ...
        expected,maxElements);
    Tnew=zeros(expected,3);
    ptr=0;
    for it=1:nt
        v=T(it,:);
        mid=markedLocal(it,:);
        active=find(mid>0);
        switch numel(active)
            case 0
                newTri=v;
            case 1
                % 旋轉頂點，使被分割的邊為 (a,b)。
                a=v(active); b=v(mod(active,3)+1); c=v(mod(active+1,3)+1);
                m=mid(active);
                newTri=[a m c; m b c];
            case 2
                % 找到連續兩條被分割的邊 (a,b)、(b,c)。
                s=active(1);
                if mid(mod(s,3)+1)==0
                    s=active(2);
                end
                a=v(s); b=v(mod(s,3)+1); c=v(mod(s+1,3)+1);
                mab=mid(s); mbc=mid(mod(s,3)+1);
                newTri=[b mbc mab; a mab c; mab mbc c];
            case 3
                a=v(1); b=v(2); c=v(3);
                mab=mid(1); mbc=mid(2); mca=mid(3);
                newTri=[a mab mca; mab b mbc; mca mbc c; mab mbc mca];
        end
        nr=size(newTri,1);
        Tnew(ptr+(1:nr),:)=newTri;
        ptr=ptr+nr;
    end
    assert(ptr==expected,'內部錯誤：三角形數量不符。');
end

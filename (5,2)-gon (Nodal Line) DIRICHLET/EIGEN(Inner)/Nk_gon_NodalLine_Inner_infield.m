%% ==========================================================
% 數值節線切割：只建立節線內部子域 F15.DAT
%
% 讀取 xspmoU2.dat: [x, y, u]，由 u=0 等值線建立內孔。
% 僅使用實際計算的數值節線；不進行圓形擬合。
%
% 注意：此節線是由離散位移場插值得到的數值近似，並非解析精確節線。
%% ==========================================================
clear; clc; close all;

%% ===== 0. 使用者設定 =====
outputFile = 'F15.DAT';                 % 直接輸出到目前工作資料夾
modeFile   = 'xspmoU2.dat';             % 三欄：x, y, u

% 外邊界：與原始 (N,k)-gon 程式完全相同
Nout     = 5;
kout     = 2;
Rout     = 2.2270327288;
nNodeOut = 200;                       % 必須為 2*Nout 的倍數

% 內邊界：沿數值零位移曲線等弧長佈點
nNodeIn = 40;
centerX = 0.0;
centerY = 0.0;
theta0Out = pi/2;

% 用於辨識中心節線的參考點（需落在內側節線圍成的區域內）
seedX = centerX;
seedY = centerY;

% 以內部位移資料插值後擷取 u=0
nGrid = 501;                          % 等值線擷取網格，越大越細
interpMethod = 'linear';               % scatteredInterpolant 方法：'linear' 或 'natural'

% 安全檢查；零位移資料精度有限，不要求所有重取樣點的 u 完全等於 0
maxAllowedResidualRatio = 0.12;        % max |u_on_cut| / max |u_input| 上限
showAllZeroContours = true;            % 圖上顯示其他找到的零等值線，供檢查

%% ===== 1. 參數與資料檢查 =====
check_nkgon_parameters(Nout,kout,Rout,nNodeOut,'外邊界');
assert(nNodeIn >= 3 && nNodeIn == floor(nNodeIn), ...
    'nNodeIn 必須為大於等於 3 的整數。');
assert(nGrid >= 101 && nGrid == floor(nGrid), ...
    'nGrid 必須為至少 101 的整數。');
if ~isfile(modeFile)
    uploadedName = 'xspmoU2(20261117-051939).dat';
    if isfile(uploadedName)
        modeFile = uploadedName;      % 支援本次上傳時的原始檔名
    else
        error('找不到模態資料檔：%s。請將它放在目前工作資料夾。',modeFile);
    end
end

A = readmatrix(modeFile);
assert(size(A,2)>=3,'模態資料至少要包含三欄：x y u。');
A = A(:,1:3);
A = A(all(isfinite(A),2),:);
assert(size(A,1)>=20,'有效模態資料不足。');

% 去除重複座標（若有相同座標，取位移平均值）
[xy,~,group] = unique(A(:,1:2),'rows');
u = accumarray(group,A(:,3),[],@mean);
assert(min(u)<0 && max(u)>0, ...
    '位移場必須含有正負值，才可尋找內部零位移節線。');

%% ===== 2. 建立五角星外邊界（沿用原程式） =====
[xOut,yOut] = make_nkgon_boundary( ...
    Nout,kout,Rout,nNodeOut,centerX,centerY,theta0Out);
xOut = xOut(:).'; yOut = yOut(:).';
if signed_polygon_area(xOut,yOut)<0
    xOut=fliplr(xOut); yOut=fliplr(yOut);
end

%% ===== 3. 以離散位移場擷取 u=0 節線 =====
F = scatteredInterpolant(xy(:,1),xy(:,2),u,interpMethod,'none');
xx = linspace(min(xy(:,1)),max(xy(:,1)),nGrid);
yy = linspace(min(xy(:,2)),max(xy(:,2)),nGrid);
[X,Y] = meshgrid(xx,yy);
U = F(X,Y);

% 外域及資料凸包外不採用；不使用 extrapolation 生成虛假節線
insideStar = inpolygon(X,Y,xOut,yOut);
U(~insideStar) = NaN;

% contourc 返回多段曲線，以 header [level; number_of_points] 分段
C = contourc(xx,yy,U,[0 0]);
curves = parse_contour_matrix(C);
assert(~isempty(curves), ...
    '沒有擷取到 u=0 等值線；請檢查 modeFile、位移正負號或 nGrid。');

% 尋找包住 seed 的「封閉」節線，並避免把外輪廓誤當內孔
isEligible = false(numel(curves),1);
areas = inf(numel(curves),1);
gridSize = max((xx(end)-xx(1))/(nGrid-1),(yy(end)-yy(1))/(nGrid-1));
for t=1:numel(curves)
    q=curves{t};
    if size(q,1)<10; continue; end
    if norm(q(1,:)-q(end,:)) > 1.5*gridSize
        continue;                       % 排除開放等值線
    end
    q(end,:) = q(1,:);                    % 僅對已幾乎閉合的曲線確保嚴格閉合
    [in,on] = inpolygon(seedX,seedY,q(:,1),q(:,2));
    if ~(in || on); continue; end
    if ~all(inpolygon(q(:,1),q(:,2),xOut,yOut)); continue; end
    a = abs(signed_polygon_area(q(:,1).',q(:,2).'));
    if a < 50*gridSize^2; continue; end
    isEligible(t)=true;
    areas(t)=a;
    curves{t}=q;
end
assert(any(isEligible), ...
    ['找不到包住參考點的封閉內部節線。請先檢查資料繪圖、' ...
     '增加 nGrid，或調整 seedX / seedY；程式不會自動以圓形取代。']);

% 若有多條同心封閉節線，優先取最小包圍面積的那條
[~,pick]=min(areas);
q=curves{pick};
q=remove_duplicate_neighbors(q);
assert(size(q,1)>=10,'選到的節線有效節點太少。');

%% ===== 4. 封閉節線依弧長等距重新佈點 =====
[xCut,yCut] = resample_closed_curve(q,nNodeIn);
% 數值節線先統一成 CCW，適用於獨立的內部子域
if signed_polygon_area(xCut,yCut)<0
    xCut=fliplr(xCut); yCut=fliplr(yCut);
end

% 幾何檢查：節線位於五角星內部且沒有自交
assert(all(inpolygon(xCut,yCut,xOut,yOut)), ...
    '部分節線節點落在五角星外，不輸出 F15.DAT。');
assert(~polyline_intersects([xOut(:);xOut(1)],[yOut(:);yOut(1)], ...
    [xCut(:);xCut(1)],[yCut(:);yCut(1)]), ...
    '節線與外邊界相交，不輸出 F15.DAT。');
assert(~polyline_self_intersects([xCut(:);xCut(1)],[yCut(:);yCut(1)]), ...
    '節線重新佈點後自交，不輸出 F15.DAT。');
assert(inpolygon(seedX,seedY,xCut,yCut), ...
    '重新佈點的內側區域沒有包住指定 seed，請檢查節線。');

% 節線位移殘差（只做品質檢查）
uOnCut=F(xCut(:),yCut(:));
assert(all(isfinite(uOnCut)), ...
    '有節線節點落在插值資料支持範圍外。');
residualRatio=max(abs(uOnCut))/max(abs(u));
if residualRatio>maxAllowedResidualRatio
    warning('節線上最大位移殘差比 %.4g 偏大，請檢查網格及插值結果。',residualRatio);
end

%% ===== 5. 僅輸出節線內部子域的 F15.DAT =====
% 內側子域：只有節線形成的封閉邊界，逆時針 CCW
% 外五角星只用於辨識與檢查節線，不寫入 F15.DAT
writeF15(outputFile,xCut(:),yCut(:),nNodeIn);
fprintf('已輸出節線內部子域：%s（%d 元素，CCW）\n',fullfile(pwd,outputFile),nNodeIn);

%% ===== 6. 邊界佈點與模態預覽 =====
figure('Color','w','Name','Numerical nodal line / F15 boundary'); hold on;
if showAllZeroContours
    for t=1:numel(curves)
        v=curves{t};
        plot(v(:,1),v(:,2),'-','Color',[.78 .78 .78], ...
            'HandleVisibility','off');
    end
end
h1=plot([xOut,xOut(1)],[yOut,yOut(1)],'-o', ...
    'LineWidth',1.2,'MarkerSize',3);
h2=plot(q(:,1),q(:,2),'--', ...
    'Color',[.85 .35 .08],'LineWidth',1.5);
h3=plot([xCut,xCut(1)],[yCut,yCut(1)],'-o', ...
    'LineWidth',1.2,'MarkerSize',4);
plot(seedX,seedY,'k+','MarkerSize',7,'HandleVisibility','off');
axis equal; box on; grid on;
xlabel('x'); ylabel('y');
title(sprintf('(%d,%d)-gon: numerical nodal line / %d-point cut', ...
    Nout,kout,nNodeIn));
legend([h1,h2,h3],{'Original star boundary (CCW)', ...
    'Extracted u=0 nodal line','Resampled nodal line (CCW)'},'Location','best');

figure('Color','w','Name','Original modal field and selected nodal line'); hold on;
contourf(X,Y,U,30,'LineStyle','none'); colorbar;
plot(q(:,1),q(:,2),'k-','LineWidth',2);
plot([xOut,xOut(1)],[yOut,yOut(1)],'k-','LineWidth',1.3);
axis equal; box on; xlabel('x'); ylabel('y');
title('Numerical mode (u) and extracted u=0 line');

fprintf('============================================\n');
fprintf('模態資料：%s（%d 個不同場點）\n',modeFile,size(xy,1));
fprintf('選中第 %d 條封閉節線；內部面積約 %.7f\n',pick,areas(pick));
fprintf('原始五角星：%d 元素（僅作節線定位參考，不輸出）\n',nNodeOut);
fprintf('內部子域節線：%d 元素，CCW\n',nNodeIn);
fprintf('節線上最大相對位移殘差：%.6g\n',residualRatio);
fprintf('============================================\n');

%% ===== 以下為本程式使用的函式 =====
function curves=parse_contour_matrix(C)
    curves={};
    idx=1;
    while idx<size(C,2)
        n=round(C(2,idx));
        if n<=0 || idx+n>size(C,2)
            error('contourc 產生的等值線格式不正確。');
        end
        q=C(:,idx+(1:n)).';
        curves{end+1,1}=q; %#ok<AGROW>
        idx=idx+n+1;
    end
end

function q=remove_duplicate_neighbors(q)
    if norm(q(1,:)-q(end,:))>1e-12
        q(end+1,:)=q(1,:);
    end
    keep=[true; hypot(diff(q(:,1)),diff(q(:,2)))>1e-11];
    q=q(keep,:);
    if norm(q(1,:)-q(end,:))>1e-12
        q(end+1,:)=q(1,:);
    end
end

function [x,y]=resample_closed_curve(q,n)
    % q 為閉合曲線，最後一點與第一點重合
    d=hypot(diff(q(:,1)),diff(q(:,2)));
    assert(all(d>0),'節線含重複相鄰點。');
    s=[0;cumsum(d)];
    total=s(end);
    assert(total>0,'節線周長必須大於零。');
    target=(0:n-1)'*total/n;
    x=interp1(s,q(:,1),target,'linear').';
    y=interp1(s,q(:,2),target,'linear').';
end

function tf=polyline_intersects(x1,y1,x2,y2)
    tf=false;
    for i=1:numel(x1)-1
        p=[x1(i),y1(i)]; r=[x1(i+1),y1(i+1)]-p;
        for j=1:numel(x2)-1
            q=[x2(j),y2(j)]; s=[x2(j+1),y2(j+1)]-q;
            if segments_intersect(p,r,q,s)
                tf=true; return;
            end
        end
    end
end

function tf=polyline_self_intersects(x,y)
    tf=false; n=numel(x)-1;
    for i=1:n
        p=[x(i),y(i)]; r=[x(i+1),y(i+1)]-p;
        for j=i+2:n
            if i==1 && j==n; continue; end % 相鄰閉合段共享端點
            q=[x(j),y(j)]; s=[x(j+1),y(j+1)]-q;
            if segments_intersect(p,r,q,s)
                tf=true; return;
            end
        end
    end
end

function tf=segments_intersect(p,r,q,s)
    cross2=@(a,b) a(1)*b(2)-a(2)*b(1);
    den=cross2(r,s);
    if abs(den)<1e-13
        tf=false; return;               % 近乎平行：當前資料預期不會重疊
    end
    t=cross2(q-p,s)/den;
    u=cross2(q-p,r)/den;
    tf=t>=0 && t<=1 && u>=0 && u<=1;
end

function [x,y] = make_nkgon_boundary( ...
    N,k,R,nNode,centerX,centerY,theta0)

    nPerSeg = nNode/(2*N);

    %% 外尖點 P
    theta = theta0 + 2*pi*(0:N-1)/N;

    Px = centerX + R*cos(theta);
    Py = centerY + R*sin(theta);

    %% 凹點 Q
    Qx = zeros(1,N);
    Qy = zeros(1,N);

    if k == 1

        for i = 1:N
            iNext = mod(i,N) + 1;

            Qx(i) = 0.5*(Px(i)+Px(iNext));
            Qy(i) = 0.5*(Py(i)+Py(iNext));
        end

    else

        for i = 1:N

            i1 = i;
            i2 = mod(i-1+k,N) + 1;
            i3 = mod(i,N) + 1;
            i4 = mod(i-k+N,N) + 1;

            x1 = Px(i1);
            y1 = Py(i1);

            dx1 = Px(i2)-Px(i1);
            dy1 = Py(i2)-Py(i1);

            x3 = Px(i3);
            y3 = Py(i3);

            dx2 = Px(i4)-Px(i3);
            dy2 = Py(i4)-Py(i3);

            den = dx1*dy2-dy1*dx2;

            if abs(den) < 1.0e-14
                error('建立 (%d,%d)-gon 凹點時兩直線平行。',N,k);
            end

            s = ((x3-x1)*dy2-(y3-y1)*dx2)/den;

            Qx(i) = x1+s*dx1;
            Qy(i) = y1+s*dy1;
        end
    end

    %% 邊界控制點 P,Q,P,Q,...
    Bx = zeros(1,2*N);
    By = zeros(1,2*N);

    for i = 1:N
        Bx(2*i-1) = Px(i);
        By(2*i-1) = Py(i);

        Bx(2*i) = Qx(i);
        By(2*i) = Qy(i);
    end

    %% 沿 2N 條邊均勻佈點
    x = zeros(1,nNode);
    y = zeros(1,nNode);

    idx = 0;

    for i = 1:2*N

        iNext = mod(i,2*N)+1;

        for j = 0:nPerSeg-1

            idx = idx+1;
            t = j/nPerSeg;

            x(idx) = (1-t)*Bx(i)+t*Bx(iNext);
            y(idx) = (1-t)*By(i)+t*By(iNext);
        end
    end
end

%% ==========================================================
% 參數檢查
%% ==========================================================
function check_nkgon_parameters(N,k,R,nNode,nameText)

    if N < 3 || N ~= floor(N)
        error('%s：N 必須為大於等於 3 的整數。',nameText);
    end

    if k < 1 || k ~= floor(k)
        error('%s：k 必須為大於等於 1 的整數。',nameText);
    end

    if k > floor((N-1)/2)
        error('%s：需滿足 k <= floor((N-1)/2)。',nameText);
    end

    if R <= 0
        error('%s：R 必須大於 0。',nameText);
    end

    if mod(nNode,2*N) ~= 0
        error('%s：nNode 必須可以被 2*N 整除。',nameText);
    end
end

%% ==========================================================
% 有號面積：>0 CCW，<0 CW
%% ==========================================================
function A = signed_polygon_area(x,y)

    x2 = [x(2:end),x(1)];
    y2 = [y(2:end),y(1)];

    A = 0.5*sum(x.*y2-x2.*y);
end

%% ==========================================================
% 輸出 F15.DAT：每條邊界各自閉合，不能誤接不同邊界
%% ==========================================================
function writeF15(filename,x,y,loopSizes)
    loopSizes=loopSizes(:).';
    nTotal=sum(loopSizes);
    assert(numel(x)==nTotal && numel(y)==nTotal, ...
        'F15.DAT 節點總數與各邊界元素數不一致。');
    assert(all(loopSizes>=3),'每條閉合邊界至少需要 3 個節點。');
    fid=fopen(filename,'wt');
    assert(fid>0,'無法建立 %s',filename);
    closer=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%6d\n%6d\n',-1,15);
    for i=1:nTotal
        fprintf(fid,'%10d%10d%10d%10d%13.5E%13.5E%13.5E\n', ...
            i,0,0,11,x(i),y(i),0);
    end
    fprintf(fid,'%6d\n%6d\n%6d\n',-1,-1,71);
    offset=0;
    for t=1:numel(loopSizes)
        n=loopSizes(t);
        for i=1:n
            id=offset+i;
            next=offset+mod(i,n)+1;
            fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n', ...
                id,1,21,1,1,7,2);
            fprintf(fid,'%10d%10d\n',id,next);
        end
        offset=offset+n;
    end
    fprintf(fid,'%6d\n',-1);
end

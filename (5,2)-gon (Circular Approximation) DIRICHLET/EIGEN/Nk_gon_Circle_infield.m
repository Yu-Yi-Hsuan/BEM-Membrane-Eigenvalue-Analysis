%% ==========================================================
% 雙聯通邊界 F15.DAT 產生程式
%
% 外邊界：一組 (N,k)-gon
% 內邊界：圓形孔洞
%
% 可分別指定：
%   外邊界：N、k、外接圓半徑 Rout、節點數 nNodeOut
%   內圓：半徑 Rin、節點數 nNodeIn
%
% 方向：
%   外邊界：逆時針 CCW
%   內圓：順時針 CW
%
% MATLAB R2024a
%% ==========================================================

clear; clc; close all;

%% ===== 使用者設定 =====
outputFile = 'F15.DAT';

% ---------- 外邊界：(N,k)-gon ----------
Nout     = 5;
kout     = 2;
Rout     = 2.2270327288;
nNodeOut = 200;

% ---------- 內邊界：圓形 ----------
Rin      = 0.4938;     % 內圓半徑
nNodeIn  = 40;         % 內圓節點數

% ---------- 共用設定 ----------
centerX = 0.0;
centerY = 0.0;

% 外 (N,k)-gon 第一個尖點方向
theta0Out = pi/2;

%% ===== 參數檢查 =====
check_nkgon_parameters(Nout,kout,Rout,nNodeOut,'外邊界');

if Rin <= 0
    error('內圓半徑 Rin 必須大於 0。');
end

if nNodeIn < 3 || nNodeIn ~= floor(nNodeIn)
    error('內圓節點數 nNodeIn 必須為大於等於 3 的整數。');
end

%% ===== 建立外邊界 =====
[xOut,yOut] = make_nkgon_boundary( ...
    Nout,kout,Rout,nNodeOut,centerX,centerY,theta0Out);

% 外邊界維持逆時針
if signed_polygon_area(xOut,yOut) < 0
    xOut = fliplr(xOut);
    yOut = fliplr(yOut);
end

%% ===== 建立內圓邊界 =====
% 為了讓內孔法向量正確，內邊界採順時針 CW
thetaIn = linspace(0,-2*pi,nNodeIn+1);
thetaIn(end) = [];

xIn = centerX + Rin*cos(thetaIn);
yIn = centerY + Rin*sin(thetaIn);

%% ===== 確認內圓完全位於外邊界內 =====
insideCheck = inpolygon(xIn,yIn,xOut,yOut);

if ~all(insideCheck)
    warning(['部分內圓節點不在外 (N,k)-gon 內。' ...
             '請減小 Rin 或調整外邊界幾何。']);
end

%% ===== 整合節點 =====
nTotal = nNodeOut + nNodeIn;

xAll = [xOut(:); xIn(:)];
yAll = [yOut(:); yIn(:)];

%% ===== 輸出 F15.DAT =====
fid = fopen(outputFile,'w');

if fid == -1
    error('無法建立 %s。',outputFile);
end

% ---------- 節點資料 ----------
fprintf(fid,'%6d\n',-1);
fprintf(fid,'%6d\n',15);

for i = 1:nTotal
    fprintf(fid,'%10d%10d%10d%10d%13.5E%13.5E%13.5E\n', ...
        i,0,0,11,xAll(i),yAll(i),0.0);
end

% ---------- 元素資料 ----------
fprintf(fid,'%6d\n',-1);
fprintf(fid,'%6d\n',-1);
fprintf(fid,'%6d\n',71);

elemID = 0;

% 外邊界元素
for i = 1:nNodeOut

    elemID = elemID + 1;

    if i < nNodeOut
        node1 = i;
        node2 = i + 1;
    else
        node1 = nNodeOut;
        node2 = 1;
    end

    fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n', ...
        elemID,1,21,1,1,7,2);

    fprintf(fid,'%10d%10d\n',node1,node2);
end

% 內圓元素
offset = nNodeOut;

for i = 1:nNodeIn

    elemID = elemID + 1;

    if i < nNodeIn
        node1 = offset + i;
        node2 = offset + i + 1;
    else
        node1 = offset + nNodeIn;
        node2 = offset + 1;
    end

    fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n', ...
        elemID,1,21,1,1,7,2);

    fprintf(fid,'%10d%10d\n',node1,node2);
end

fprintf(fid,'%6d\n',-1);
fclose(fid);

%% ===== 繪製邊界佈點圖 =====
figure;
hold on;

hOut = plot([xOut,xOut(1)],[yOut,yOut(1)],'-o', ...
    'MarkerSize',3,'LineWidth',1.2);

hIn = plot([xIn,xIn(1)],[yIn,yIn(1)],'-o', ...
    'MarkerSize',3,'LineWidth',1.2);

axis equal;
box on;
grid on;

xlabel('x');
ylabel('y');

title(sprintf( ...
    'Double-connected domain: outer (%d,%d)-gon, inner circle', ...
    Nout,kout));

legend([hOut,hIn], ...
    'Outer boundary (CCW)', ...
    'Inner circle (CW)', ...
    'Location','best');

allX = [xOut,xIn];
allY = [yOut,yIn];

margin = 0.15*max([max(allX)-min(allX), max(allY)-min(allY)]);

xlim([min(allX)-margin,max(allX)+margin]);
ylim([min(allY)-margin,max(allY)+margin]);

%% ===== 顯示結果 =====
fprintf('============================================\n');
fprintf('%s 已成功建立。\n',outputFile);

fprintf('\n外邊界：\n');
fprintf('  (N,k)-gon = (%d,%d)\n',Nout,kout);
fprintf('  Rout      = %.10f\n',Rout);
fprintf('  節點數    = %d\n',nNodeOut);
fprintf('  方向      = 逆時針 CCW\n');

fprintf('\n內邊界：\n');
fprintf('  圓半徑 Rin = %.10f\n',Rin);
fprintf('  節點數     = %d\n',nNodeIn);
fprintf('  方向       = 順時針 CW\n');

fprintf('\n總節點數 = %d\n',nTotal);
fprintf('總元素數 = %d\n',elemID);
fprintf('============================================\n');

%% ==========================================================
% 建立單一 (N,k)-gon 邊界
%% ==========================================================
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

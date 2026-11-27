%% ==========================================================
% L 型薄膜邊界
% MATLAB R2024a
%% ==========================================================
clear; clc; close all;

%% ===== 使用者設定 =====
outputFile = 'F15.DAT';
nNode = 160;       % 總邊界節點數等於常數邊界元素總數
                   % 必須能被 8 整除（周長為 8）

%% ===== 文獻幾何設定（為便於比較，請勿縮放） =====
% 邊界逆時針排列，各邊長依序為 2、1、1、1、1、2
V = [-1 -1;
      1 -1;
      1  0;
      0  0;
      0  1;
     -1  1];

nEdges = size(V,1);
Vnext = V([2:end,1],:);
edgeLengths = hypot(Vnext(:,1)-V(:,1),Vnext(:,2)-V(:,2));
perimeter = sum(edgeLengths);
areaDomain = polyarea(V(:,1),V(:,2));

assert(isscalar(nNode) && isfinite(nNode) && nNode>=8 && ...
       nNode==fix(nNode) && mod(nNode,8)==0, ...
       'nNode 必須為大於等於 8 且可被 8 整除的整數。');
assert(abs(perimeter-8)<1e-12 && abs(areaDomain-3)<1e-12, ...
       'L 型幾何與 Amore（2008）的設定不一致。');

% 所有直線邊均採用相同的實際元素長度。
% 長度為 2 的邊，其元素數為長度為 1 的邊的兩倍。
nPerUnit = nNode/8;
nElem = round(edgeLengths*nPerUnit);
assert(sum(nElem)==nNode, '邊界元素數分配錯誤。');

%% ===== 產生邊界節點（不重複各邊終點） =====
x = zeros(nNode,1);
y = zeros(nNode,1);
idx = 0;
for s = 1:nEdges
    A = V(s,:);
    B = Vnext(s,:);
    m = nElem(s);
    for j = 0:m-1
        idx = idx+1;
        t = j/m;
        point = (1-t)*A+t*B;
        x(idx) = point(1);
        y(idx) = point(2);
    end
end
assert(idx==nNode);

segmentLength = hypot(diff([x;x(1)]), diff([y;y(1)]));
assert(all(segmentLength>0),'偵測到長度為零的元素。');
assert(max(segmentLength)-min(segmentLength)<1e-10, ...
       '各邊界元素長度不一致。');

%% ===== 輸出 F15.DAT（維持原正方形程式格式） =====
fid = fopen(outputFile,'w');
if fid==-1
    error('無法建立 %s。',outputFile);
end
fprintf(fid,'%6d\n',-1);
fprintf(fid,'%6d\n',15);
for i=1:nNode
    fprintf(fid,'%10d%10d%10d%10d%13.5E%13.5E%13.5E\n', ...
        i,0,0,11,x(i),y(i),0.0);
end
fprintf(fid,'%6d\n',-1);
fprintf(fid,'%6d\n',-1);
fprintf(fid,'%6d\n',71);
for i=1:nNode
    nextNode = mod(i,nNode)+1;
    fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n', ...
        i,1,21,1,1,7,2);
    fprintf(fid,'%10d%10d\n',i,nextNode);
end
fprintf(fid,'%6d\n',-1);
fclose(fid);

%% ===== 繪製 L 型邊界佈點圖 =====
figure('Color','w');
plot([x;x(1)],[y;y(1)],'-o','LineWidth',1, ...
     'MarkerSize',3.5);
hold on;
plot(V(:,1),V(:,2),'ks','MarkerFaceColor','k','MarkerSize',5);
axis equal;
xlim([-1.2,1.2]); ylim([-1.2,1.2]);
grid on; box on;
xlabel('x'); ylabel('y');
title(sprintf('Amore（2008）L 型邊界佈點圖：%d 個元素',nNode));

%% ===== 顯示計算資訊 =====
fprintf('%s 已成功建立。\n',outputFile);
fprintf('幾何：[-1,1] × [-1,1] 移除右上角單位正方形。\n');
fprintf('面積 = %.6f；周長 = %.6f\n',areaDomain,perimeter);
fprintf('各邊元素數 = %s\n',mat2str(nElem.'));
fprintf('總邊界元素數 = %d\n',nNode);
fprintf('均勻元素長度 = %.8f\n',segmentLength(1));
fprintf('邊界方向：逆時針；內凹角位置：(0,0)。\n');

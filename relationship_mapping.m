% =========================================================================
% 相机像素 ↔ 3D世界坐标 ↔ 投影仪列号 映射关系建立
%
% 流程:
%   1. 加载相机内参 K
%   2. 检测标定板图像上的棋盘格角点 → 计算外参 R, t
%   3. 射线-平面求交: 相机像素 (u,v) → 标定板上 3D 坐标 (X,Y,0)
%   4. 加载绝对相位 → 投影仪列号
%   5. 输出映射表: [u, v, X, Y, Z, projector_col]
%
% 参考文献: Zhang (2000), 结构光投影仪标定方法
% =========================================================================
clear; close all; clc;

%% ======================== 参数设置 ========================
% 棋盘格参数（必须与 neican.m 中的标定板一致）
squareSize = 25;        % 棋盘格方格边长（mm）
numRows    = 8;         % 棋盘格内角点行数
numCols    = 11;        % 棋盘格内角点列数

% 图像文件路径
boardImageFile    = 'calib_board.png';      % 标定板图像（用于检测角点）
projectedImageFile = 'projected_on_board.png'; % 投影图案到标定板上的图像（用于可视化）

% 无效像素过滤阈值
minReliability = 0.05;   % 相位可靠性阈值（B/A），低于此值视为不可靠

%% ======================== 加载相机内参 ========================
fprintf('加载相机内参...\n');
if exist('camera_intrinsics.mat', 'file')
    Kdata = load('camera_intrinsics.mat');
    K = Kdata.K;
    if isfield(Kdata, 'cameraParams')
        cameraParams = Kdata.cameraParams;   % 含畸变系数，用于去畸变
    else
        cameraParams = [];                    % 旧标定文件，跳过畸变校正
    end
else
    error('未找到 camera_intrinsics.mat，请先运行 neican.m');
end
fprintf('  fx = %.2f, fy = %.2f, cx = %.2f, cy = %.2f\n', K(1,1), K(2,2), K(1,3), K(2,3));

%% ======================== 加载绝对相位结果 ========================
fprintf('加载绝对相位结果...\n');
phaseData = load('absolute_phase_result.mat');
absColumn  = phaseData.absColumn;    % 亚像素投影仪列号
validMask  = phaseData.validMask;    % 有效像素掩码（来自格雷码）
reliability = phaseData.reliability; % 相位可靠性
[imgH, imgW] = size(absColumn);

%% ======================== 检测标定板角点 + 计算外参 ========================
fprintf('检测标定板角点...\n');
boardImg = imread(boardImageFile);
if size(boardImg, 3) == 3
    boardGray = rgb2gray(boardImg);
else
    boardGray = boardImg;
end

% 检测棋盘格角点
[imagePtsRaw, boardSizeDetected] = detectCheckerboardPoints(boardGray);
if ~isempty(cameraParams)
    imagePts = undistortPoints(imagePtsRaw, cameraParams);   % 去畸变角点用于解外参
else
    imagePts = imagePtsRaw;
end

if size(imagePts, 1) ~= numRows * numCols
    warning('检测到 %d 个角点，期望 %d 个。请检查棋盘格参数。', ...
        size(imagePts,1), numRows*numCols);
end

% 角点排序：detectCheckerboardPoints 返回的是按行列排序的角点
% 左上角为第1个，按行优先排列

% 构建世界坐标点（Z=0 平面）
[Xw, Yw] = meshgrid(0:numCols-1, 0:numRows-1);
worldPts = [Xw(:) * squareSize, Yw(:) * squareSize];

fprintf('  检测到 %d 个角点\n', size(imagePts, 1));

%% ======================== 计算外参 R, t ========================
% 方法: DLT 计算单应性矩阵 H，然后从 H 和 K 分解出 R, t
% 参考: Zhang (2000), Appendix A

fprintf('计算相机外参...\n');

% 计算单应性矩阵（世界 Z=0 平面 → 图像平面）
H = computeHomography_rm(worldPts, imagePts);

% 从 H 和 K 分解外参
Kinv = inv(K);
h1 = H(:,1); h2 = H(:,2); h3 = H(:,3);

% 缩放因子
lambda = 1 / norm(Kinv * h1);

r1 = lambda * Kinv * h1;
r2 = lambda * Kinv * h2;
r3 = cross(r1, r2);
t  = lambda * Kinv * h3;

% SVD 正交化 R（确保 R 是有效的旋转矩阵）
Rmat = [r1, r2, r3];
[U, ~, V] = svd(Rmat);
Rmat = U * V';

% 确保 det(R) = 1（不是 -1 的反射）
if det(Rmat) < 0
    Rmat(:,3) = -Rmat(:,3);
end

fprintf('  旋转矩阵 R:\n');
disp(Rmat);
fprintf('  平移向量 t:\n');
disp(t);

%% ======================== 构建棋盘格掩码 ========================
% 用最外圈角点构成的多边形定义棋盘格有效区域
% detectCheckerboardPoints 返回角点顺序: 左上→右上→...→左下→右下
% 取四个角: 第1个(左上), 第numCols个(右上), 倒数第numCols个(左下), 最后1个(右下)

% 构建棋盘格角点凸包掩码
cornerPoly = [
    imagePtsRaw(1, :);                          % 左上
    imagePtsRaw(numCols, :);                    % 右上
    imagePtsRaw(end, :);                        % 右下
    imagePtsRaw(end - numCols + 1, :)           % 左下
];

% 略微收缩边界（向内缩5像素，避免包含棋盘格边缘外的背景）
center = mean(cornerPoly, 1);
shrinkRatio = 0.98;  % 收缩比例
cornerPoly = center + (cornerPoly - center) * shrinkRatio;

% 生成棋盘格多边形掩码
boardMask = poly2mask(cornerPoly(:,1), cornerPoly(:,2), imgH, imgW);

fprintf('  棋盘格区域像素数: %d\n', sum(boardMask(:)));

%% ======================== 三层无效像素过滤 ========================
fprintf('过滤无效像素...\n');

% 第①层: 只保留棋盘格内的像素
% 第②层: 只保留格雷码解码有效的像素（validMask）
% 第③层: 只保留相位可靠性足够高的像素
reliablePhase = (reliability > minReliability);

% 综合有效掩码
finalMask = boardMask & validMask & reliablePhase;

nTotal   = numel(finalMask);
nInBoard = sum(boardMask(:));
nValid   = sum(finalMask(:));

fprintf('  总像素: %d\n', nTotal);
fprintf('  棋盘格内: %d (%.1f%%)\n', nInBoard, 100*nInBoard/nTotal);
fprintf('  最终有效: %d (%.1f%% of total)\n', nValid, 100*nValid/nTotal);

%% ======================== 射线-平面求交 ========================
% 对棋盘格内每个有效像素，计算其在标定板上的 3D 坐标
%
% 标定板平面（世界坐标 Z=0）在相机坐标系中的表示:
%   平面法向量: n_c = R * [0;0;1] = R(:,3)
%   平面上一点: t  (平移向量即平面上世界原点在相机坐标系中的位置)
%   平面方程:  n_c' * (P - t) = 0
%
% 从像素 (u,v) 出发的射线（相机坐标系）:
%   d = K^(-1) * [u; v; 1]
%   P = λ * d
%
% 交点: n_c' * (λ*d - t) = 0  →  λ = (n_c' * t) / (n_c' * d)

fprintf('射线-平面求交...\n');

% 平面参数
n_c = Rmat(:, 3);        % 世界 Z 轴在相机坐标系中的方向
t_c = t;                 % 世界原点移到相机坐标系中的位置

% 获取所有有效像素的坐标
[validRows, validCols] = find(finalMask);
nPts = length(validRows);

fprintf('  待求交点: %d 个像素\n', nPts);

% 预分配
worldXYZ = zeros(nPts, 3);
projCols  = zeros(nPts, 1);
camUV     = [validCols, validRows];   % 原始像素坐标
if ~isempty(cameraParams)
    camUVund = undistortPoints(camUV, cameraParams);   % 去畸变后的理想像素坐标
else
    camUVund = camUV;
end

% 批量计算
ndott = n_c' * t_c;

for i = 1:nPts
    u = camUVund(i, 1);   % 去畸变坐标（x 方向）
    v = camUVund(i, 2);   % 去畸变坐标（y 方向）
    
    % 射线方向（相机坐标系）
    d = Kinv * [u; v; 1];
    
    % 交点参数 λ
    ndotd = n_c' * d;
    if abs(ndotd) < 1e-10
        continue;  % 射线平行于平面，跳过
    end
    lambda = ndott / ndotd;
    
    if lambda <= 0
        continue;  % 平面在相机后方，跳过
    end
    
    % 相机坐标系中的交点
    P_c = lambda * d;
    
    % 转换到世界坐标系
    P_w = Rmat' * (P_c - t_c);
    
    worldXYZ(i, :) = P_w';
    projCols(i)    = absColumn(validRows(i), validCols(i));   % 相位图按原始像素索引
end

% 移除未计算成功的点（projCols 仍为 0）
validIdx = projCols > 0;
worldXYZ  = worldXYZ(validIdx, :);
projCols  = projCols(validIdx);
camUV     = camUV(validIdx, :);       % 映射表保留原始像素坐标（project.m 用作索引）
camUVund  = camUVund(validIdx, :);   % 去畸变坐标另存，供 Xc/Yc 使用
nPtsValid = sum(validIdx);

fprintf('  有效交点: %d 个\n', nPtsValid);
fprintf('  世界坐标 X 范围: [%.1f, %.1f] mm\n', min(worldXYZ(:,1)), max(worldXYZ(:,1)));
fprintf('  世界坐标 Y 范围: [%.1f, %.1f] mm\n', min(worldXYZ(:,2)), max(worldXYZ(:,2)));
fprintf('  投影仪列号范围: [%.2f, %.2f]\n', min(projCols), max(projCols));

%% ======================== 重投影验证 ========================
% 将计算出的 3D 点反投影回图像，检查误差

fprintf('重投影验证...\n');
worldHomog = [worldXYZ, ones(nPtsValid, 1)];
projHomog  = (K * (Rmat * worldHomog(:,1:3)' + t_c))';
projPts    = projHomog(:, 1:2) ./ projHomog(:, 3);

reprojErrors = sqrt(sum((projPts - camUVund).^2, 2));   % 与去畸变坐标对比
meanReprojErr = mean(reprojErrors);
maxReprojErr  = max(reprojErrors);

fprintf('  重投影误差 均值: %.4f px, 最大值: %.4f px\n', meanReprojErr, maxReprojErr);

%% ======================== 结果可视化 ========================
figure('Name', '映射关系可视化', 'Position', [50, 50, 1500, 700]);

% (1) 标定板原图 + 检测角点
subplot(2, 3, 1);
imshow(boardGray); hold on;
plot(imagePtsRaw(:,1), imagePtsRaw(:,2), 'r+', 'MarkerSize', 8);
plot(cornerPoly([1:end,1], 1), cornerPoly([1:end,1], 2), 'g-', 'LineWidth', 2);
title('标定板 + 棋盘格角点');

% (2) 投影图案图（如果有）
subplot(2, 3, 2);
if exist(projectedImageFile, 'file')
    projImg = imread(projectedImageFile);
    imshow(projImg); hold on;
    plot(cornerPoly([1:end,1], 1), cornerPoly([1:end,1], 2), 'g-', 'LineWidth', 2);
end
title('投影图案 + 棋盘格边界');

% (3) 综合掩码
subplot(2, 3, 3);
imagesc(finalMask);
axis image; colormap(gca, [0 0 0; 0 1 0]);
title(sprintf('最终有效像素 (%d)', nValid));

% (4) 3D 重建点
subplot(2, 3, 4);
% 降采样显示（避免点太密）
ds = max(1, round(nPtsValid / 5000));
scatter3(worldXYZ(1:ds:end,1), worldXYZ(1:ds:end,2), worldXYZ(1:ds:end,3), 3, projCols(1:ds:end), '.');
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('重建 3D 点 (颜色=投影仪列号)');
axis equal; colorbar; view(3);

% (5) 投影仪列号 vs 世界 X 坐标
subplot(2, 3, 5);
scatter(worldXYZ(1:ds:end,1), projCols(1:ds:end), 3, '.');
xlabel('世界 X (mm)'); ylabel('投影仪列号');
title('映射: X → 投影仪列号');
grid on;

% (6) 重投影误差分布
subplot(2, 3, 6);
histogram(reprojErrors, 50);
xlabel('重投影误差 (px)'); ylabel('频数');
title(sprintf('重投影误差 (均值=%.4f px)', meanReprojErr));

sgtitle('相机像素 → 3D坐标 → 投影仪列号 映射关系');

%% ======================== 部分映射点展示 ========================
fprintf('\n========== 映射表示例（前20个点） ==========\n');
fprintf('  u(px)    v(px)    X(mm)    Y(mm)    Z(mm)    投影仪列号\n');
fprintf('  ------   ------   ------   ------   ------   ----------\n');
nShow = min(20, nPtsValid);
for i = 1:nShow
    fprintf('  %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %8.4f\n', ...
        camUV(i,1), camUV(i,2), worldXYZ(i,1), worldXYZ(i,2), worldXYZ(i,3), projCols(i));
end

%% ======================== 保存结果 ========================
mappingTable = table(camUV(:,1), camUV(:,2), ...
    worldXYZ(:,1), worldXYZ(:,2), worldXYZ(:,3), projCols, ...
    'VariableNames', {'u_px', 'v_px', 'X_mm', 'Y_mm', 'Z_mm', 'projector_col'});

save('pixel_to_projector_mapping.mat', 'mappingTable', 'camUV', 'camUVund', 'worldXYZ', ...
    'projCols', 'finalMask', 'reprojErrors', 'Rmat', 't_c', 'K');
fprintf('\n映射表已保存至 pixel_to_projector_mapping.mat\n');
fprintf('  映射表变量 mappingTable: %d 行 × 6 列\n', nPtsValid);
fprintf('  列: [u_px, v_px, X_mm, Y_mm, Z_mm, projector_col]\n');


% =========================================================================
%                         辅助函数
% =========================================================================

%% 函数：计算单应性矩阵 H (DLT 方法)
function H = computeHomography_rm(worldPts, imagePts)
    N = size(worldPts, 1);
    A = zeros(2*N, 9);
    
    for i = 1:N
        X = worldPts(i, 1);
        Y = worldPts(i, 2);
        u = imagePts(i, 1);
        v = imagePts(i, 2);
        
        A(2*i-1, :) = [X, Y, 1, 0, 0, 0, -u*X, -u*Y, -u];
        A(2*i,   :) = [0, 0, 0, X, Y, 1, -v*X, -v*Y, -v];
    end
    
    [~, ~, V] = svd(A);
    h = V(:, 9);
    H = reshape(h, [3, 3])';
    H = H / H(3,3);
end

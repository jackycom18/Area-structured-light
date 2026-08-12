% =========================================================================
% 张正友标定法 — 仅求解内参矩阵 K
% 参考文献: Z. Zhang, "A Flexible New Technique for Camera Calibration", PAMI 2000
% =========================================================================
clear; close all; clc;

%% ======================== 参数设置 ========================
% 棋盘格参数
squareSize = 25;        % 棋盘格方格边长（mm，按实际修改）
numRows    = 8;         % 棋盘格内角点行数
numCols    = 11;        % 棋盘格内角点列数

% 图像路径（支持 jpg/png/bmp）
imageDir  = 'calib_images/';   % 标定图片所在文件夹
imageExt  = '*.jpg';            % 图片格式

%% ======================== 准备世界坐标 ========================
% 棋盘格在世界坐标系 Z=0 平面上
[Xw, Yw] = meshgrid(0:numCols-1, 0:numRows-1);
worldPts = [Xw(:) * squareSize, Yw(:) * squareSize, zeros(numRows*numCols, 1)];

%% ======================== 读取图像并检测角点 ========================
imageFiles = dir(fullfile(imageDir, imageExt));
if isempty(imageFiles)
    error('在 %s 中未找到 %s 文件', imageDir, imageExt);
end

fprintf('共找到 %d 张图片\n', length(imageFiles));

validIdx = [];
allImagePts = {};   % 每张图检测到的角点

for i = 1:length(imageFiles)
    imgPath = fullfile(imageDir, imageFiles(i).name);
    img = imread(imgPath);
    
    % 检测棋盘格角点
    [imagePts, boardSize] = detectCheckerboardPoints(img);
    
    % 验证是否检测到正确数量的角点
    if size(imagePts,1) == numRows * numCols
        validIdx = [validIdx, i];
        allImagePts{end+1} = imagePts;
        fprintf('  [%2d] %s  ✓ 检测到 %d 个角点\n', i, imageFiles(i).name, size(imagePts,1));
    else
        fprintf('  [%2d] %s  ✗ 角点数不足 (检测到 %d, 期望 %d)\n', ...
            i, imageFiles(i).name, size(imagePts,1), numRows*numCols);
    end
end

numViews = length(allImagePts);
fprintf('\n有效视角数: %d\n', numViews);

if numViews < 3
    error('至少需要 3 个有效视角才能求解内参');
end
%% ======================== 标定内参与畸变系数（张正友法） ========================
% 使用 MATLAB 自带 estimateCameraParameters：
%   内部完成 DLT 单应 + 闭式解 + LM 非线性优化，直接输出 K 与畸变系数
fprintf('\n标定内参 K 与畸变系数...\n');

imagePoints = cat(3, allImagePts{:});   % Mx2xN，所有视角的角点
worldPts2D  = worldPts(:, 1:2);         % 棋盘格平面 Z=0，只取 X,Y

[cameraParams, ~, estimationErrors] = estimateCameraParameters(imagePoints, worldPts2D);

K = cameraParams.IntrinsicMatrix';      % MATLAB 内参矩阵按转置存储
k = cameraParams.RadialDistortion;      % 径向畸变 [k1 k2 k3]
p = cameraParams.TangentialDistortion;  % 切向畸变 [p1 p2]
meanError = estimationErrors.MeanReprojectionError;

fprintf('\n========== 内参矩阵 K ==========\n');
fprintf('  fx = %.4f  (焦距 x 方向，像素)\n', K(1,1));
fprintf('  fy = %.4f  (焦距 y 方向，像素)\n', K(2,2));
fprintf('  cx = %.4f  (主点 x 坐标，像素)\n', K(1,3));
fprintf('  cy = %.4f  (主点 y 坐标，像素)\n', K(2,3));
fprintf('  s  = %.4f  (倾斜因子)\n', K(1,2));
fprintf('\n内参矩阵:\n');
disp(K);
fprintf('\n========== 畸变系数 ==========\n');
fprintf('  径向畸变 k = [%.6f, %.6f, %.6f]\n', k(1), k(2), k(3));
fprintf('  切向畸变 p = [%.6f, %.6f]\n', p(1), p(2));

%% ======================== 重投影误差 ========================
fprintf('\n========== 重投影误差 ==========\n');
fprintf('  平均误差: %.4f 像素\n', meanError);

%% ======================== 保存内参与畸变 ========================
save('camera_intrinsics.mat', 'K', 'k', 'p', 'cameraParams');
fprintf('\n内参 K、畸变系数 k/p、cameraParams 已保存至 camera_intrinsics.mat\n');
fprintf('后续 relationship_mapping.m 将用 cameraParams 对像素坐标去畸变\n');

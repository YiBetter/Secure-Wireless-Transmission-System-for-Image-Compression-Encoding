% 读取原图和压缩后的图像文件
original_image = imread('World.png');  % 原图
compressed_image = imread('compressed_image.jpg');  % 压缩后的图像

% 获取文件的大小（单位：字节）
original_info = dir('World.png');  % 获取原图文件的信息
compressed_info = dir('compressed_image.jpg');  % 获取压缩后图像文件的信息

original_size = original_info.bytes;  % 原图大小
compressed_size = compressed_info.bytes;  % 压缩后图像大小

% 计算压缩率
compression_ratio = original_size / compressed_size;

% 显示压缩率
disp(['压缩比为: ', num2str(compression_ratio)]);

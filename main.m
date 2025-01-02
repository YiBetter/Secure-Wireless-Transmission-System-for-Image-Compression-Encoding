%% 主程序
% jpeg压缩彩色数字图像传输系统
clc;
clear;
close all;

% 让用户选择图像文件
[filename, pathname] = uigetfile({'*.bmp;*.jpg;*.png', 'Image Files (*.bmp, *.jpg, *.png)'}, '选择图像');
if filename == 0
    disp('未选择图像文件');
    return;
end

% 读取选定的图像
img = imread(fullfile(pathname, filename));  % 读取用户选择的图像

tic;

% 信源编码（JPEG压缩）
imwrite(img, 'compressed_image.jpg', 'Quality', 75);  % 进行JPEG压缩，质量为75
compressed_img = imread('compressed_image.jpg');  % 读取压缩后的图像

% 将压缩图像分成RGB三个通道
R_comp = compressed_img(:,:,1);
G_comp = compressed_img(:,:,2);
B_comp = compressed_img(:,:,3);
% 拼接成一个字节串
compressed_bytes = [R_comp(:)',G_comp(:)',B_comp(:)'];

% 生成当前时间戳并使用其生成256位随机密钥
currentTimestamp = posixtime(datetime('now'));  % 获取当前时间戳
rng(currentTimestamp, 'twister');  % 使用时间戳设置随机数种子
key = uint8(rand(1,32)*255);  % 生成256位随机密钥（32字节）
key1 = uint8(rand(1,32)*255);  % （作为错误密钥解密）

% 对数据进行加密
encrypted_comp_bytes = RC4(key, compressed_bytes);

encrypted_data_bin = uint8_to_logical_binSeq(encrypted_comp_bytes);

% 定义数据块的大小
blockSize = 4; % 每个数据块4个数据位

% 计算数据块数量：每个编码块由7位（4数据位 + 3校验位）组成
numBlocks = length(encrypted_data_bin) / blockSize; % 总块数
decodedData = zeros(1,numBlocks * blockSize); % 每个数据块解码后有4个数据位

%总误码比特数
bit_error = 0;

% 信噪比
SNR = 10;

% 流式传输处理
for i = 1:numBlocks
    % 获取当前数据块（4位数据位）
    dataBlock = encrypted_data_bin((i-1)*blockSize + 1:i*blockSize);

    % 汉明编码（每4位数据编码为7位）
    encodedData = hamming_encode(dataBlock);

    % BPSK调制
    bpskMod = @(x) 2*x - 1;  % BPSK: 0 -> -1, 1 -> 1
    txSig = bpskMod(encodedData);  % 使用BPSK调制

    % 信道传输（添加高斯噪声）
    rxSig = awgn(txSig, SNR, 'measured');  % 添加AWGN噪声

    % BPSK解调
    bpskDemod = @(x) (x >= 0);  % 解调：信号 >= 0 -> 1，信号 < 0 -> 0
    rxData = bpskDemod(rxSig);  % 使用BPSK解调

    % 汉明解码
    [decodedBlock, error] = hamming_decode(rxData);
    if(error)
        %disp(['第',num2str(i),'块数据出错，已纠正。']);
        bit_error = bit_error + 1;
    end
    
    % 将解码后的数据块存储到预分配的数组中
    decodedData((i-1)*blockSize + 1:i*blockSize) = decodedBlock;
end

%比特误码率
ber = bit_error / (numBlocks * 7);
disp(['比特误码率为：',num2str(ber)]);

% 解密操作：对接收到的二进制数据进行解密并转换为像素值
decrypted_data_bin = RC4(key1, logical_binSeq_to_uint8(logical(decodedData),8));

m = size(R_comp,1);
n = size(R_comp,2);
% 由于图像RGB通道矩阵维度一致，所以可以全部按照R通道处理
R_rx = decrypted_data_bin(1:m * n);  % 取 R 通道的像素值
G_rx = decrypted_data_bin(m * n + 1 : m * n + m * n);  % 取 G 通道的像素值
B_rx = decrypted_data_bin(m * n + m * n + 1 : end);  % 取 B 通道的像素值
R_rx = reshape(R_rx,m,n);
G_rx = reshape(G_rx,m,n);
B_rx = reshape(B_rx,m,n);

% 合并RGB通道
img_rx = cat(3, R_rx, G_rx, B_rx);
imwrite(img_rx,'received_image.jpg');

toc;% 数据传输总共耗时

% 显示原始图像和接收图像
figure;
subplot(1, 2, 1);
imshow(img);
title('原始图像');
subplot(1, 2, 2);
imshow(img_rx);
title('接收图像');

%% 辅助函数
function binSeq = uint8_to_logical_binSeq(uint8Vec)
    % 将 uint8 类型的行向量转换为逻辑值（二进制序列）
    % uint8Vec: uint8 类型的行向量
    % binSeq: 转换后的逻辑值二进制序列
    
    % 获取行向量的长度
    numElements = length(uint8Vec);
    
    % 初始化逻辑值二进制序列
    binSeq = false(1, numElements * 8);  % 预分配空间，长度是 uint8Vec 的 8 倍
    
    % 遍历每个 uint8 元素
    idx = 1;  % 二进制序列的索引
    for i = 1:numElements
        % 将当前 uint8 元素转换为逻辑值（0 或 1）
        binStr = dec2bin(uint8Vec(i), 8) == '1';  % 转换为逻辑值数组
        
        % 将转换后的逻辑值二进制序列存入 binSeq
        binSeq(idx:idx+7) = binStr;
        
        % 更新索引
        idx = idx + 8;
    end
end

function bytes_sequence = logical_binSeq_to_uint8(binStr,bitsPerElement)
    % 将逻辑值（二进制序列）转换为 uint8 类型的行向量
    numElements = length(binStr) / bitsPerElement;
    bytes_sequence = zeros(1,numElements, 'uint8');
    for i = 0:numElements-1
        startIdx = i * bitsPerElement + 1;
        binNum = binStr(startIdx:startIdx+bitsPerElement-1);
        bytes_sequence(i+1) = bit2dec(binNum,0);  % 转换回十进制
    end
end

function [ y ] = bit2dec(x ,dir)
    % 8位二进制数转换成 十进制数
    % y为dec十进制输出,x为输入的8位二进制数组
    % d7 d6 d5 d4 d3 d2 d1 d0
    % dir = 1(高位是d0,低位是d7)
    %     = 0(高位是d7,低位是d0)

    y = 0;

    for i = 1:8
        if( dir == 1)
            y = y + x(i)*2^(i-1) ;
        else
            y = y + x(i)*2^(8-i);
        end
    end
end

% 汉明编码函数
function encodedData = hamming_encode(data)
    % 假设输入是一个4位二进制数据
    % 数据按位置插入检验位
    data = data(:); % 确保数据是列向量
    encodedData = zeros(7, 1);  % 编码后的7位数据
    
    % 数据位置: 3、5、6、7
    encodedData([3, 5, 6, 7]) = data;
    
    % 计算检验位
    encodedData(1) = mod(encodedData(3) + encodedData(5) + encodedData(7), 2);  % P1
    encodedData(2) = mod(encodedData(3) + encodedData(6) + encodedData(7), 2);  % P2
    encodedData(4) = mod(encodedData(5) + encodedData(6) + encodedData(7), 2);  % P4
    
    encodedData = encodedData';
end

% 汉明解码的错误检测和纠正
function [decodedBlock, errorDetected] = hamming_decode(block)
    % 计算检验位
    P1 = mod(block(1) + block(3) + block(5) + block(7), 2);  % 计算第1检验位
    P2 = mod(block(2) + block(3) + block(6) + block(7), 2);  % 计算第2检验位
    P4 = mod(block(4) + block(5) + block(6) + block(7), 2);  % 计算第4检验位
    
    % 错误位置（1基索引），如果检验位的和为1，则对应的位发生错误
    errorPos = P1 * 1 + P2 * 2 + P4 * 4;
    errorDetected = (errorPos ~= 0);  % 如果errorPos不为0，则检测到错误
    
    if errorDetected
        % 校正错误
        block(errorPos) = mod(block(errorPos) + 1, 2);
    end
    
    % 提取数据位（即第3、5、6、7位是数据位）
    decodedBlock = block([3, 5, 6, 7]);
end

function cipherText = RC4(key, data)
    % RC4 加密/解密函数
    % key: 密钥，字节序列
    % data: 要加密或解密的字节数据
    % return: 加密/解密后的字节数据
    
    % 密钥调度算法（KSA）
    S = 0:255;  % 初始化 S 数组
    j = 0;      % 初始化 j
    key_len = length(key);  % 获取密钥长度
    
    for i = 0:255
        j = mod(j + S(i + 1) + key(mod(i, key_len) + 1), 256);
        % 交换 S[i] 和 S[j]
        temp = S(i + 1);
        S(i + 1) = S(j + 1);
        S(j + 1) = temp;
    end
    
    % 伪随机生成算法（PRGA）
    i = 0;
    j = 0;
    cipherText = uint8(zeros(1, length(data)));  % 初始化密文数组
    
    for k = 1:length(data)
        i = mod(i + 1, 256);
        j = mod(j + S(i + 1), 256);
        
        % 交换 S[i] 和 S[j]
        temp = S(i + 1);
        S(i + 1) = S(j + 1);
        S(j + 1) = temp;
        
        % 生成伪随机字节流
        keystream = S(mod(S(i + 1) + S(j + 1), 256) + 1);
        
        % 将明文与密钥流异或，得到密文（或解密）
        cipherText(k) = bitxor(data(k), keystream);
    end
end

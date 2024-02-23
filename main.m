%{
    Copyright (c) 2023 BINGBING DAN

    Author: Bingbing Dan
    Email: danbingbing20@mails.ucas.ac.cn
    Affiliation: University of Chinese Academy of Sciences

    Corresponding Publication:
    Bingbing Dan, et al. "Dynamic Weight Guided Smooth-Sparse Decomposition 
    for Small Target Detection against Strong Vignetting Background"
    IEEE Transactions on Instrumentation and Measurement, 2023

    Description:
    This code is a supplementary material for the above-mentioned publication. It implements the DW-SSD model
    described therein and provides a practical example of the concepts presented.

%}

clc;
clear;
close all;

len_of_seq = 6;
N_k = [];
for ii = 1:len_of_seq

    img = imread(['images\' num2str(ii) '.bmp']);
    
    alg = SSD;
    alg.O = mat2gray(img);
    alg.k = 5;
    alg.M = N_k;
    alg = alg.process();
    N_k = cat(3,N_k,alg.O-alg.B);

    if ii > alg.k
        figure;imshow(alg.O,[]);title('Original Image')
        figure;imshow(alg.B,[]);title('Vignetting Backround')
        figure;imshow(alg.T,[]);title('Small Target')
    end
end
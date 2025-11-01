%A multi-level image thresholding algorithm based on Kapur's Entropy using a
%metaheuristic algorithm (SSNMRA)
% --------------------------------------------------------------- %
function [xBestR,BestvalR,Convergence_curve,time]=SSNMRA_kapur
tic;
% Number of search agents
n=30;
N_Iter=1000;
% Number of different solutions
nd=3;
I=imread('225017.jpg');
figure(1)
%hold on
subplot(1,3,1)
imshow(I);
title('Input test image')

a='225017'
I=rgb2gray(I);
[n_countR, x_valueR] = imhist(I(:,:,1));
Nt = size(I,1) * size(I,2); 
Lmax = 255;   %256 different maximum levels are considered in an image (i.e., 0 to 255)

% Distribution histogram 0 - 256 
for i = 1:Lmax
        probR(i) = n_countR(i) / Nt;
end
 N_IterTotalR=N_Iter;
%Lower bounds and Upper bounds
LbR=zeros(1,nd); 
UbR=Lmax*ones(1,nd);  %(here it is from 0 to 255)
fitnessR=zeros(n,1);  %zeros
[xBestR,BestvalR,Convergence_curve]=SSNMRA(n,N_IterTotalR,LbR,UbR,fitnessR,probR,nd);
%% Same as above for RGB images treating each components seperately and initializing
%% Displaying segmented output
    gBestR = sort(xBestR);
    Iout = imageGRAY(I,gBestR);
rgbImage = ind2rgb(Iout, jet(256));


%Show results on images
figure(1)
%hold on
subplot(1,3,2)
    imshow(Iout);
    
    title('Gray scale theshold image')
figure(1)

subplot(1,3,3)
imshow(rgbImage);
title('Theshold image with psedo color')

    Iout2 = mat2gray(Iout);
    %show results
    intensity = gBestR;     %intensity  value of the best elemento (TH)
    PSNRV = PSNR(I, Iout);  %PSNR between original image I and the segmented image Iout
    
    
    
    %plot the threslhol over the histogram
    figure(2)
    plot(probR)
    xlabel('Grey level')
    ylabel('Frequency')
    legend('Histogram')
    
    %hold on
    vmax = max(probR);
    for i = 1:nd
        color=rand(1,3);
        line([intensity(i), intensity(i)],[0 vmax],[1 1],'Color',color,'Marker','.','LineStyle','-','LineWidth',2);  
    end
    legend('Histogram','Threshold 1','Threshold 2','Threshold 3','Threshold 4','Threshold 5','Threshold 6','Threshold 7')    
    
    figure(3)
    hold on      
    plot(Convergence_curve,'DisplayName','PO','Color', 'r','LineWidth',2);
    xlabel('Iteration')
    ylabel('Objective function value')
     

time=toc 
%close all

 
 

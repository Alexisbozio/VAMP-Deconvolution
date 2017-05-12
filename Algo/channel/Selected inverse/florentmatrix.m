clear

seed=int64(14753);
rng(seed);                  %control the randomness 
rho=0.01;
delta=1e-6;

%Generate image
n=64;m=32;
k=int64(ceil(rho * n*n));     %k=sparsity
x = zeros(n*n,1);             %size of the image
beq = randperm(n*n);          %distributes randomly positions of non zero component in x
supp = beq(1:k);            %build a vector with k first elements of beq (bequille)
x(supp) = 1.;               %put 1 in every position given by supp 

PSF = fspecial('gaussian', [n n], 2);

PSF = PSF/max(PSF(:));

 

% Define Matrix BLUR

BLUR=zeros(n*n,n*n);

 

for i=1:n

    for j=1:n

        %here I need a PSF centered in (i,j)

        IM=circshift(PSF,[i-1-n/2,j-1-n/2]);

        num=i+(j-1)*n;        

        BLUR(num,:)=reshape(IM,1,n*n);

    end

end

    



 

IMAGE_BLURED=reshape(BLUR*x,n,n);

 

% Cool, now we need to downsize

v = ones(1,alpha);

M = [];

for I = 1:m

    M = blkdiag(M,v);

end

 

FINAL_IMAGE=M*IMAGE_BLURED*M';

figure(1); 
    subplot(221); 
    imx = imshow(IMAGE);
    
    title('original');
    
    subplot(222); 
    imz = imshow(IMAGE_BLURED);
    
    title('convoluted');
    
    subplot(223); 
    imy = imshow(FINAL_IMAGE);
    
    title('transformed');
%{
Of course, I could also do as such:

 

FINAL_OPERATOR=zeros(m*m,n*n);

for i=1:n*n

    FINAL_OPERATOR(:,i)=reshape((M*reshape(BLUR(:,i),n,n)*M'),m*m,1);

end

FINAL_IMAGE2=reshape(FINAL_OPERATOR*reshape(IMAGE,n*n,1),m,m);

 

Y=reshape(FINAL_IMAGE,m*m,1);

Y=Y-mean2(Y);

FINAL_OPERATOR=FINAL_OPERATOR-mean2(FINAL_OPERATOR);

 

tic;[U,S,V]=svd(FINAL_OPERATOR);toc;

%}

 


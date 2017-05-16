%Super_resolution reconstruction with Vamp on pictures by treating the
%image by block

clear
seed=int64(14753);

addpath(genpath('..'));


profile on

denoiser='BM3D';%Available options are NLM, Gauss, Bilateral, BLS-GSM, BM3D, fast-BM3D, and BM3D-SAPCA 
filename='lena.png';
ImIn=double(imread(filename));

%Parameters of the problem
delta=1e-6;                  %variance of noise
height=200;                  %length side of original image
width=height;
downsizing=2;               % downsizing factor
blur_var=1;                  % variance of Gaussian kernel
blur_width=blur_var*6-1;     % width of Gaussian kernel (3.03 deviation from the mean enclose 99% of the pdf)

down_height=height/downsizing;
n=height^2;
m=down_height^2;

%Create the image
init_1=200;
init_2=200;
x_0=ImIn(init_1:(init_1-1+height),init_2:(init_2-1+width));
Imsize=size(x_0);
x0=x_0(:);
dr = (max(x0)-min(x0));
x0=(x0-min(x0))./dr;

%Rescale the image in range [0,1]
x_0=x_0./dr;      

%Plot
figure(1);subplot(221);imagesc(x_0);title('Original Image');
    axis image;colormap gray;colorbar();pause(0.1);
    
%Generate Super-resolution Matrix F
F= proj_sr( n, downsizing, blur_width, blur_var );

% Transformed image y
y = F*x0 + sqrt(delta) .* randn(m,1) ;

y_square=reshape(y, [down_height, down_height]);

%Plot
figure(1);subplot(222);
imy_sq = imagesc(y_square);
colormap gray;axis image;colorbar();pause(0.1);title('Transformed');


opts.channel_prmts=delta;

opts.t_max=60;
opts.eps_conv=1e-13;
opts.damp_meas = 0.5;

r=@channel_awgn;            %Using Selected Inversion
s=@channel_wood;            %Using Woodsbury Identity
opts.channel=r;

x_hat1=zeros(n,1);
A=zeros(n,1);
B=zeros(n,1);


%Block configuration
n_blocks=5;
size_block=down_height/n_blocks;


j_axis=randperm(n_blocks);
i_axis=randperm(n_blocks);

for j=j_axis
    for i=i_axis
            
             %Select the block
            rows=((i-1)*size_block+1:i*size_block);
            column=((j-1)*size_block+1:j*size_block);
            block=zeros(down_height);
            block(rows,column)=1;
            sel_el_m=logical(reshape(block,[],1));
            y_vec=y(sel_el_m,:);
            
            %Build the matrix that is talking to the block
            F_ij=F(sel_el_m,:);
            sel_el_n=logical(sum(logical(F_ij)));
            m_red=nnz(sel_el_m);
            n_red=nnz(sel_el_n);
                    
            F_red = F_ij(:, sel_el_n);
            Ap=A(sel_el_n,:);
            Bp=B(sel_el_n,:);
            
            x_red=x0;
            x_red(~sel_el_n)=0;
             
            x_trek=logical(reshape(x_red, [height, height]));
            height_red=nnz(sum(x_trek));
            weight_red=nnz(sum(x_trek,2));
      
            opts.height=height_red;
            opts.width=weight_red;
            opts.x0=x0(sel_el_n);
            fprintf('o Running V-AMP...\n')
            tic;
            [x_hat_red,mse, A1,B1] = vamp_image(y_vec, F_red, Ap, Bp, opts);
            toc 
            
            %Update priors
            A(sel_el_n)=Ap+A1;
            B(sel_el_n)=Bp+B1;
            x_hat1(sel_el_n)=x_hat_red;
            
                              
            %Print
            figure(3);subplot(221);ima = imagesc(reshape(x_hat1, [height, height]));
            colormap gray;axis image;pause(0.1);colorbar();title('estimate with vamp');
            
           
    end
end    



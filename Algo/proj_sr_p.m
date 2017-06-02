function F = proj_sr_p( n, block_width,pa ,fwhm)
%PROJ_SR Build (sparse) super-resolution matrix.
%
%   n:              number of pixels (i.e., n = l * l)
%   block_width:    size of block to be decimated (e.g. 2 for 2x2 block)
%   blur_width:     width of Gaussian kernel (e.g. 9 for 9x9 kernel)
%   blur_var:       variance of Gaussian kernel
%
len = sqrt(n);
if len - floor(len) ~= 0
    error('sqrt(n) must be an integer.')
end

%{ 
Lets see if we already have this projector stored away somewhere...
proj_file_name = sprintf('imageproc/saved_proj/H_%d_%d_%d_%0.5e.mat',sqrt(n),block_width,blur_width,blur_var);
if exist(proj_file_name,'file')
	load(proj_file_name);
	return;
else
    %}




    pixelPitch = sqrt(pa);
    srPixelPitch = pixelPitch / block_width;   
    fwhmInPixels = fwhm / srPixelPitch ;
    psfSigmaInPixels = 2 * sqrt(2*log(2)) * fwhmInPixels;
    cutoff=1e-4;
    H = fspecial('gaussian', 129, psfSigmaInPixels);
    edgeIdx=find(sum(H>cutoff,1));
    edgeIdx=edgeIdx(1);
    c = ceil(129 / 2) ;     
    radius  = c - edgeIdx;
    
    if rem(radius,2) == 1
        radius=radius+1;
    end

    H = H(c-radius:c+radius, c-radius:c+radius);
    H(H < cutoff) = 0.0;     

  
	% Generate convolution matrix for averaging/blurring kernel
	
	sum_kern = ones(block_width);
	kern = conv2(H,sum_kern);
	% kern = conv2(gauss_blur, fspecial('average', block_width));

	F = conv_mat(kern, len);

	% Remove appropriate rows from convolution matrix, emulating decimation
	rows = [];
	for i = 0:ceil(len / block_width - 1)
	   rows = [rows, i * len * block_width + (1:block_width:len) + ...
	       + (1 - mod(block_width, 2)) * (.5 * block_width - 1) * (len + 1)];
	end
	F = F(rows, :);


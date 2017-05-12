function C = generate_convolution_matrix( n, blur_width, blur_var )
%PROJ_SR Build (sparse) super-resolution matrix.
%
%   n:              number of pixels (i.e., n = l * l)
%   block_width:    size of block to be decimated (e.g. 2 for 2x2 block)
%   blur_width:     width of Gaussian kernel (e.g. 9 for 9x9 kernel)
%   blur_var:       variance of Gaussian kernel
%
%block_width=1;
len = sqrt(n);
if len - floor(len) ~= 0
    error('sqrt(n) must be an integer.')
end

	kern = fspecial('gaussian', blur_width, blur_var);
	

	C = conv_mat(kern, len);

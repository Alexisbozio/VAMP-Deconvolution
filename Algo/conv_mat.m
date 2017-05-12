
function mat = conv_mat( filter, image_len )
%CONV_MAT Generate 2D convolution matrix, as in mode 'same' of conv2.
%
mat = convmtx2(filter, image_len, image_len);

% Remove rows appropriately so as to have a n x n matrix
pad = ceil(.5 * (length(filter) - 1.));
im = padarray(ones(image_len), [pad pad]);
if ~mod(length(filter), 2)
    im = im(1:(end - 1), 1:(end - 1));
end
rows = logical(im(:));
mat = mat(rows, :);
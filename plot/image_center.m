function output = image_center(input, window_size)

output = zeros([window_size window_size]);
size1 = size(input,1)-1;
size2 = size(input,2)-1;
offset1 = round((window_size-size1)/2);
offset2 = round((window_size-size2)/2);
output(offset1:(offset1+size1), offset2:(offset2+size2)) = input;

end

function [filtered_image] = applyHannFilt (image,factor)

if ndims(image)>10
    error('too many dimensions!')
end

if nargin < 2
    factor = 1;
end

filtered_image = zeros(size(image));
for i1=1:size(image,4)
    for i2=1:size(image,5)
        for i3=1:size(image,6)
            for i4=1:size(image,7)
                for i5=1:size(image,8)
                    for i6=1:size(image,9)
                        for i7=1:size(image,10)
                            for i8=1:size(image,11)
                                for i9=1:size(image,12)
                                    kspace =ifftshift(ifftn(fftshift(image(:,:,:,i1,i2,i3,i4,i5,i6,i7,i8,i9))));
                                    filtered_kspace = hannFilt(kspace, 'sep', factor);
                                    filtered_image(:,:,:,i1,i2,i3,i4,i5,i6,i7,i8,i9) = fftshift(fftn(ifftshift(filtered_kspace)));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Line Of Sight Propagation Map Only
Rx.LosElectricField = 10*log10(abs(Rx.LosElectricField));
Rx.LosElectricField(find(isinf(Rx.LosElectricField) == 1)) = 0;
imageEF = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end
for i = 1:mesh_.zNodeNum
    Rx.TotalEFLayer(:,i) = (Rx.LosElectricField(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):(i*...
        (mesh_.xNodeNum.*mesh_.yNodeNum))));
    imageEF = mat2gray(reshape(Rx.TotalEFLayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageEFScaled = imresize(imrotate(imageEF,90),imageEFScale);
    imageEFScaledOverlayed = imoverlay(imageEFScaled,structImage,[0,0,0]);
    imshow(rgb2gray(imageEFScaledOverlayed));
    colorbarLabels = min((Rx.TotalEFLayer(:,i))) + (0:5) .* ((max(Rx.TotalEFLayer(:,i))-...
        min(Rx.TotalEFLayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title(['EF LOS Only']);
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end
    
%% Reflection Propagation Map Only
Rx.ReflectEF = 10*log10(abs(RxReflectEF(1,:)));
% Rx.ReflectEF = RxReflectEF(1,:);
Rx.ReflectEF(find(isinf(Rx.ReflectEF) == 1)) = 0;
Rx.ReflectEF(find((Rx.ReflectEF) == 0)) = min(min(Rx.ReflectEF));
imageEF = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end
for i = 1:mesh_.zNodeNum
    Rx.TotalEFLayer(:,i) = (Rx.ReflectEF(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):...
        (i*(mesh_.xNodeNum.*mesh_.yNodeNum))));
    imageEF = mat2gray(reshape(Rx.TotalEFLayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageEFScaled = imresize(imrotate(imageEF,90),imageEFScale);
    imageEFScaledOverlayed = imoverlay(imageEFScaled,structImage,[0,0,0]);
    imshow(rgb2gray(imageEFScaledOverlayed));
    colorbarLabels = min((Rx.TotalEFLayer(:,i))) + (0:5) .* ((max(Rx.TotalEFLayer(:,i))-...
        min(Rx.TotalEFLayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title('First Reflections Only');
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end

%% Second Reflection Propagation Map Only
Rx.ReflectEF = 10*log10(abs(RxReflectEF(2,:)));
Rx.ReflectEF(find(isinf(Rx.ReflectEF) == 1)) = 0;
Rx.ReflectEF(find((Rx.ReflectEF) == 0)) = min(min(Rx.ReflectEF));
imageEF = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end
for i = 1:mesh_.zNodeNum
    Rx.TotalEFLayer(:,i) = (Rx.ReflectEF(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):(i*...
        (mesh_.xNodeNum.*mesh_.yNodeNum))));
    imageEF = mat2gray(reshape(Rx.TotalEFLayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageEFScaled = imresize(imrotate(imageEF,90),imageEFScale);
    % To Take care of the Blackage
    blackageMask = (imageEF == 1);
    blackageMaskScaled = imresize(imrotate(blackageMask,90),imageEFScale);
    blackageMaskScaled = blackageMaskScaled | structImage;
    imageEFScaledOverlayed = imoverlay(imageEFScaled,blackageMaskScaled,[0,0,0]);
    imshow(rgb2gray(imageEFScaledOverlayed));
    colorbarLabels = min((Rx.TotalEFLayer(:,i))) + (0:5) .* ((max(Rx.TotalEFLayer(:,i))-...
        min(Rx.TotalEFLayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title('Second Reflections Only');
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end
    
%% Reflection & Line Of Signt Propagation Map
Rx.TotalEF = 10*log10(abs(Rx.LosElectricField + (reflectExaggerationFac * Rx.ReflecEF(1,:)) +...
 (reflectExaggerationFac *Rx.ReflecEF(2,:)) ));
imageEF = zeros(mesh_.yNodeNum,mesh_.xNodeNum,3);
if mesh_.zNodeNum ~= 1
    zplaneHeight =  linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum);
end

for i = 1:mesh_.zNodeNum
    Rx.TotalEFLayer(:,i) = (Rx.TotalEF(((i-1)*(mesh_.xNodeNum.*mesh_.yNodeNum)+1):(i*...
    (mesh_.xNodeNum.*mesh_.yNodeNum))));
    imageEF = mat2gray(reshape(Rx.TotalEFLayer(:,i),mesh_.yNodeNum,mesh_.xNodeNum));
    figure 
    colormap gray
    imageEFScaled = imresize(imrotate(imageEF,90),imageEFScale);
    imageEFScaledOverlayed = imoverlay(imageEFScaled,structImage,[0,0,0]);
    imshow(rgb2gray(imageEFScaledOverlayed));
    colorbarLabels = min((Rx.TotalEFLayer(:,i))) + (0:5) .* ((max(Rx.TotalEFLayer(:,i))-...
        min(Rx.TotalEFLayer(:,i)))./5);
    colorbar('YTickLabel',num2str((colorbarLabels'),'%10.1f'));
    title(['Z-Plane height of ',num2str(zplaneHeight(i),'%10.2f'),'; LOS = ',...
        num2str(losFlag),'; Reflec = ',...
        num2str(reflectionFlag),'; Tx Power = ',num2str(Tx.power')]);
    if grayScaleImage == 0
        colormap(gca,'jet');
    end
end

figure
surf(X,Y,imageEF)
colorbar

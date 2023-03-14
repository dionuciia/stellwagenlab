warning('off','all')
tic
dirs = dir('output*');
for j = 1:length(dirs)
    DIR(j,:) = dirs(j).name;
end
for k = 1:length(dirs)
    cd(DIR(k,:))

    fovie = dir('*filtered.mat');
    binary = dir("Binary_*.mat");
    load(fovie.name);
    load(binary.name);
    mkdir filtered_movie\
    mkdir binary\
    disp([' ']);
    disp(['Extracting pre-processed movie ',sprintf('%d',k),'...']);

    cd filtered_movie\
    try
        for i=1:length(A_dFoF)
        filename = ['fovie_frame',sprintf('%04d',i),'.jpg'];
        imwrite(A_dFoF(:,:,i), filename);
        end
    catch 
        i = i-1;
        disp(['End of ',num2str(i),' filtered frames.']);
    end 
    cd ..\binary\
    try
        for i=1:length(BW_ppA)
        filename = ['binary_frame',sprintf('%04d',i),'.jpg'];
        imwrite(BW_ppA(:,:,i), filename);
        end
    catch 
       i = i-1;
     disp(['End of ',num2str(i),' binary frames.']);
    end 

    cd('../../');
end
toc
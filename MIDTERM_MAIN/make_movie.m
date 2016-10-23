
file_name_core = 'matlab_data_file';
i = 1;
while (exist(strcat(file_name_core,num2str(i),'.txt'),'file')==2)

    filename=strcat(file_name_core,num2str(i),'.txt');
    delimiter = ' ';
    fid = fopen(filename,'r');
    tLines = fgets(fid);
    numCols = length(strsplit(tLines,' '))-2;
    fclose(fid);

    fid = fopen(filename,'r');
    C=textscan(fid,[repmat('%d',[1,numCols])],'CollectOutput',1);

    a = C{1};
    image(a.*255);
    F(i) = getframe;
    colormap(winter);
    set(gcf,'color','w');

    if(exist('PNG','dir')~=7)
        mkdir('PNG')
    end
    names = strsplit(filename,'/');
    figurename = names(length(names));
    figurename = figurename(1:end-4);
    saveas(gcf,strcat('PNG/',filename,'.png'));
    fclose(fid);
    i = i+1;
end

movie(F)
movie2avi(F,strcat(file_name_core,'.avi'));
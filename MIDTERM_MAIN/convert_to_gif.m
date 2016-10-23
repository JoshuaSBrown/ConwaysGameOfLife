function convert_to_gif(filename)
delimiter = ' ';
fid = fopen(filename,'r');
tLines = fgets(fid);
numCols = length(strsplit(tLines,' '))-2;
fclose(fid);

fid = fopen(filename,'r');
C=textscan(fid,[repmat('%d',[1,numCols])],'CollectOutput',1);

a = C{1};
image(a.*255);
colormap(winter);
set(gcf,'color','w');

mkdir('PNG')
names = strsplit(filename,'/');
figurename = names(length(names));
figurename = figurename(1:end-4);
saveas(gcf,strcat('PNG/',filename,'.png'));
fclose(fid);

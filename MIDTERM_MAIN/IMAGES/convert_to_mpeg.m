
Manifold=zeros(900,900,3,51);

name='fileout';

for x=0:50
   file = strcat(name,num2str(x),'.jpeg');
   Manifold(:,:,1,x+1)=imread(file);
end

Manifold = Manifold./255;
% create the video writer with 1 fps
 writerObj = VideoWriter('myVideo.avi');
 writerObj.FrameRate = 1;
 % set the seconds per image
 secsPerImage = (linspace(0,50,51)+1)/2;
 % open the video writer
 open(writerObj);
 % write the frames to the video
 for u=1:50
     % convert the image to a frame
     frame = im2frame(Manifold(:,:,:,u));
     for v=1:secsPerImage(u) 
         writeVideo(writerObj, frame);
     end
 end
 % close the writer object
 close(writerObj);

function spinodal_decomp(D,gamma,options)

arguments
   D double = 10
   gamma double = 5
   options.dt double = 0.005
   options.GridSize double = 200
   options.NumIterations double = 10000
   options.FrameSpacing double = 10
   options.CaptureMode char = 'standard'
   options.ImgStyle char = 'binary'
   options.Colormap char = 'pink'
   options.FileName char = 'spinodal_decomposition'
end

% Random Starting Concentrations (-1 and 1 represent different species)
u = 2*randi(2,options.GridSize)-3;

writer = VideoWriter(options.FileName);
open(writer);
figure(1)
colormap(options.Colormap)
image(u,'CDataMapping','scaled');
frame = getframe(1);
writeVideo(writer,frame);

% Variables for incremental capture
count = 1;
frameStep = 1;

for i = 1:options.NumIterations
   u = iterate(u,D,gamma,options.dt);
   % Incremental video mode
   if (strcmp(options.CaptureMode,'incremental'))
       if (count == frameStep)
           if (strcmp(options.ImgStyle,'true'))
               image(u,'CDataMapping','scaled');
           else
               uMod = round((u+1)/2); % Binarizes image
               image(uMod,'CDataMapping','scaled');
           end
           frame = getframe(1);
           writeVideo(writer,frame);
           count = 0;
           frameStep = frameStep + 1;
       end
       count = count+1;
   % Standard video mode
   else
       if (mod(i,options.FrameSpacing) == 0)
           if (strcmp(options.ImgStyle,'true'))
               image(u,'CDataMapping','scaled');
           else
               uMod = round((u+1)/2);
               image(uMod,'CDataMapping','scaled');
           end
           frame = getframe(1);
           writeVideo(writer,frame);
       end
   end
end

close(writer);

fprintf('Done!\n');

% Forward Euler Method iteration of model
function uOut = iterate(u,D,gamma,dt) 
    % Calculates laplacian of concentration field
    uLaplace = laplacian(u);
    % Calculates chemical potentials
    uMu = u.^3 - u - gamma*uLaplace;
    % Laplacian of chemical potentials
    muLaplace = laplacian(uMu);
    % Cahn-Hilliard Equation
    duT = D*muLaplace;
    % Foreward Euler Method
    uOut = u + dt*duT;
end 

end

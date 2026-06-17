figure(3);
theme(gcf,"light");
set(gca,"NextPlot","replacechildren");

times=arrayfun(@(s) md.hist{s}.t,1:length(md.hist));

T0=0;
T1=18;
N=20;

video=VideoWriter("movie.avi");
video.FrameRate=300/10;
%video.FrameRate=10;
open(video);

mxal=md.hist{end}.al{2,1};
mxal=mxal(:);
mxal=max(mxal);

for t=linspace(T0,T1,N)
    [~,s]=min(abs(times-t));
    h=md.hist{s};

    clf;

    geo.plotfield(h.al{2,1},2,@(al) al,displ=h.displ,scale=0);
    shading interp;
    pbaspect([1,1,1]);
    daspect([1,1,1]);
    axis off;
    %xlim([-1.1*2,1.1*2]);
    %ylim([-1.1*1,1.1*1]);
    colormap parula(30);
    colorbar;
    %clim([0,0.06]);
    title(sprintf("t=%.5f",t));

    drawnow;
    frame=getframe(gcf);
    writeVideo(video,frame);
end

close(video);
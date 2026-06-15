figure(3);
theme(gcf,"light");
set(gca,"NextPlot","replacechildren");

times=arrayfun(@(s) md.hist{s}.t,1:length(md.hist));

T0=0;
T1=33.64;
N=200;

video=VideoWriter("clip.avi");
video.FrameRate=300/10;
%video.FrameRate=10;
open(video);

mxal=md.hist{end}.al{2,1};
mxal=mxal(:);
mxal=max(mxal);

mxte=max(arrayfun(@(s) max(reshape(pagenorm(md.hist{s}.Te,"fro"),[],1)),1:length(md.hist)));

for t=linspace(T0,T1,N)
    [~,s]=min(abs(times-t));
    h=md.hist{s};

    clf;

    subplot(1,2,1);

    geo.plotfield(h.al{2,1},2,@(al) al,displ=h.displ,scale=10);
    shading interp;
    pbaspect([1,1,1]);
    daspect([1,1,1]);
    axis off;
    xlim([-1.1*2,1.1*2]);
    ylim([-1.1*1,1.1*1]);
    colormap parula(30);
    colorbar;
    clim([0,mxal]);
    title(sprintf("t=%.5f",t));

    subplot(1,2,2);

    geo.plotfield(h.Te,2,@(Te) pagenorm(Te,"fro"),displ=h.displ,scale=10);
    shading interp;
    pbaspect([1,1,1]);
    daspect([1,1,1]);
    axis off;
    xlim([-1.1*2,1.1*2]);
    ylim([-1.1*1,1.1*1]);
    colormap parula(30);
    colorbar;
    clim([0,mxte]);
    title(sprintf("t=%.5f",t));

    drawnow;
    frame=getframe(gcf);
    writeVideo(video,frame);
end

close(video);
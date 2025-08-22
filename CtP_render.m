%% Code adapted from Inverted Cart Pole repository on github
function CtP_render(t,X,m,M,L)
	for j=1:length(t)
		
		x = X(j,1);
		th = X(j,3);

		% dimensions
		% L = 2;  % pendulum length
		W = 1*sqrt(M/5);  % cart width
		H = .5*sqrt(M/5); % cart height
		wr = .2; % wheel radius
		mr = .3*sqrt(m); % mass radius

		% positions
		y = wr/2+H/2; % cart vertical position
		w1x = x-.9*W/2;
		w1y = 0;
		w2x = x+.9*W/2-wr;
		w2y = 0;

		px = x + L*sin(th);
		py = y + L*cos(th);

		plot([-10 10],[0 0],'k','LineWidth',1) % ground
		hold on
		rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[0.8 0.8 0.8]); % cart
		rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[0 0 0]) % wheels
		rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[0 0 0])

		plot([x px],[y py],'k','LineWidth',4)  % pole

		rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[1 0.1 .1])
		xlim([-5 5]);
		ylim([-2 3.5]);
		% set(gcf,'Position',[100 500 1000 600])
		drawnow
		hold off
		pause(.01);
	end
end

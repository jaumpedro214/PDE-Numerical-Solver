output_sol = load('output_solution.out');
plot_info = load('plot_info.out');
output_mash = load('output_mash.out');

nodes_x = plot_info(1,1);
nodes_y = plot_info(2,1);
ts = plot_info(3,1);
nn = nodes_x*nodes_y;

zmin = min( output_sol(:,1) );
zmax = max( output_sol(:,1) );

for i = 1:1:ts;
  sol = output_sol( (1+(i-1)*nn):1:(i)*nn , 1 );
  sol = reshape(sol,nodes_x,nodes_y);
  x = output_mash(1:nodes_x);
  y = output_mash(nodes_x+1:nodes_y+nodes_x);
  
  surf(x,y,sol)
  zlim([ zmin zmax ])
  title("","interpreter", "tex");
  
  if( mod(i,10) == 0 )
    %print( ["plot_" int2str(i)], "-dpng" );
  endif
  pause(0.01)
  
 endfor
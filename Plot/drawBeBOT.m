clear all
close all

x = [1 2 3 4 5 6 7];
y = [1 2 3 4 5 6 7];
r_trn = linspace(1,length(x),length(x));
tN = linspace(0,1,length(x));
time = linspace(0,1,1000);

plot(BernsteinPoly(x,time),BernsteinPoly(y,time),'Linewidth',2,'Color','b'); hold on
plot(x, y, 'ok', 'MarkerSize', 10, 'Color','b');
plot([x(1); x(2)], [y(1); y(2)], '--','Color','r')
plot([x(end-1); x(end)], [y(end-1); y(end)], '--','Color','r')
grid on; 
strValues = strtrim(cellstr(num2str(r_trn(:))));
text(x,y,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8 );
axis equal
grid on

while 1

    select = ginput(1);


    for i = 1 : length(x)
    
        if norm(select-[x(i) y(i)]) < 0.1
    
            sel = i;
            plot(x(i), y(i), 'ok', 'MarkerSize', 10, 'Color','g'); 

            p = ginput(1);
            hold off
            
            x(i) = p(1); y(i) = p(2); 
            
            plot(BernsteinPoly(x,time),BernsteinPoly(y,time),'Linewidth',2,'Color','b'); hold on
            plot(x, y, 'ok', 'MarkerSize', 10, 'Color','b'); 
            plot([x(1); x(2)], [y(1); y(2)], '--','Color','r')
            plot([x(end-1); x(end)], [y(end-1); y(end)], '--','Color','r')
            grid on;
            strValues = strtrim(cellstr(num2str(r_trn(:))));
            text(x,y,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8 );
            axis equal
            grid on
    
        end

    end 

end
% 
% for k=1:5
%     p=ginput(1);
%     I(round(p(1,2)),round(p(1,1)))=100;
%     set(h,'CData',I);
%     drawnow;
% end
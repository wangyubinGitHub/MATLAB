function [result] = NUM_interpolate_sinc(y,nx,x,maxDepth)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

midleft = floor(x);
midright = midleft + 1;
mm=0.0;

if (nx < 1) 
	result= -1; 
    return;
end
                
if (x > nx) 
	result =  y(nx);
    return;
end
if (x < 1) 
    result = y(1);
    return;
end
if (x == midleft) 
    result = y(midleft);
    return;
end
if (maxDepth > midright - 1) 
    maxDepth = midright - 1;
end
if (maxDepth > nx - midleft) 
    maxDepth = nx - midleft;
end
if (maxDepth <= 0) 
    result = y(floor (x + 0.5)); 
    return;
end
% if (maxDepth == 1)    
%     result = y(midleft) + (x - midleft) * (y(midright) - y(midleft)); 
%     return;
% end
% if (maxDepth == 2) 
%     yl = y (midleft);
%     yr = y (midright); 
% 	dyl = 0.5 * (yr - y(midleft - 1));
%     dyr = 0.5 * (y(midright + 1) - yl);
%     fil = x - midleft;
%     fir = midright - x; 
% 	result = yl * fir + yr * fil - fil * fir * (0.5 * (dyr - dyl) + (fil - 0.5) * (dyl + dyr - 2 * (yr - yl)));
%     return;
% end
left = midright - maxDepth;
right = midleft + maxDepth;
    
a = pi * (x - midleft);%
halfsina = 0.5 * sin (a);
aa = a / (x - left + 1);
daa = pi / (x - left + 1);

for ix = left:midleft
        iix = midleft + left - ix;
        d = halfsina / a * (1.0 + cos (aa));
		mm = mm + y(iix) * d;
		a = a + pi;
		aa = aa + daa;
		halfsina = - halfsina;
end

	a = pi * (midright - x);
	halfsina = 0.5 * sin (a);
	aa = a / (right - x + 1);
	daa = pi / (right - x + 1); 
for ix = midright: right
		d = halfsina / a * (1.0 + cos (aa));
		mm = mm + y(ix) * d;
		a = a + pi;
		aa = aa + daa;
		halfsina = - halfsina;
end

result = mm;
return;
end


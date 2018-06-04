function [ymid,xmid] = NUMimproveMaximum (y, nx, ixmid,interpolation)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if (ixmid <= 1) 
    xmid = 1; 
    ymid = y(1);
    return ;
end

if (ixmid >= nx)
    xmid = nx; 
    ymid = y(nx);
    return;
end


if (interpolation <= 0)
    xmid = ixmid;
    ymid = y(ixmid);
    return;  
end


if (interpolation == 1) 
    dy = 0.5 * (y(ixmid + 1) - y(ixmid - 1));
	d2y = 2 * y(ixmid) - y(ixmid - 1) - y(ixmid + 1);
	xmid = ixmid + dy / d2y;
	ymid =  y(ixmid) + 0.5 * dy * dy / d2y;
    return;
end
    if(interpolation == 3)
        [xmid,result] = NUMminimize_brent(ixmid - 1, ixmid + 1, 1e-10, y, 70, nx, 1);
    else
        [xmid,result] = NUMminimize_brent(ixmid - 1, ixmid + 1, 1e-10, y, 700, nx, 1);
    end
    ymid = - result;
    return;
end


% Determine positions of RIS elements
function pos = RISPosition(N1,N2,dist,center)                 
    d1 = (0:N1-1)-(N1-1)/2;
    d2 = (0:N2-1)-(N2-1)/2;
    pos{1} = center(1)+d1*dist;
    pos{2} = center(2)+d2*dist;
end 

function angle = checkTangent(result,x,y)
    if x<0 && y>0
        % segon quadrant
        angle = pi + result;
    elseif  x<0 && y<0
        % tercer quadrant
        angle = result + pi;
    else
        angle = result;
    end
end

function angle = checkTangent(result,num,den)
    if num > 0 && den <0
        % segon quadrant
        angle = pi + result;
    elseif  num < 0 && den <0
        % tercer quadrant
        angle = result + pi;
    elseif num < 0 && den > 0
        % quart quadrant
        angle = result + 2*pi;
    else
        angle = result;
    end
    rad2deg(angle)
end

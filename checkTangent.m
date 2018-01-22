function angle = checkTangent(result,num,den)
    if den<0 && num>0
        % segon quadrant
        angle = pi + result;
    elseif  den<0 && num<0
        % tercer quadrant
        angle = result + pi;
    else
        angle = result;
    end
end

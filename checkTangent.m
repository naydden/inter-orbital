function angle = checkTangent(result,num,den)
% Function that rectifies the angle obtained with atan() as a function
% of the numerator and the denominator
% OUTPUT  angle: rectified angle [rad]
% INPUTS  result: angle obtained with atan() [rad]
% num: numerator  den: denominator
    if den<0 && num>0
        % second quadrant
        angle = pi + result;
    elseif  den<0 && num<0
        % third quadrant
        angle = result + pi;
    else
        angle = result;
    end
end

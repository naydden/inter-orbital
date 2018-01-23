function JDN=JulianDate(year,month,day)

a = floor((14 - month) / 12);
y = year + 4800 - a;
m = month + 12*a - 3;

JDN = day + floor((153*m + 2)/5) + 365*y + floor(y/4) 
      - floor(y/100) + floor(y/400) - 32045;
end
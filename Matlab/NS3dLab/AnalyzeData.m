if(mod(t,HowOftenSave)==0)
    bmax = 1;
    bmin = 0;
    binvals = ( linspace(bmin, bmax, 260));
    HIST(counter,:) = hist(T(find(T>=bmin & T <= bmax)), binvals); counter = counter + 1;
end
    
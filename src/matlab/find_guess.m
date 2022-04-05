% Find point with different sign as f(l) and f(r) by guessing
% f(l) and f(r) should have the same sign
function guess = find_guess(f, l, r)
    if sign(f(l)) ~= sign(f(r))
        error('Using this function wrong\n');
    end

    % Try midway point
    if sign(f((l+r)/2)) ~= sign(f(l))
        guess = (l+r)/2;
        return;
    end

    % Try quarterway point
    if sign(f((l+r)/4)) ~= sign(f(l))
        guess = (l+r)/4;
        return;
    end

    maxtry = 10000; % Try 10000 guesses
    it = 0;
    guess = rand(1)*(r-l)+l;
    while(sign(f(guess)) == sign(f(l)))
        guess = rand(1)*(r-l)+l;
        it = it + 1;
        if(it >= maxtry)
            fprintf('Could not find a good guess\n');
            guess = -1;
            return;
        end
    end
end
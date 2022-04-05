function [mu, fmu, Ej, fEj, Estar, fEstar, N, Ez, M11,M12,traceM, m, a] = ist_spectra(Emin, Emax, numE, xz, u)

Ez = linspace(Emin, Emax, numE);
kz = sqrt(Ez);
M11 = ones(length(kz), 1);  
M12 = ones(length(kz), 1);
S = zeros(length(kz), 1);
traceM = ones(length(kz),1);
M11prev = 1;
M11signflipEz = [];
M11signflipEzidx = [];
g = @(E) M12f(E,xz,u,1);


[M,~] = mmat_mex(Ez(1), xz, u);
M11(1) = M(1,1);
M12(1) = M(1,2);
traceM(1) = 0.5*(M(1,1)+M(2,2));
M11prev = M11(1);
i = 2;
for E = Ez(2:end)
    [M,~] = mmat_mex(E, xz, u);
    M11(i) = M(1,1);
    M12(i) = M(1,2);
    traceM(i) = 0.5*(M(1,1)+M(2,2));
    delS = abs(sign(M11(i))-sign(M11prev))/2;
    if(delS > 0) 
        M11signflipEz = [M11signflipEz E];
        M11signflipEzidx = [M11signflipEzidx i];
    end
    if (isnan(delS))
        S(i) = 0;
    else
        S(i) = S(i-1) + delS;
    end
    M11prev = M11(i);
    i = i+1;
end

lastnan = sum(isnan(M11));
M11 = M11(lastnan+1:end);
M12 = M12(lastnan+1:end);
S   = S(lastnan+1:end);
traceM = traceM(lastnan+1:end);
Ez  = Ez(lastnan+1:end);
Emin = Ez(1);

M12signflips = sort([find(M12(1:end-1)>0 & M12(2:end) < 0); find(M12(1:end-1)<0 & M12(2:end) > 0)]);
M12signfliplocs = Ez(M12signflips);

% Step I: Search for Estar which seperate degrees of freedom (N of these)
N = length(M11signflipEz);
f = @(E) M11f(E,xz,u,1);
M11fabs = @(E) abs(M11f(E,xz,u,1));
Estar = zeros(N,1);
fEstar = zeros(N,1);
tolX = 1e-300; % For custom bisection search
opt = optimset('FunValCheck','on','TolX',eps, 'TolFun', 1e-5);%, 'Display','iter-detailed');
[Ezero,fval,exitflag,output] = fzero(f,M11signflipEz(1),opt);
if exitflag < 1
    error('Something went wrong finding zero\n');
end
Estar(1) = Ezero;
fEstar(1) = fval;

% Check for later to move Estar to left by epsilon if it's a bit too far to
% the right for the Step II pass to find the first zero of g()
if (g(Estar(1)) < 0) % Can't be true ever - must be numerical problem
    Estar(1) = Estar(1) - 2*eps(Estar(1));
    fEstar(1) = f(Estar(1));
end

for i = 2:N-1
    % TODO: Replace with fmincon to find zeros which only touch but don't
    % change sign
    lb = (M11signflipEz(i-1)+M11signflipEz(i))/2;
    ub = (M11signflipEz(i+1)+M11signflipEz(i))/2;
    x0 = M11signflipEz(i);
    [Ezero,fval,exitflag,output] = fmincon(M11fabs,x0,[],[],[],[],lb,ub); %fzero(f,M11signflipEz(i),opt);
    if exitflag < 1
        error('Something went wrong finding zero\n');
    end
    Estar(i) = Ezero;
    fEstar(i) = fval;
end

[Ezero,fval,exitflag,output] = fzero(f,Ez(end),opt);
if (Ezero > Emax)
    [Ezero,fval,exitflag,output] = fzero(f,[M11signflipEz(end-1)+10*eps Emax],opt);
end
Estar(end) = Ezero;
fEstar(end) = fval;

% Step II: The Auxilary Spectrum
mu = zeros(N-1,1);
fmu = zeros(N-1,1);

midx = 1;
for i = 2:N
    leftbnd = Estar(i-1);
    rghtbnd = Estar(i);
    if (i==2)
        leftbnd = Emin;
    end
    
    signdiff = sign(g(leftbnd)) ~= sign(g(rghtbnd));
    if (signdiff) 
        [Ezero,fval,~,~] = BetterIntervalSearch_M12(xz, u, leftbnd, rghtbnd, tolX);
        mu(midx) = Ezero;
        fmu(midx) = fval;
        midx = midx + 1;
    else % Possibly more than one zero crossing
        guess = find_guess(g, leftbnd, rghtbnd);
        if sign(g(leftbnd)) == sign(g(guess))
            fprintf('something wrong here, bad guess\n');
        end

        if guess == -1
            continue; % Skip this loop - zero will probably be picked up on the next one
        else
            [Ezero,fval,~,~] = BetterIntervalSearch_M12(xz, u, leftbnd, guess, tolX);
        end

        if exitflag < 0
            error('could not find zero with double crossing\n');
        end

        if Ezero > Estar(i) || Ezero < leftbnd % outside interval
            continue
        else
            mu(midx) = Ezero;
            fmu(midx) = fval;
            midx = midx + 1;

            % Double crossing
            if sign(g(guess)) == sign(g(rghtbnd))
                fprintf('something wrong here\n');
            end
            if guess == -1
                continue;
            else
                [Ezero,fval,~,~] = BetterIntervalSearch_M12(xz, u, guess, rghtbnd, tolX);
            end
            mu(midx) = Ezero;
            fmu(midx) = fval;
            midx = midx + 1;
        end
    end
end

% Step III: The Main Spectrum
Ej = zeros(2*N-1,1);
fEj = zeros(2*N-1,1);
pwr = 1;
tolX = 1e-20;
[Ezero,fval,~,~] = BetterIntervalSearch_PlusOne(xz, u, Ez(1), mu(1), tolX);
%l,Ez(1)-extra_rng,mu(1)+extra_rng,opts);
Ej(1) = Ezero;
fEj(1) = (fval)^(1/pwr)-1;
[Ezero,fval,~,~] = BetterIntervalSearch_MinusOne(xz, u, Ez(1), mu(1), 1e-20);
Ej(2) = Ezero;
fEj(2) = (fval)^(1/pwr)+1;

for i = 1:N-2
    
    loffset = 0.0;
    % Move away from left boundary
    if (abs(fmu(i)) > 1) % Very steep boundary
        loffset = 2*eps(mu(i));
    end
    [Ezero,fval,~,~] = BetterIntervalSearch_MinusOne(xz, u, mu(i)+loffset, mu(i+1), tolX);
    %fminbnd(h,mu(i)-extra_rng,mu(i+1)+extra_rng,opts);
%     if(exitflag<0)
%         continue
%     else
        Ej(2*i+1) = Ezero;
        fEj(2*i+1) = (fval)^(1/pwr)+1;
%     end
    [Ezero,fval,~] = BetterIntervalSearch_PlusOne(xz, u, mu(i)+loffset, mu(i+1), tolX);
    %[Ezero,fval,exitflag,output] = fminbnd(l,mu(i)-extra_rng,mu(i+1)+extra_rng,opts);
%     if(exitflag<0)
%         continue
%     else
        Ej(2*i+2) = Ezero;
        fEj(2*i+2) = (fval)^(1/pwr)-1;
%     end
end

Ej_up_till_now = sort(Ej(1:end-1));
[M,~] = mmat_mex(Ej_up_till_now(end), xz, u);
traceatlastEj = trace(M)/2;

if( traceatlastEj < 0) % tr M @ last Ej evaluated
    [Ezero,fval,~,~] = BetterIntervalSearch_PlusOne(xz, u, mu(end), Ez(end), tolX);
    %[Ezero,fval,exitflag,~] = fminbnd(l,mu(end)-extra_rng, Ez(end),opts);
    %if(exitflag<0)
    %    fprintf('Missed a Ezero when searching for last Ej\n');
    %    fval = -2;
    %    Ezero = -1;
    %else
    fval = (fval)^(1/pwr) -1;
    %end
    
else
    %[Ezero,fval,exitflag,~] = fminbnd(h,mu(end)-extra_rng, Ez(end),opts);
    [Ezero,fval,~,~] = BetterIntervalSearch_MinusOne(xz, u, mu(end), Ez(end), tolX);
    %if(exitflag<0)
    %    fprintf('Missed a Ezero when searching for last Ej\n');
    %    fval = -2;
    %    Ezero = -1;
    %else
        fval = (fval)^(1/pwr) +1;
    %end
end
Ej(end) = Ezero;
fEj(end) = fval;

[Ej, SI] = sort(Ej);
fEj = fEj(SI);

m = zeros(N-1,1);
a = zeros(N-1,1);
Ejsort = sort(Ej);
for j = 1:1:N-1
    m(j) = (Ejsort(2*j+1) - Ejsort(2*j))/(Ejsort(2*j+1) - Ejsort(2*j-1));
    a(j) = Ejsort(2*j+1) - Ejsort(2*j);
end

end
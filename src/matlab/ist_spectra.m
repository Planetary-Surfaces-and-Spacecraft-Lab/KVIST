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

% Which method do I use to compute mu?


%if N-1 < length(M12signflips) % Lost an Estar, use alt method
%     N = length(M12signflips)+1;
%     mu = zeros(length(M12signflips),1);
%     fmu = zeros(length(M12signflips),1);
%     for i = 1:length(M12signflips)
%         [Ezero,fval,exitflag,output] = fzero(g,M12signfliplocs(i),opt);
%         if exitflag < 0
%             error('Something wrong at mu step');
%         end
%         mu(i) = Ezero;
%         fmu(i) = fval;
%     end
%else % Number of Estar is correct, use Estar to set bisection method
    mu = zeros(N-1,1);
    fmu = zeros(N-1,1);
    skipnext = 0;

    for i = 2:N
        if skipnext
            skipnext = 0;
            continue
        end

        % One zero crossing
         %fprintf('i=%i\ng(Estar(i-1)) = %g\ng(Estar(i)) = %g\n', ...
         %        i, sign(g(Estar(i-1))), sign(g(Estar(i))));
         if sign(g(Estar(i-1))) ~= sign(g(Estar(i)))
          %   fprintf('Signs not same\n');
           %  fprintf('g(%g) = %g\ng(%g) = %g\n', Estar(i-1), g(Estar(i-1)), ...
            %                                     Estar(i), g(Estar(i)));
            [Ezero,fval,~,~] = BetterIntervalSearch_M12(xz, u, Estar(i-1), Estar(i), tolX);
            %[Ezero,fval,exitflag,output] = fzero(g,[Estar(i-1) Estar(i)],opt);
            mu(i-1) = Ezero;
            fmu(i-1) = fval;
        else % Possibly more than one zero crossing
            %fprintf('Signs same at i = %i\n', i);
            guess = Estar(i-1)+(Estar(i)-Estar(i-1))/4;
            [Ezero,fval,exitflag,output] = fzero(g,guess,opt);
            if exitflag < 0
                error('could not find zero with double crossing\n');
            end

            if Ezero > Estar(i) || Ezero < Estar(i-1) % outside interval
                continue
            else
                mu(i-1) = Ezero;
                fmu(i-1) = fval;

                fprintf('Double crossing\n');
                guess = Estar(i-1)+3*(Estar(i)-Estar(i-1))/4;
                [Ezero,fval,exitflag,output] = fzero(g,guess,opt);
                mu(i) = Ezero;
                fmu(i) = fval;
                skipnext = 1;
            end
        end
    end
%end


% Last one
% Eguess = (Estar(end)+Ez(end))/2;
% [Ezero,fval,exitflag,output] = fzero(g,Eguess);
% mu(end) = Ezero;
% [M,~] = mmat_mex(mu(end), xz, u);
% fmu(end) = M(1,2);

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
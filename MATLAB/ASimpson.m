function [R,C,fC] = ASimpson(f,a,b,tol,level,level_max,par,fpar)
    level = level+1;
    h = b-a;
    c = (a+b)/2;
        d = (a+c)/2;
    e = (b+c)/2;
    C = par; 
    
    if isempty(fpar)
        fa = f(a);
        fb = f(b);
        fC = containers.Map([a,b],[fa,fb]); 
    else
        fC = fpar; 
    end
    
    if isKey(fC,c)
        fc = fC(c);
    else
        fc = f(c); C= [C,c]; fC(c) = fc;
    end
    if isKey(fC,d)
        fd = values(fC,d);
    else
        fd = f(d); C = [C,d]; fC(d) = fd;
    end
    if isKey(fC,e)
        fe = fC(e);
    else
        fe = f(e); C = [C,e]; fC(e) = fe;
    end
    if isKey(fC,a)
        fa = fC(a);
    else
        fa = f(a); C = [C,a]; fC(a) = fa;
    end
    if isKey(fC,b)
        fb = fC(b);
    else
        fb = f(b); C = [C,b]; fC(b) = fb;
    end
    
    one_simp = h*(fa+4*fc+ fb)/6;
    two_simp = h*(fa+4*fd+2*fc+4*fe+ fb)/12;
    if(level >= level_max)
        R = two_simp;
    elseif(abs(two_simp-one_simp) < 15*tol)
        R = two_simp+(two_simp-one_simp)/15;
    else 
        [l_simp,p,fp] = ASimpson(f,a,c,tol/2,level,level_max,C,fC);
        [r_simp,q,fq] = ASimpson(f,c,b,tol/2,level,level_max,p,fp);
        fC = fq;
        C = q;
        R = l_simp+r_simp;
    end
end


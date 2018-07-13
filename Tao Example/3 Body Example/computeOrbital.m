function [ res ] = computeOrbital( qq, vv )
% Cartesian -> osculating orbital elements
% qq,vv ~ 2*1 or 3*1, m: mass of the planet; Gcst & m0: inherited globally
% res = [a;e;i;omega;Omega;E]
global Gcst m0

    threshold=1e-6; % equatorial orbit threshold

    if length(qq)==2
        qq=[qq;0];
        vv=[vv;0];
    end
    mu=Gcst*m0;

    % Cartesian -> osculating orbital elements: a,e,I,omega,Omega,Eanomaly
        % see https://en.wikibooks.org/wiki/Astrodynamics/Classical_Orbit_Elements and https://en.wikipedia.org/wiki/Argument_of_periapsis
    a=1/(2/norm(qq)-vv'*vv/mu);
    hvec=cross(qq,vv);
    nvec=cross([0;0;1],hvec);
    c=norm(hvec);
    evec=(vv'*vv/mu-1/norm(qq))*qq-qq'*vv/mu*vv;
    e=norm(evec);
    I=acos(hvec(3)/c);
    
	if (abs(sin(I))>threshold)          % non equatorial orbit, i.e., I>0
        Omega=real(acos(nvec(1)/norm(nvec)));
        if nvec(2)<0
            Omega=2*pi-Omega;
        end
        if (e>threshold)
            omega=real(acos(nvec'*evec/norm(nvec)/norm(evec)));
            if evec(3)<0                % noncircular orbit
                omega=2*pi-omega;
            end
        else                            % circular orbit, adopt convention of omega=0
            omega=0;
        end
    else                                % equatorial orbit; node not defined, adopt convention of Omega=0
        Omega=0;
        if (e>threshold)                % noncircular orbit
            omega=atan2(evec(2),evec(1));
            if hvec(3)<0
                omega=2*pi-omega;
            end
        else                            % circular orbit, adopt convention of omega=0
            omega=0;
        end
    end
        
    temp=[cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1]*[1 0 0; 0 cos(I) sin(I); 0 -sin(I) cos(I)]*[cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1]*qq;
    Eanomaly=atan2(temp(2)/a/sqrt(1-e^2),temp(1)/a+e);

    res=[a;e;I;omega;Omega;Eanomaly];
end
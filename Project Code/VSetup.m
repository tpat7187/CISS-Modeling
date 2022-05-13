%V Setup

V1 = [];
V2 = [];
V3 = [];

for m  = 1:total_sites
    phi(m,1) = (m-1)*2*pi/ions_per_lap;
    r(m,1) = a*cos(phi(m));
    r(m,2) = a*sin(phi(m));
    r(m,3) = (m-1)*c/(lap*ions_per_lap-1);
end

%general cases
for k = 1:total_sites-2
    dm1 = (r(k,:)-r(k+1,:))/norm(r(k,:)-r(k+1,:));
    dm2 = (r(k,:)-r(k+2,:))/norm(r(k,:)-r(k+2,:));
    V = cross(dm1,dm2);
    V1(k,k+2) = V(1);
    V2(k,k+2) = V(2);
    V3(k,k+2) = V(3);
end

%boundary conditions (10->2)
dm1 = (r(total_sites,:) - r(1,:))/norm(r(total_sites,:)-r(1,:));
dm2 = (r(total_sites,:) - r(2,:))/norm(r(total_sites,:)-r(2,:));
V = cross(dm1,dm2);
V1(total_sites,2) = V(1);
V2(total_sites,2) = V(2);
V3(total_sites,2) = V(3);

%boundary conditions (9->1)
dm1 = (r(total_sites-1,:) - r(total_sites,:))/norm(r(total_sites-1,:)-r(total_sites,:));
dm2 = (r(total_sites-1,:) - r(1,:))/norm(r(total_sites,:)-r(1,:));
V = cross(dm1,dm2);
V1(total_sites-1,1) = V(1);
V2(total_sites-1,1) = V(2);
V3(total_sites-1,1) = V(3);


% FK model in 0d

% some found values, 100 particles, 5 iterations
%tr =  33.8854, tsi =  29.4032, twp=836.6229, td =   0.2674 error=57.10
%tr =  33.9775, tsi =  29.6192, twp=857.9661, td =   0.2653 error=190.53
%tr =  32.7071, tsi =  28.3530, twp=850.5890, td =   0.2088 error=35.85
%tr =  34.2110, tsi =  29.7101, twp=867.0858, td =   0.2739 error=79.96
%tr =  34.0868, tsi =  29.5414, twp=851.1914, td =   0.2550 error=97.61
% CORRECT: 
%tr =  33.33,   tsi =  29.0,    twp=870,      td =   0.25

% parameter values
%tvp=3.33;
%tv1m=19.6;
%tv2m=1250.;
%twp=870.0;
% twm=41.0;
%td=0.25;
% to=12.5;
%tr=33.33;
%tsi=29.0;
% xk=10.0;
% ucsi=0.85;
% uc=0.13;
% uv=0.04;
stimmag=0.1;

usave0 = load('../data/zebrafish_onecl.txt'); %400ms
usave0 = usave0/max(usave0(:));

tr_true = 33.33; tsi_true = 29; twp_true = 870; td_true = 0.25;
tvp_true=3.33; tv1m_true=19.6; tv2m_true=1250.;
twm_true=40.0; to_true=12.5; xk_true=10.0; 
ucsi_true=0.85; uc_true=0.13; uv_true=0.04;
% particle params
nparticles = 1000;
trmin = 25;         trmax = 200;
tsimin = 10;        tsimax = 300;
twpmin = 50;        twpmax = 500;
tdmin = 0.15;       tdmax = 0.4;
tvpmin = 1;         tvpmax = 20;
tv1mmin = 10;       tv1mmax = 50;
tv2mmin = 500;      tv2mmax = 1500;
twmmin = 5;         twmmax = 100;
tomin = 5;          tomax = 50;
xkmin = 1;          xkmax = 15;
ucsimin = 0.2;      ucsimax = 0.9;
ucmin = 0.1;        ucmax = 0.25;
uvmin = 0.005;      uvmax = 0.05;

% original; note crazy tv2mmax 
% trmin = 25;         trmax = 36;
% tsimin = 25;        tsimax = 30;
% twpmin = 500;       twpmax = 1000;
% tdmin = 0.15;       tdmax = 0.3;
% tvpmin = 1;         tvpmax = 10;
% tv1mmin = 10;       tv1mmax = 50;
% tv2mmin = 500;      tv2mmax = 15000;
% twmmin = 10;        twmmax = 100;
% tomin = 5;          tomax = 50;
% xkmin = 1;          xkmax = 100;
% ucsimin = 0.3;      ucsimax = 0.95;
% ucmin = 0.1;        ucmax = 0.25;
% uvmin = 0.09;       uvmax = 0.01;

% initialize particles
tr_particle = trmin+(trmax-trmin)*rand(1,nparticles);
tsi_particle = tsimin+(tsimax-tsimin)*rand(1,nparticles);
twp_particle = twpmin+(twpmax-twpmin)*rand(1,nparticles);
td_particle = tdmin+(tdmax-tdmin)*rand(1,nparticles);
tvp_particle = tvpmin+(tvpmax-tvpmin)*rand(1,nparticles);
tv1m_particle = tv1mmin+(tv1mmax-tv1mmin)*rand(1,nparticles);
tv2m_particle = tv2mmin+(tv2mmax-tv2mmin)*rand(1,nparticles);
twm_particle = twmmin+(twmmax-twmmin)*rand(1,nparticles);
to_particle = tomin+(tomax-tomin)*rand(1,nparticles);
xk_particle = xkmin+(xkmax-xkmin)*rand(1,nparticles);
ucsi_particle = ucsimin+(ucsimax-ucsimin)*rand(1,nparticles);
uc_particle = ucmin+(ucmax-ucmin)*rand(1,nparticles);
uv_particle = uvmin+(uvmax-uvmin)*rand(1,nparticles);

% initialize directions (not sure what to set these to!
vinit_scale = 0.00005;
tr_v = vinit_scale*rand(1,nparticles);
tsi_v = vinit_scale*rand(1,nparticles);
twp_v = vinit_scale*rand(1,nparticles);
td_v = vinit_scale*rand(1,nparticles);
tvp_v = vinit_scale*rand(1,nparticles);
tv1m_v = vinit_scale*rand(1,nparticles);
tv2m_v = vinit_scale*rand(1,nparticles);
twm_v = vinit_scale*rand(1,nparticles);
to_v = vinit_scale*rand(1,nparticles);
xk_v = vinit_scale*rand(1,nparticles);
ucsi_v = vinit_scale*rand(1,nparticles);
uc_v = vinit_scale*rand(1,nparticles);
uv_v = vinit_scale*rand(1,nparticles);
test = rand(1,nparticles);

% set up best storage for each particle
tr_best = tr_particle;
tsi_best = tsi_particle;
twp_best = twp_particle;
td_best = td_particle;
tvp_best = td_particle;
tv1m_best = td_particle;
tv2m_best = td_particle;
twm_best = twm_particle;
to_best = to_particle;
xk_best = xk_particle;
ucsi_best = ucsi_particle;
uc_best = uc_particle;
uv_best = uv_particle;run_simulations_solver
dist_best = 100000*ones(size(tr_particle));

% numerical parameters
dt = 0.02; % size of time step in ms
period=400.0;
nperiod = ceil(period/dt);
nbeats=1;
endtime=nbeats*period;
nsteps = ceil(endtime/dt);
stimdur= 0.5;
nstimdur = ceil(stimdur/dt);
t = 0:dt:endtime;

v0=1;
w0=1;

% first, solve the model
particle_iterations = 8
best_error = zeros(particle_iterations,1);

debug=0;

for k=1:particle_iterations
    fprintf("Current avg tr velocity:\t%f\n",mean(tr_v));
    fprintf('iteration %d\n',k)
    if(debug)
        fprintf('***** k = %d\n',k)
        [tr_particle; tsi_particle; twp_particle; td_particle; tvp_particle; tv1m_particle; tv2m_particle; twm_particle; to_particle; xk_particle; ucsi_particle; uc_particle; uv_particle]
    end
    % first, solve the model

    % initial values for state variables
    u=zeros(1,nparticles);  
    v=v0*ones(1,nparticles);
    w=w0*ones(1,nparticles);
    nsave=length(usave0);
    usave = zeros(nsave,nparticles);
    vsave = usave; wsave=usave;
    tsave = (1:nsave)';
    ktime=0;
%    usave(1,:) = u;
%    vsave(1,:) = v;
%    wsave(1,:) = w;

    for ntime=1:nsteps
        
%         p=0; % heaviside functions
%         q=0;
%         if(u>=uc)
%             p=1;
%         end
%         if(u>=uv)
%             q=1.;
%         end
        p=(u>=uc_particle);
        q=(u>=uv_particle);
        dv=(1.-p).*(1.-v)./((1.-q).*tv1m_particle+tv2m_particle.*q) - p.*v./tvp_particle;
        dw=(1.-p).*(1.-w)./twm_particle - p.*w./twp_particle;
        v=v+dt*dv; % solving/updating v and w
        w=w+dt*dw;
        jfi=-v.*p.*(u-uc_particle).*(1.-u)./td_particle;
        jso=u.*(1.-p)./to_particle+p./tr_particle;
        jsi=-w.*(1.+tanh(xk_particle.*(u-ucsi_particle)))./(2.*tsi_particle);
        
        % stimulus
        istim=0;
        if((mod(ntime,period/dt)<7/dt)&&(mod(ntime,period/dt)>2/dt))
            istim=stimmag;
        end
        
        u=u - (jfi+jso+jsi-istim)*dt; % updating/integrating u
        
        if(mod(ntime,round(1/dt))==0)
            ktime=ktime+1;
            usave(ktime,:) = u;
            vsave(ktime,:) = v;
            wsave(ktime,:) = w;
        end
    end
    
%    plot(t,usave,'k',t,vsave,'r','linewidth',2)
    figure(1),subplot(2,1,1)
    plot(tsave,usave0,'k',tsave,usave,'r','linewidth',2),ylim([-0.1 1.6])
    xlabel('Time'),title(['iteration = ' num2str(k)])
    drawnow
%    q=input('press any key to continue')

    % second, measure distance as SSE for each particle
    dist=min(100000,sum((usave0-usave).^2));
    if(debug)
        dist
    end
    
    % third, identify local_best and global_best
    % compare each particle's current and best distance
    % if current distance is smaller, replace best params and best distance
    % then, after updating all the local bests, 
    % find index of smallest local best to identify global best
    prev_dist_best = dist_best;
    dist_best=(dist<dist_best).*dist+(dist>=dist_best).*dist_best;
    if(debug)
        dist_best
    end
    % look for any problems
    change = prev_dist_best - dist_best;
    uhoh = find(change<0);
    if(~isempty(uhoh))
        fprintf('help! there is a problem in iteration %d\n',k)
    end

    figure(1),subplot(2,1,2)
    plot(prev_dist_best,'ko'),hold on, plot(dist_best,'ro'),hold off
    
    tr_temp=(dist<prev_dist_best).*tr_particle+(dist>=prev_dist_best).*tr_best;
    tsi_temp=(dist<prev_dist_best).*tsi_particle+(dist>=prev_dist_best).*tsi_best;
    twp_temp=(dist<prev_dist_best).*twp_particle+(dist>=prev_dist_best).*twp_best;
    td_temp=(dist<prev_dist_best).*td_particle+(dist>=prev_dist_best).*td_best;
    tvp_temp=(dist<prev_dist_best).*tvp_particle+(dist>=prev_dist_best).*tvp_best;
    tv1m_temp=(dist<prev_dist_best).*tv1m_particle+(dist>=prev_dist_best).*tv1m_best;
    tv2m_temp=(dist<prev_dist_best).*tv2m_particle+(dist>=prev_dist_best).*tv2m_best;
    twm_temp=(dist<prev_dist_best).*twm_particle+(dist>=prev_dist_best).*twm_best;
    to_temp=(dist<prev_dist_best).*to_particle+(dist>=prev_dist_best).*to_best;
    xk_temp=(dist<prev_dist_best).*xk_particle+(dist>=prev_dist_best).*xk_best;
    ucsi_temp=(dist<prev_dist_best).*ucsi_particle+(dist>=prev_dist_best).*ucsi_best;
    uc_temp=(dist<prev_dist_best).*uc_particle+(dist>=prev_dist_best).*uc_best;
    uv_temp=(dist<prev_dist_best).*uv_particle+(dist>=prev_dist_best).*uv_best;
    
    tr_best=tr_temp;
    tsi_best=tsi_temp;
    twp_best=twp_temp;
    td_best=td_temp;
    tvp_best=tvp_temp;
    tv1m_best=tv1m_temp;
    tv2m_best=tv2m_temp;
    twm_best=twm_temp;
    to_best=to_temp;
    xk_best=xk_temp;
    ucsi_best=ucsi_temp;
    uc_best=uc_temp;
    uv_best=uv_temp;

    if(debug)
        [tr_best; tsi_best; twp_best; td_best; tvp_best; tv1m_best; tv2m_best; twm_best; to_best; xk_best; ucsi_best; uc_best; uv_best]
    end
    [dist_global,ind_global]=min(dist_best(:));
    if(k>1)
        if(best_error(k-1)<dist_best(ind_global))
            fprintf('there is a problem in iteration %d\n',k)
            fprintf(' prev min dist = %f8.2\n',best_error(k-1))
            fprintf(' new min dist = %f8.2\n',dist_best(ind_global))
            fprintf(' dist_global = %f8.2\n',dist_global)
        end
    end
    tr_global=tr_best(ind_global);
    tsi_global=tsi_best(ind_global);
    twp_global=twp_best(ind_global);
    td_global=td_best(ind_global);
    tvp_global=tvp_best(ind_global);
    tv1m_global=tv1m_best(ind_global);
    tv2m_global=tv2m_best(ind_global);
    twm_global=twm_best(ind_global);
    to_global=to_best(ind_global);
    xk_global=xk_best(ind_global);
    ucsi_global=ucsi_best(ind_global);
    uc_global=uc_best(ind_global);
    uv_global=uv_best(ind_global);
% only update if the best changed on this iteration
    if(min(dist(:))<min(prev_dist_best(:)))
        u_global = usave(:,ind_global);
    end
    best_error(k)=dist_best(ind_global);
    fprintf('iteration %d, min error = %8.2f\n',k,dist_global); 
    fprintf('  tr=%8.3f tsi=%8.3f twp=%8.3f td=%8.3f\n',tr_global,tsi_global,twp_global,td_global);
    fprintf('  tvp=%8.3f tv1m=%8.3f tv2m=%8.3f\n',tvp_global,tv1m_global,tv2m_global);
    fprintf('  twm=%8.3f to=%8.3f xk=%8.3f\n',twm_global,to_global,xk_global);
    fprintf('  ucsi=%8.3f uc=%8.3f uv=%8.3f\n',ucsi_global,uc_global,uv_global);

    if(debug)
        [tr_global; tsi_global; twp_global; td_global; tvp_global; tv1m_global; tv2m_global; twm_global; to_global; xk_global; ucsi_global; uc_global; uv_global]
    end
    
%     figure(1)
%     subplot(2,1,1)
%     hold on,plot(t,usave(:,ind_best))
    
    % fourth, update each particle
    phi_local=2.05; phi_global=2.05; phi=phi_local+phi_global;
    chi = 0.05 * 2/(phi-2+sqrt(phi*(phi-4)));
    chi_tr = chi*(trmax+trmin)/2;
    chi_tsi = chi*(tsimax+tsimin)/2;
    chi_twp = chi*(twpmax+twpmin)/2;
    chi_td = chi*(tdmax+tdmin)/2;
    chi_tvp = chi*(tvpmax+tvpmin)/2;
    chi_tv1m = chi*(tv1mmax+tv1mmin)/2;
    chi_tv2m = chi*(tv2mmax+tv2mmin)/2;
    chi_twm = chi*(twmmax+twmmin)/2;
    chi_to = chi*(tomax+tomin)/2;
    chi_xk = chi*(xkmax+xkmin)/2;
    chi_ucsi = chi*(ucsimax+ucsimin)/2;
    chi_uc = chi*(ucmax+ucmin)/2;
    chi_uv = chi*(uvmax+uvmin)/2;
    tr_v = chi_tr * (tr_v + phi_local*rand(size(tr_particle)) .* (tr_best-tr_particle) ...
        + phi_global*rand(size(tr_particle)) .* (tr_global-tr_particle));
    tsi_v = chi_tsi * (tsi_v + phi_local*rand(size(tsi_particle)) .* (tsi_best-tsi_particle) ...
        + phi_global*rand(size(tsi_particle)) .* (tsi_global-tsi_particle));
    twp_v = chi_twp * (twp_v + phi_local*rand(size(twp_particle)) .* (twp_best-twp_particle) ...
        + phi_global*rand(size(twp_particle)) .* (twp_global-twp_particle));
    td_v = chi_twp * (td_v + phi_local*rand(size(td_particle)) .* (td_best-td_particle) ...
        + phi_global*rand(size(td_particle)) .* (td_global-td_particle));
    tvp_v = chi_tvp * (tvp_v + phi_local*rand(size(tvp_particle)) .* (tvp_best-tvp_particle) ...
        + phi_global*rand(size(tvp_particle)) .* (tvp_global-tvp_particle));
    tv1m_v = chi_tv1m * (tv1m_v + phi_local*rand(size(tv1m_particle)) .* (tv1m_best-tv1m_particle) ...
        + phi_global*rand(size(tv1m_particle)) .* (tv1m_global-tv1m_particle));
    tv2m_v = chi_tv2m * (tv2m_v + phi_local*rand(size(tv2m_particle)) .* (tv2m_best-tv2m_particle) ...
        + phi_global*rand(size(tv2m_particle)) .* (tv2m_global-tv2m_particle));
    twm_v = chi_twm * (twm_v + phi_local*rand(size(twm_particle)) .* (twm_best-twm_particle) ...
        + phi_global*rand(size(twm_particle)) .* (twm_global-twm_particle));
    to_v = chi_to * (to_v + phi_local*rand(size(to_particle)) .* (to_best-to_particle) ...
        + phi_global*rand(size(to_particle)) .* (to_global-to_particle));
    xk_v = chi_xk * (xk_v + phi_local*rand(size(xk_particle)) .* (xk_best-xk_particle) ...
        + phi_global*rand(size(xk_particle)) .* (xk_global-xk_particle));
    ucsi_v = chi_ucsi * (ucsi_v + phi_local*rand(size(ucsi_particle)) .* (ucsi_best-ucsi_particle) ...
        + phi_global*rand(size(ucsi_particle)) .* (ucsi_global-ucsi_particle));
    uc_v = chi_uc * (uc_v + phi_local*rand(size(uc_particle)) .* (uc_best-uc_particle) ...
        + phi_global*rand(size(uc_particle)) .* (uc_global-uc_particle));
    uv_v = chi_uv * (uv_v + phi_local*rand(size(uv_particle)) .* (uv_best-uv_particle) ...
        + phi_global*rand(size(uv_particle)) .* (uv_global-uv_particle));
    
    tr_particle = tr_particle + tr_v;
    tsi_particle = tsi_particle + tsi_v;
    twp_particle = twp_particle + twp_v;
    td_particle = td_particle + td_v;
    tvp_particle = tvp_particle + tvp_v;
    tv1m_particle = tv1m_particle + tv1m_v;
    tv2m_particle = tv2m_particle + tv2m_v;
    twm_particle = twm_particle + twm_v;
    to_particle = to_particle + to_v;
    xk_particle = xk_particle + xk_v;
    ucsi_particle = ucsi_particle + ucsi_v;
    uc_particle = uc_particle + uc_v;
    uv_particle = uv_particle + uv_v;
    
    if(debug)
        [tr_particle; tsi_particle; twp_particle; td_particle; tvp_particle; tv1m_particle; tv2m_particle; twm_particle; to_particle; xk_particle; ucsi_particle; uc_particle; uv_particle]
    end
    
%     figure(2)
%     subplot(3,1,1)
%     plot(a_particle,'ro'),hold on,plot(amax*ones(size(a_particle)),'k'),plot(amin*ones(size(a_particle)),'k'),hold off
%     subplot(3,1,2)
%     plot(b_particle,'ro'),hold on,plot(bmax*ones(size(b_particle)),'k'),plot(bmin*ones(size(b_particle)),'k'),hold off
%     subplot(3,1,3)
%     plot(eps_particle,'ro'),hold on,plot(epsmax*ones(size(eps_particle)),'k'),plot(epsmin*ones(size(eps_particle)),'k'),hold off
    
    % check for any updated param being out of bounds
    tr_particle = (tr_particle<trmax).*tr_particle + (tr_particle>=trmax).*(trmin+0.75*(trmax-trmin)*rand(1,nparticles));
    tr_particle = (tr_particle>trmin).*tr_particle + (tr_particle<=trmin).*(trmin+0.25*(trmax-trmin)+0.75*(trmax-trmin)*rand(1,nparticles));
    tsi_particle = (tsi_particle<tsimax).*tsi_particle + (tsi_particle>=tsimax).*(tsimin+0.75*(tsimax-tsimin)*rand(1,nparticles));
    tsi_particle = (tsi_particle>tsimin).*tsi_particle + (tsi_particle<=tsimin).*(tsimin+0.25*(tsimax-tsimin)+0.75*(tsimax-tsimin)*rand(1,nparticles));
    twp_particle = (twp_particle<twpmax).*twp_particle + (twp_particle>=twpmax).*(twpmin+0.75*(twpmax-twpmin)*rand(1,nparticles));
    twp_particle = (twp_particle>twpmin).*twp_particle + (twp_particle<=twpmin).*(twpmin+0.25*(twpmax-twpmin)+0.75*(twpmax-twpmin)*rand(1,nparticles));
    td_particle = (td_particle<tdmax).*td_particle + (td_particle>=tdmax).*(tdmin+0.75*(tdmax-tdmin)*rand(1,nparticles));
    td_particle = (td_particle>tdmin).*td_particle + (td_particle<=tdmin).*(tdmin+0.25*(tdmax-tdmin)+0.75*(tdmax-tdmin)*rand(1,nparticles));
    tvp_particle = (tvp_particle<tvpmax).*tvp_particle + (tvp_particle>=tvpmax).*(tvpmin+0.75*(tvpmax-tvpmin)*rand(1,nparticles));
    tvp_particle = (tvp_particle>tvpmin).*tvp_particle + (tvp_particle<=tvpmin).*(tvpmin+0.25*(tvpmax-tvpmin)+0.75*(tvpmax-tvpmin)*rand(1,nparticles));
    tv1m_particle = (tv1m_particle<tv1mmax).*tv1m_particle + (tv1m_particle>=tv1mmax).*(tv1mmin+0.75*(tv1mmax-tv1mmin)*rand(1,nparticles));
    tv1m_particle = (tv1m_particle>tv1mmin).*tv1m_particle + (tv1m_particle<=tv1mmin).*(tv1mmin+0.25*(tv1mmax-tv1mmin)+0.75*(tv1mmax-tv1mmin)*rand(1,nparticles));
    tv2m_particle = (tv2m_particle<tv2mmax).*tv2m_particle + (tv2m_particle>=tv2mmax).*(tv2mmin+0.75*(tv2mmax-tv2mmin)*rand(1,nparticles));
    tv2m_particle = (tv2m_particle>tv2mmin).*tv2m_particle + (tv2m_particle<=tv2mmin).*(tv2mmin+0.25*(tv2mmax-tv2mmin)+0.75*(tv2mmax-tv2mmin)*rand(1,nparticles));
    twm_particle = (twm_particle<twmmax).*twm_particle + (twm_particle>=twmmax).*(twmmin+0.75*(twmmax-twmmin)*rand(1,nparticles));
    twm_particle = (twm_particle>twmmin).*twm_particle + (twm_particle<=twmmin).*(twmmin+0.25*(twmmax-twmmin)+0.75*(twmmax-twmmin)*rand(1,nparticles));
    to_particle = (to_particle<tomax).*to_particle + (to_particle>=tomax).*(tomin+0.75*(tomax-tomin)*rand(1,nparticles));
    to_particle = (to_particle>tomin).*to_particle + (to_particle<=tomin).*(tomin+0.25*(tomax-tomin)+0.75*(tomax-tomin)*rand(1,nparticles));
    xk_particle = (xk_particle<xkmax).*xk_particle + (xk_particle>=xkmax).*(xkmin+0.75*(xkmax-xkmin)*rand(1,nparticles));
    xk_particle = (xk_particle>xkmin).*xk_particle + (xk_particle<=xkmin).*(xkmin+0.25*(xkmax-xkmin)+0.75*(xkmax-xkmin)*rand(1,nparticles));
    ucsi_particle = (ucsi_particle<ucsimax).*ucsi_particle + (ucsi_particle>=ucsimax).*(ucsimin+0.75*(ucsimax-ucsimin)*rand(1,nparticles));
    ucsi_particle = (ucsi_particle>ucsimin).*ucsi_particle + (ucsi_particle<=ucsimin).*(ucsimin+0.25*(ucsimax-ucsimin)+0.75*(ucsimax-ucsimin)*rand(1,nparticles));
    uc_particle = (uc_particle<ucmax).*uc_particle + (uc_particle>=ucmax).*(ucmin+0.75*(ucmax-ucmin)*rand(1,nparticles));
    uc_particle = (uc_particle>ucmin).*uc_particle + (uc_particle<=ucmin).*(ucmin+0.25*(ucmax-ucmin)+0.75*(ucmax-ucmin)*rand(1,nparticles));
    uv_particle = (uv_particle<uvmax).*uv_particle + (uv_particle>=uvmax).*(uvmin+0.75*(uvmax-uvmin)*rand(1,nparticles));
    uv_particle = (uv_particle>uvmin).*uv_particle + (uv_particle<=uvmin).*(uvmin+0.25*(uvmax-uvmin)+0.75*(uvmax-uvmin)*rand(1,nparticles));
    
    if(debug)
        [tr_particle; tsi_particle; twp_particle; td_particle; tvp_particle; tv1m_particle; tv2m_particle; twm_particle; to_particle; xk_particle; ucsi_particle; uc_particle; uv_particle]
    end
    
%     figure(2)
%     subplot(3,1,1)
%     plot(a_particle,'ro'),hold on,plot(amax*ones(size(a_particle)),'k'),plot(amin*ones(size(a_particle)),'k'),hold off
%     subplot(3,1,2)
%     plot(b_particle,'ro'),hold on,plot(bmax*ones(size(b_particle)),'k'),plot(bmin*ones(size(b_particle)),'k'),hold off
%     subplot(3,1,3)
%     plot(eps_particle,'ro'),hold on,plot(epsmax*ones(size(eps_particle)),'k'),plot(epsmin*ones(size(eps_particle)),'k'),hold off
    
end

fprintf('tr = %8.4f, tsi = %8.4f, twp=%8.4f, td = %8.4f\n',tr_global, tsi_global, twp_global,td_global)
fprintf('tvp = %8.4f, tv1m = %8.4f, tv2m=%8.4f\n',tvp_global, tv1m_global, tv2m_global)
fprintf('twm = %8.4f, to = %8.4f, xk=%8.4f\n',twm_global, to_global, xk_global)
fprintf('ucsi = %8.4f, uc = %8.4f, uv=%8.4f\n',ucsi_global, uc_global, uv_global)

% copying and testing
ktime=0;
u=zeros(1,1);
v=v0*ones(1,1);
w=w0*ones(1,1);
usave1 = zeros(nsave,nparticles);
vsave1 = usave1; wsave1=usave1;

% tr_global
% tsi_global
% twp_global

for ntime=1:nsteps
    
        p=0; % heaviside functions
        q=0;
        p=(u>=uc_global);
        q=(u>=uv_global);
        dv=(1.-p).*(1.-v)./((1.-q).*tv1m_global+tv2m_global.*q) - p.*v./tvp_global;
        dw=(1.-p).*(1.-w)./twm_global - p.*w./twp_global;
        v=v+dt*dv; % solving/updating v and w
        w=w+dt*dw;
        jfi=-v.*p.*(u-uc_global).*(1.-u)./td_global;
        jso=u.*(1.-p)./to_global+p./tr_global;
        jsi=-w.*(1.+tanh(xk_global.*(u-ucsi_global)))./(2.*tsi_global);
        
        % stimulus
        istim=0;
%        if(mod(ntime,period/dt)<1/dt)
        if((mod(ntime,period/dt)<7/dt)&&(mod(ntime,period/dt)>2/dt))
            istim=stimmag;
        end
        
        u=u - (jfi+jso+jsi-istim)*dt; % updating/integrating u
        
        if(mod(ntime,round(1/dt))==0)
            ktime=ktime+1;
            usave1(ktime,:) = u;
            vsave1(ktime,:) = v;
            wsave1(ktime,:) = w;
        end

end
test_dist=sum((usave0-usave1).^2);
%fprintf('new checked dist = %8.2f\n',test_dist)
fprintf('should agree checked dist = %8.2f\n',sum((usave0-u_global).^2))

figure(2)
%global_error_check = sum((usave0-u_global).^2);
subplot(2,1,1)
plot(tsave,usave0,'k',tsave,u_global,'b',tsave,usave1,'r','linewidth',2),ylim([-0.1 1.6])
xlabel('Time'),title(['iteration = ' num2str(k)])

subplot(2,1,2)
plot(best_error,'ko')
xlabel('Iteration'),ylabel('Error')


%%
% From the local best calculated, find those whose errors meet some
% threshold and form a distribution from them for each parameter
% use dist_best for the distances (errors)
% use tr_best, tsi_best, twp_best, td_best for the corresponding param vals

figure(10)
plot(sort(dist_best),'x')
ylabel('Error')

% Define threshold for error to be in the distribution
%threshold_error = 2*min(dist_best); % 2 * the best value
%threshold_error = 0.5*max(dist_best); % 1/2 * the worst value
%threshold_error = 0.1*max(dist_best); % 0.1 * the worst value
%threshold_error = 0.05*max(dist_best);
threshold_error = 0.5;

% define the population 
indPop = find(dist_best<threshold_error);

fprintf('population size = %d\n',length(indPop));
fprintf('tr:  mean = %8.3f\n',...
    mean(tr_best(indPop)));
fprintf('to:  mean = %8.3f\n',...
    mean(to_best(indPop)));
fprintf('tsi: mean = %8.3f\n',...
    mean(tsi_best(indPop)));
fprintf('twp: mean = %8.3f\n',...
    mean(twp_best(indPop)));
fprintf('twm: mean = %8.3f\n',...
    mean(twm_best(indPop)));
fprintf('td:  mean = %8.3f\n',...
    mean(td_best(indPop)));
fprintf('tvp: mean = %8.3f\n',...
    mean(tvp_best(indPop)));
fprintf('tv1m: mean = %8.3f\n',...
    mean(tv1m_best(indPop)));
fprintf('tv2m: mean = %8.3f\n',...
    mean(tv2m_best(indPop)));
fprintf('xk:  mean = %8.3f\n',...
    mean(xk_best(indPop)));
fprintf('ucsi:  mean = %8.3f\n',...
    mean(ucsi_best(indPop)));
fprintf('uc:  mean = %8.3f\n',...
    mean(uc_best(indPop)));
fprintf('uv:  mean = %8.3f\n',...
    mean(uv_best(indPop)));

figure(11)
subplot(2,2,1)
histogram(tr_best(indPop)),hold on, 
ylimhere = get(gca,'ylim'); 
plot([tr_global tr_global],[0 max(ylimhere)],'r'),hold off
xlabel('tr')
subplot(2,2,2)
histogram(tsi_best(indPop)),xlabel('tsi'),hold on, 
ylimhere = get(gca,'ylim'); 
plot([tsi_global tsi_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,3)
histogram(twp_best(indPop)),xlabel('twp'),hold on, 
ylimhere = get(gca,'ylim'); 
plot([twp_global twp_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,4)
histogram(td_best(indPop)),xlabel('td'),hold on, 
ylimhere = get(gca,'ylim'); 
plot([td_global td_global],[0 max(ylimhere)],'r'),hold off

figure(12)
subplot(2,2,1)
histogram(tvp_best(indPop)),xlabel('tvp'),hold on, 
ylimhere = get(gca,'ylim'); 
plot([tvp_global tvp_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,2)
histogram(tv1m_best(indPop)),xlabel('tv1m'),hold on, 
ylimhere = get(gca,'ylim'); 
plot([tv1m_global tv1m_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,3)
histogram(tv2m_best(indPop)),xlabel('tv2m'),hold on, 
ylimhere = get(gca,'ylim'); 
plot([tv2m_global tv2m_global],[0 max(ylimhere)],'r'),hold off

figure(13)
subplot(2,2,1)
histogram(twm_best(indPop)),xlabel('twm'), hold on
ylimhere = get(gca,'ylim'); 
plot([twm_global twm_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,2)
histogram(to_best(indPop)),xlabel('to'), hold on
ylimhere = get(gca,'ylim'); 
plot([to_global to_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,3)
histogram(xk_best(indPop)),xlabel('xk'), hold on
ylimhere = get(gca,'ylim'); 
plot([xk_global xk_global],[0 max(ylimhere)],'r'),hold off

figure(14)
subplot(2,2,1)
histogram(ucsi_best(indPop)),xlabel('ucsi'), hold on
ylimhere = get(gca,'ylim'); 
plot([ucsi_global ucsi_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,2)
histogram(uc_best(indPop)),xlabel('uc'), hold on
ylimhere = get(gca,'ylim'); 
plot([uc_global uc_global],[0 max(ylimhere)],'r'),hold off
subplot(2,2,3)
histogram(uv_best(indPop)),xlabel('uv'), hold on
ylimhere = get(gca,'ylim'); 
plot([uv_global uv_global],[0 max(ylimhere)],'r'),hold off

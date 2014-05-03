% isabel
W=500; H=500; D=100; T=48;

% gen seeds
STEP=100;
gen_seed('seeds.tmp.txt', STEP/2:STEP:W, STEP/2:STEP:H, D/2);

% run 
cmd = 'parallelPathline_omp /data/flow/isabel_all/all.list -loader=1 -saver=2 -limit=3 -omp_nproc=7 -stepsize=.125 -seedfile=seeds.tmp.txt -out=true.out' 
system(cmd);

cmd = 'parallelPathline_omp /data/flow/isabel_all/fitted/fitted8/all_bezier_rms.list -loader=bezier -saver=2 -limit=3 -omp_nproc=7 -stepsize=0.015625 -seedfile=seeds.tmp.txt -out=sampled8.out' 
system(cmd);

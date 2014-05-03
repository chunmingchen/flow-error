% output: trace_list in cell. each cell size: 4 x traces

filename='/data/flow/isabel_all/2d/pathlines_seeds.out'

fin = fopen(filename, 'rb')
bmin = fread(fin, [4,1], 'float32')
bmax = fread(fin, [4,1], 'float32')

traces_list=[];
while(1)
    traces=fread(fin, 1, 'int');
    if traces==-1
        break;
    else
        traces_list(end+1) = traces;
    end
        
end

trace_list = {};
for traces = traces_list
    trace_list{end+1} = fread(fin, [4, traces], 'float32');
end

fclose(fin);
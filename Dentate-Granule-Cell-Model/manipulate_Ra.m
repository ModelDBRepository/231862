function neuron = manipulate_Ra(neuron,boolean,regions)
% changes Ra in a tree region by a factor of 1e7 in order to exclude it from any
% current exchange (boolean = 1 ). Boolean = 0 reverses the change 

if nargin < 3
    regions = {'all'};
end
if ~iscell(regions)
    regions = {regions};
end
if nargin < 2
    boolean = 1;
end

if boolean
    scale = 1000000;
else
    scale = 1/1000000;
end

for t = 1:numel(neuron.mech)
    fields = fieldnames(neuron.mech{t});
    if ~any(strcmp(regions,'all'))
        fields = intersect(fields,regions);  % only take the specified regions
    end
    for f1 = 1:numel(fields)
        if ~isfield(neuron.mech{t}.(fields{f1}),'pas')
            if isfield(neuron.mech{t}.all,'pas')
                neuron.mech{t}.(fields{f1}).pas = neuron.mech{t}.all.pas;
            else
                error('cannot manipulate Ra because mechanism "pas" does neither exist in this region nor defined for all regions')
            end
        end
        neuron.mech{t}.(fields{f1}).pas.Ra = neuron.mech{t}.(fields{f1}).pas.Ra * scale;
    end
end

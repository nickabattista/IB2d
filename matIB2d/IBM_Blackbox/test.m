%%
% Get input file without braces.
% Outputs grid info:
params = read_input_file('input2d.IB2D');
Nx = params{find(strcmp({params{:,1}},'Nx')),2}
Ny = params{find(strcmp({params{:,1}},'Ny')),2}
Lx = params{find(strcmp({params{:,1}},'Lx')),2}
Ly = params{find(strcmp({params{:,1}},'Ly')),2}
supp = params{find(strcmp({params{:,1}},'supp')),2}
d_springs = params{find(strcmp({params{:,1}},'damped_springs')),2}
str_name = params{find(strcmp({params{:,1}},'string_name')),2}


params

%%
% Get input file with braces.
% Outputs grid info
params = read_input_file('input2d.IB2D');
grid_Info = params{find(strcmp({params{:,1}},'GridParameters')),2}

grid_Info

%modelStruct = mStruct;
eval( sprintf( '%s = mStruct;', structVar ) );

jacFileName = [ parserPath '/outputC/model/aJac.c' ];

disp( 'Computing RHS derivatives ...' );
x = sym('x_%d', [1, nStates]);
if isfield( mStruct, 'u' )
    u = sym('u_%d', [1, numel( fieldnames( mStruct.u ))]);
else
    u = sym([]);
end
S = sym('s_%d_%d', [nStates + nPars, nStates]).';
p = sym('p_%d', [nPars,1]).';

m_file = textread( odeInput, '%s', 'delimiter', '\n' );

m_file = sprintf( '%s\n', m_file{2:end} );
eval( m_file );

dfdy = jacobian(dx, x).';
dfdp = jacobian(dx, p);

for a = 1 : nStates
    dfdp( a, a+nPars ) = 0;
end

if ( aJac == 1 )
    
    %% Generate the analytical sensitivity equations
    disp( 'Generating sensitivity equations ...' );
    ccode( ( dfdy.' * S + dfdp ).', 'file', jacFileName );

    disp( 'Reading sensitivity equation C-file' );
    c_file = textread( jacFileName, '%s', 'delimiter', '\n' );

    disp( 'Replacing state names' );
    c_file = regexprep( c_file, 'MatrixWithNoName\[(\d*)\]', 'NV_DATA_S(ySdot[$1])' );
    c_file = regexprep( c_file, 'A0\[(\d*)\]', 'NV_DATA_S(ySdot[$1])' );
    c_file = regexprep( c_file, 'p_(\d*)', 'data->p[$1-1]' );
    c_file = regexprep( c_file, 'u_(\d*)', 'data->u[$1-1]' );
    c_file = regexprep( c_file, 'x_(\d*)', 'stateVars[$1-1]' );
    c_file = regexprep( c_file, 's_(\d*)_(\d*)', 'NV_DATA_S(yS[$1-1])[$2-1]' );
    c_file = regexprep( c_file, '^t0 = ', '' );

    fullstring = sprintf( '%s\n', c_file{:} );

    [startIndex, endIndex, tokIndex, matchStr, tokenStr] = regexp( fullstring, 't(\d*)' );

    maxT = 0;
    for a = 1 : length( tokenStr )
        val = str2num( tokenStr{a}{1} );
        maxT = max( [ val, maxT ] );
    end

    preString = sprintf('int sensRhs (int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {\n\n');
    for a = 1 : maxT
        preString = sprintf( '%srealtype t%d;', preString, a );
    end
    postString = sprintf( '\nstruct mData *data = ( struct mData * ) user_data;\nrealtype *stateVars;\nstateVars = NV_DATA_S(y);\n\n' );
    sensRHS = sprintf( '%s\n%s\n%s\nreturn 0;\n};\n', preString, postString, fullstring );

else
    sensRHS = '';
end

if ( fJac == 1 )
    %% Analytical Jacobian
    disp( 'Generating jacobian ...' );
    %ccode( (dfdy.').', 'file', jacFileName );
    ccode( (dfdy.').', 'file', jacFileName );
    %ccode( (sym(0).').', 'file', jacFileName );

    disp( 'Reading jacobian C-file' );
    c_file = textread( jacFileName, '%s', 'delimiter', '\n' );

    disp( 'Replacing state names' );
    c_file = regexprep( c_file, 'MatrixWithNoName\[(\d*)\]', 'DENSE_COL(Jac,$1)' );
    c_file = regexprep( c_file, 'A0\[(\d*)\]', 'DENSE_COL(Jac,$1)' );
    c_file = regexprep( c_file, 'p_(\d*)', 'data->p[$1-1]' );
    c_file = regexprep( c_file, 'u_(\d*)', 'data->u[$1-1]' );
    c_file = regexprep( c_file, 'x_(\d*)', 'stateVars[$1-1]' );
    c_file = regexprep( c_file, '^t0 = ', '' );

    fullstring = sprintf( '%s\n', c_file{:} );

    [startIndex, endIndex, tokIndex, matchStr, tokenStr] = regexp( fullstring, 't(\d*)' );

    maxT = 0;
    for a = 1 : length( tokenStr )
        val = str2num( tokenStr{a}{1} );
        maxT = max( [ val, maxT ] );
    end

    preString = sprintf('int fJac (long int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {\n\n');
    for a = 1 : maxT
        preString = sprintf( '%srealtype t%d;', preString, a );
    end
    postString = sprintf( '\nstruct mData *data = ( struct mData * ) user_data;\nrealtype *stateVars;\nstateVars = NV_DATA_S(y);\n\n' );
    sensRHS = sprintf( '%s%s\n%s\n%s\nreturn 0;\n};\n', sensRHS, preString, postString, fullstring );
end

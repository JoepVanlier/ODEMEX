#include "../dxdt.h"



int rhs( realtype t, N_Vector y, N_Vector ydot, void *f_data ) {

	struct mData *data = ( struct mData * ) f_data;

	realtype a1 , a2 , alf , order , rho ;
	int k1 ;

	realtype *stateVars;
	realtype *ydots;

	stateVars = NV_DATA_S(y);
	ydots = NV_DATA_S(ydot);

	
	/*Parameters:*/
	a1  =data->p[0];
	a2  =data->p[1];
	alf =data->p[2];
	
	rho     =2;
	order   =24;
	
	ydots[0]      = (a1 / (1 +a2 *intPow(stateVars[23],rho) ) ) -alf *stateVars[0];
	
	for(k1=1;k1<=order;k1++){
	ydots[(1 +k1 ) -1] =data->p[(4 +k1 ) -1] *stateVars[(1 +k1 -1 ) -1] -alf *stateVars[(1 +k1 ) -1];
	}


	#ifdef NON_NEGATIVE
		return 0;
	#else
		return 0;
	#endif

};
int fJac (long int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {



struct mData *data = ( struct mData * ) user_data;
realtype *stateVars;
stateVars = NV_DATA_S(y);


DENSE_COL(Jac,0)[0] = -data->p[2];
DENSE_COL(Jac,0)[1] = data->p[4];
DENSE_COL(Jac,1)[1] = -data->p[2];
DENSE_COL(Jac,1)[2] = data->p[5];
DENSE_COL(Jac,2)[2] = -data->p[2];
DENSE_COL(Jac,2)[3] = data->p[6];
DENSE_COL(Jac,3)[3] = -data->p[2];
DENSE_COL(Jac,3)[4] = data->p[7];
DENSE_COL(Jac,4)[4] = -data->p[2];
DENSE_COL(Jac,4)[5] = data->p[8];
DENSE_COL(Jac,5)[5] = -data->p[2];
DENSE_COL(Jac,5)[6] = data->p[9];
DENSE_COL(Jac,6)[6] = -data->p[2];
DENSE_COL(Jac,6)[7] = data->p[10];
DENSE_COL(Jac,7)[7] = -data->p[2];
DENSE_COL(Jac,7)[8] = data->p[11];
DENSE_COL(Jac,8)[8] = -data->p[2];
DENSE_COL(Jac,8)[9] = data->p[12];
DENSE_COL(Jac,9)[9] = -data->p[2];
DENSE_COL(Jac,9)[10] = data->p[13];
DENSE_COL(Jac,10)[10] = -data->p[2];
DENSE_COL(Jac,10)[11] = data->p[14];
DENSE_COL(Jac,11)[11] = -data->p[2];
DENSE_COL(Jac,11)[12] = data->p[15];
DENSE_COL(Jac,12)[12] = -data->p[2];
DENSE_COL(Jac,12)[13] = data->p[16];
DENSE_COL(Jac,13)[13] = -data->p[2];
DENSE_COL(Jac,13)[14] = data->p[17];
DENSE_COL(Jac,14)[14] = -data->p[2];
DENSE_COL(Jac,14)[15] = data->p[18];
DENSE_COL(Jac,15)[15] = -data->p[2];
DENSE_COL(Jac,15)[16] = data->p[19];
DENSE_COL(Jac,16)[16] = -data->p[2];
DENSE_COL(Jac,16)[17] = data->p[20];
DENSE_COL(Jac,17)[17] = -data->p[2];
DENSE_COL(Jac,17)[18] = data->p[21];
DENSE_COL(Jac,18)[18] = -data->p[2];
DENSE_COL(Jac,18)[19] = data->p[22];
DENSE_COL(Jac,19)[19] = -data->p[2];
DENSE_COL(Jac,19)[20] = data->p[23];
DENSE_COL(Jac,20)[20] = -data->p[2];
DENSE_COL(Jac,20)[21] = data->p[24];
DENSE_COL(Jac,21)[21] = -data->p[2];
DENSE_COL(Jac,21)[22] = data->p[25];
DENSE_COL(Jac,22)[22] = -data->p[2];
DENSE_COL(Jac,22)[23] = data->p[26];
DENSE_COL(Jac,23)[0] = data->p[0]*data->p[1]*stateVars[23]*1.0/pow(data->p[1]*(stateVars[23]*stateVars[23])+1.0,2.0)*-2.0;
DENSE_COL(Jac,23)[23] = -data->p[2];
DENSE_COL(Jac,23)[24] = data->p[27];
DENSE_COL(Jac,24)[24] = -data->p[2];

return 0;
};


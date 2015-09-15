#include "../dxdt.h"

realtype addTwo( realtype in1 , realtype  in2 ) {

	realtype y ;
	;
	
	y =in1 +in2;
	/**/
	/*Joep Vanlier, 2011*/
	/**/
	/*Licensing:*/
	/*Copyright (C) 2009-2011 Joep Vanlier. All rights*/
	/*reserved.*/
	/**/
	/*Contact:joep.vanlier@gmail.com*/
	/**/
	/*This file is part of the puaMAT.*/
	/**/
	/*puaMAT is free software: you can redistribute it*/
	/*and/or modify it under the terms of the GNU General*/
	/*Public License as published by the Free Software*/
	/*Foundation, either version 3 of the License, or (at*/
	/*your option) any later version.*/
	/**/
	/*puaMAT is distributed in the hope that it will be*/
	/*useful, but WITHOUT ANY WARRANTY; without even the*/
	/*implied warranty of MERCHANTABILITY or FITNESS FOR A*/
	/*PARTICULAR PURPOSE.  See the GNU General Public*/
	/*License for more details.*/
	/**/
	/*You should have received a copy of the GNU General*/
	/*Public License along with puaMAT.  If not, see*/
	/*http://www.gnu.org/licenses/*/
	/**/
	return y;

}



int rhs( realtype t, N_Vector y, N_Vector ydot, void *f_data ) {

	struct mData *data = ( struct mData * ) f_data;

	realtype c1 , dx , p1 , p2 , p3 , temp1 , temp2 , u1 , x1 , x2 ;

	realtype *stateVars;
	realtype *ydots;

	stateVars = NV_DATA_S(y);
	ydots = NV_DATA_S(ydot);

	
	p1 =data->p[0];
	p2 =data->p[1];
	p3 =data->p[2];
	x1 =stateVars[0];
	x2 =stateVars[1];
	u1 =data->u[0];
	c1 =1;
	
	temp1 =u1 -p1 *x1;
	temp2 =p2 *x2;
	
	ydots[0] =temp1 +temp2;
	ydots[1] =p1 *x1 -addTwo((p3/c1),p2) *x2;
	
	
	/**/
	/*Joep Vanlier, 2011*/
	/**/
	/*Licensing:*/
	/*Copyright (C) 2009-2011 Joep Vanlier. All rights*/
	/*reserved.*/
	/**/
	/*Contact:joep.vanlier@gmail.com*/
	/**/
	/*This file is part of the puaMAT.*/
	/**/
	/*puaMAT is free software: you can redistribute it*/
	/*and/or modify it under the terms of the GNU General*/
	/*Public License as published by the Free Software*/
	/*Foundation, either version 3 of the License, or (at*/
	/*your option) any later version.*/
	/**/
	/*puaMAT is distributed in the hope that it will be*/
	/*useful, but WITHOUT ANY WARRANTY; without even the*/
	/*implied warranty of MERCHANTABILITY or FITNESS FOR A*/
	/*PARTICULAR PURPOSE.  See the GNU General Public*/
	/*License for more details.*/
	/**/
	/*You should have received a copy of the GNU General*/
	/*Public License along with puaMAT.  If not, see*/
	/*http://www.gnu.org/licenses/*/
	/**/


	#ifdef NON_NEGATIVE
		return 0;
	#else
		return 0;
	#endif

};
int sensRhs (int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {

realtype t1;realtype t2;

struct mData *data = ( struct mData * ) user_data;
realtype *stateVars;
stateVars = NV_DATA_S(y);


t2 = data->p[1]+data->p[2];
NV_DATA_S(ySdot[0])[0] = -stateVars[0]-data->p[0]*NV_DATA_S(yS[0])[0]+data->p[1]*NV_DATA_S(yS[0])[1];
NV_DATA_S(ySdot[0])[1] = stateVars[0]+data->p[0]*NV_DATA_S(yS[0])[0]-NV_DATA_S(yS[0])[1]*t2;
NV_DATA_S(ySdot[1])[0] = stateVars[1]-data->p[0]*NV_DATA_S(yS[1])[0]+data->p[1]*NV_DATA_S(yS[1])[1];
NV_DATA_S(ySdot[1])[1] = -stateVars[1]+data->p[0]*NV_DATA_S(yS[1])[0]-NV_DATA_S(yS[1])[1]*t2;
NV_DATA_S(ySdot[2])[0] = -data->p[0]*NV_DATA_S(yS[2])[0]+data->p[1]*NV_DATA_S(yS[2])[1];
NV_DATA_S(ySdot[2])[1] = -stateVars[1]+data->p[0]*NV_DATA_S(yS[2])[0]-NV_DATA_S(yS[2])[1]*t2;
NV_DATA_S(ySdot[3])[0] = -data->p[0]*NV_DATA_S(yS[3])[0]+data->p[1]*NV_DATA_S(yS[3])[1];
NV_DATA_S(ySdot[3])[1] = data->p[0]*NV_DATA_S(yS[3])[0]-NV_DATA_S(yS[3])[1]*t2;
NV_DATA_S(ySdot[4])[0] = -data->p[0]*NV_DATA_S(yS[4])[0]+data->p[1]*NV_DATA_S(yS[4])[1];
NV_DATA_S(ySdot[4])[1] = data->p[0]*NV_DATA_S(yS[4])[0]-NV_DATA_S(yS[4])[1]*t2;

return 0;
};


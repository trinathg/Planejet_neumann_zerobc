#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

//////////////////////Boundary conditions on the top, bottom and side walls 

void bcs(vertex * node, fval * fvar)
{	
	int ind; 
	
	tdma1y(fvar,0); 

	tdma2x(fvar,0);  //d2u/dx2  //Computing the second derivatives needed for pressure BCs
	tdma2y(fvar,0);  //d2u/dy2  //Computing the second derivatives needed for pressure BCs
	
	tdma2x(fvar,1);  //d2v/dx2  //Computing the second derivatives needed for pressure BCs
	tdma2y(fvar,1);  //d2v/dy2 

	//Top boundary 

	for(int i=1;i<nx;i++)
	{
		ind = i + ny*str_x; 					
		
		fvar[ind].u[0] = - B2*fvar[ind-str_x].u[0] - B3*fvar[ind-2*str_x].u[0] - B4*fvar[ind-3*str_x].u[0]; 
		fvar[ind].u[1] = 0.0; 	
		fvar[ind].u[2] = -(1./B1)*dy*( (1./Re)*fvar[ind].uyy[1] )- B2*fvar[ind-str_x].u[2] - B3*fvar[ind-2*str_x].u[2] - B4*fvar[ind-3*str_x].u[2];
	}

	//Right boundary 

	for(int j=0;j<=ny;j++)
	{
		ind = nx + j*str_x; 			

		fvar[ind].u[0] = - B2*fvar[ind-1].u[0] - B3*fvar[ind-2].u[0] - B4*fvar[ind-3].u[0]; 
		fvar[ind].u[1] = - B2*fvar[ind-1].u[1] - B3*fvar[ind-2].u[1] - B4*fvar[ind-3].u[1];  
		//fvar[ind].u[2] = -(1./B1)*dx*( - fvar[ind].u[1]*fvar[ind].uy[0] + (1./Re)*( fvar[ind].uxx[0] + fvar[ind].uyy[0] ) )- B2*fvar[ind-1].u[2] - B3*fvar[ind-2].u[2] - B4*fvar[ind-3].u[2]; 

                fvar[ind].u[2] = - B2*fvar[ind-1].u[2] - B3*fvar[ind-2].u[2] - B4*fvar[ind-3].u[2]; 

	}

	//Bottom boundary

	for(int i=1;i<=nx;i++)
	{			
		fvar[i].u[0] = - (B2*fvar[i+str_x].u[0] + B3*fvar[i+2*str_x].u[0] + B4*fvar[i+3*str_x].u[0]);
		fvar[i].u[1] = 0.0;
		fvar[i].u[2] = (1./B1)*dy*( (1./Re)*fvar[i].uyy[1]) - ( B2*fvar[i+str_x].u[2] + B3*fvar[i+2*str_x].u[2] + B4*fvar[i+3*str_x].u[2] );
	}

	//Left boundary

	for(int j=0;j<=ny;j++)
	{
	 	ind = j*str_x;  	 		
	
		fvar[ind].u[0] = 1.0/( 1.+sinh(node[ind].x[1]-l_y/2.)*sinh(node[ind].x[1]-l_y/2.) );;	
		fvar[ind].u[1] = 0.0;	
		fvar[ind].u[2] = (1./B1)*dx*(1./Re)*( fvar[ind].uxx[0] + fvar[ind].uyy[0] ) - ( B2*fvar[ind+1].u[2] + B3*fvar[ind+2].u[2] + B4*fvar[ind+3].u[2]) ;
	}
}

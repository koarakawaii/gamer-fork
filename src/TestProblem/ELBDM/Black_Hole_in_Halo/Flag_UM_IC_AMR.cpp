#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_User_Template
// Description :  Template of user-defined flag criteria
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_UM_IC_AMR( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

// put your favorite flag criteria here
// ##########################################################################################################

// Define the AMR refinement similiar to the reconstructed halo UM_IC
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   // flag cells within the target region [Threshold ... BoxSize-Threshold]
   const double EdgeL = Threshold[0];
   const double EdgeR = amr->BoxSize[0]-Threshold[0];    // here we have assumed a cubic box

   bool Flag;

   if (  Pos[0] >= (EdgeL+0.5*dh)  &&  Pos[0] < (EdgeR-0.5*dh)  &&
         Pos[1] >= (EdgeL+0.5*dh)  &&  Pos[1] < (EdgeR-0.5*dh)  &&
         Pos[2] >= (EdgeL+0.5*dh)  &&  Pos[2] < (EdgeR-0.5*dh)     )
      Flag = true;

   else
      Flag = false;
   
// ##########################################################################################################


   return Flag;

} // FUNCTION : Flag_User_Template

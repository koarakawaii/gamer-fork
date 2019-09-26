#include "GAMER.h"

void Init_ByFunction_AssignData( const int lv, const bool SetFlu, const bool SetMag );




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFunction
// Description :  Set up the initial condition by invoking Init_ByFunction_AssignData()
//
// Note        :  1. Invoke the alternative function LB_Init_ByFunction() when LOAD_BALANCE is adopted
//-------------------------------------------------------------------------------------------------------
void Init_ByFunction()
{

// invoke the alternative load-balance function
#  ifdef LOAD_BALANCE
   LB_Init_ByFunction();
   return;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// construct levels 0 ~ NLEVEL-1
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing level %d ... ", lv );

      if ( lv == 0 )
      Init_BaseLevel();

      const bool SetFlu_Yes = true;
      const bool SetMag_Yes = true;
      Init_ByFunction_AssignData( lv, SetFlu_Yes, SetMag_Yes );

      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_NO );

      if ( lv != TOP_LEVEL )
      {
         Flag_Real( lv, USELB_NO );

         MPI_ExchangeBoundaryFlag( lv );

         Flag_Buffer( lv );

         Init_Refine( lv );
      }

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

   } // for (int lv=0; lv<NLEVEL; lv++)


// restrict all variables to be consistent with the finite volume scheme
   if ( OPT__INIT_RESTRICT )
   {
      for (int lv=TOP_LEVEL-1; lv>=0; lv--)
      {
         Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv], NULL_INT, NULL_INT, _TOTAL, _MAG );

         Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_NO );
      } // for (int lv=NLEVEL-2; lv>=0; lv--)
   } // if ( OPT__INIT_RESTRICT )


// get the total number of real patches at all ranks
   for (int lv=0; lv<NLEVEL; lv++)     Mis_GetTotalPatchNumber( lv );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_ByFunction



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFunction_AssignData
// Description :  Construct the initial condition in different models
//
// Note        :  1. Work for the option "OPT__INIT == INIT_BY_FUNCTION"
//                   --> But it is also used for "OPT__INIT == INIT_BY_FILE" if
//                       "OPT__INIT_MAG == INIT_MAG_BY_FUNCITON" or vice versa
//                2. Invoked by Init_ByFunction(), LB_Init_ByFunction(), and Init_ByFile()
//
// Parameter   :  lv     : Target refinement level
//                SetFlu : Set the fluid field
//                SetMag : Set the magnetic field
//-------------------------------------------------------------------------------------------------------
void Init_ByFunction_AssignData( const int lv, const bool SetFlu, const bool SetMag )
{

#  if   ( MODEL == HYDRO )
   Hydro_Init_ByFunction_AssignData( lv, SetFlu, SetMag );

#  elif ( MODEL == ELBDM )
   ELBDM_Init_ByFunction_AssignData( lv, SetFlu, SetMag );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

} // FUNCTION : Init_ByFunction_AssignData

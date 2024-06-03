#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_ScatterParticleData
// Description :  Scatter particle data from the root rank to all ranks
//
// Note        :  1. Typically called by a particle initial condition routine like Par_Init_ByFunction_*()
//                2. Do not put this file in "Particle/LoadBalance" since we still need to copy data from
//                   Data_Send_Flt[] to Data_Recv_Flt[] in the serial mode
//
// Parameter   :  NPar_ThisRank : Number of particles to be received by this MPI rank
//                NPar_AllRank  : Total number of particles in all MPI ranks
//                FltAttBitIdx  : Bitwise indices of the target particle attributes (e.g., _PAR_MASS | _PAR_VELX)
//                                --> A user-defined attribute with an integer index FltAttIntIdx returned by
//                                    AddParticleAttributeFlt() can be converted to a bitwise index by BIDX(FltAttIntIdx)
//                Data_Send_Flt : Pointer array for all particle attributes to be sent
//                                --> Dimension = [PAR_NATT_FLT_TOTAL][NPar_AllRank]
//                                --> Target particle attributes are set by "FltAttBitIdx"
//                Data_Recv_Flt : Pointer array for all particle attributes to be received
//
// Return      :  Data_Recv_Flt
//-------------------------------------------------------------------------------------------------------
void Par_ScatterParticleData( const long NPar_ThisRank, const long NPar_AllRank, const long FltAttBitIdx,
                              real_par *Data_Send_Flt[PAR_NATT_FLT_TOTAL], real_par *Data_Recv_Flt[PAR_NATT_FLT_TOTAL] )
{

// check integer overflow in MPI
   if ( NPar_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "Total number of particles (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                 NPar_AllRank, (long)__INT_MAX__ );


// get the number of particles in each rank and set the corresponding offsets
   int NSend_Flt[MPI_NRank], SendDisp_Flt[MPI_NRank];
   int NPar_ThisRank_int = NPar_ThisRank;    // (i) convert to "int" and (ii) remove the "const" declaration
                                             // --> (ii) is necessary for OpenMPI version < 1.7

   MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend_Flt, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      SendDisp_Flt[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp_Flt[r] = SendDisp_Flt[r-1] + NSend_Flt[r-1];
   }


// send particle attributes (one at a time) from the root rank to all ranks
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
   {
      if ( FltAttBitIdx & (1L<<v) )    MPI_Scatterv( Data_Send_Flt[v], NSend_Flt, SendDisp_Flt, MPI_GAMER_REAL_PAR,
                                                     Data_Recv_Flt[v], NPar_ThisRank,           MPI_GAMER_REAL_PAR,
                                                     0, MPI_COMM_WORLD );
   }

} // FUNCTION : Par_ScatterParticleData



#endif // #ifdef PARTICLE

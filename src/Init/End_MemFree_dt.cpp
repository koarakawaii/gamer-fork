#include "GAMER.h"

#ifndef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_dt
// Description :  Free memory previously allocated by Init_MemAllocate_dt()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_dt()
{

   for (int t=0; t<2; t++)
   {
      delete [] h_dt_Array_T    [t];   h_dt_Array_T    [t] = NULL;
      delete [] h_Flu_Array_T   [t];   h_Flu_Array_T   [t] = NULL;
      delete [] h_Corner_Array_T[t];   h_Corner_Array_T[t] = NULL;
#     ifdef GRAVITY
      delete [] h_Pot_Array_T   [t];   h_Pot_Array_T   [t] = NULL;
#     endif
   } // for (int t=0; t<2; t++)

} // FUNCTION : End_MemFree_dt



#endif // #ifndef GPU

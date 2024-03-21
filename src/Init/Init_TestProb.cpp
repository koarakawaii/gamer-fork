#include "GAMER.h"


// ******************************************************************
// add the new test problem function prototypes here
// ******************************************************************
void Init_TestProb_Hydro_BlastWave();
void Init_TestProb_Hydro_AcousticWave();
void Init_TestProb_Hydro_Bondi();
void Init_TestProb_Hydro_ClusterMerger();
void Init_TestProb_Hydro_AGORA_IsolatedGalaxy();
void Init_TestProb_Hydro_Caustic();
void Init_TestProb_Hydro_SphericalCollapse();
void Init_TestProb_Hydro_KelvinHelmholtzInstability();
void Init_TestProb_Hydro_Riemann();
void Init_TestProb_Hydro_Jet();
void Init_TestProb_Hydro_Plummer();
void Init_TestProb_Hydro_Gravity();
void Init_TestProb_Hydro_MHD_ABC();
void Init_TestProb_Hydro_MHD_OrszagTangVortex();
void Init_TestProb_Hydro_MHD_LinearWave();
void Init_TestProb_Hydro_JeansInstability();
void Init_TestProb_Hydro_ParEqmIC();
void Init_TestProb_Hydro_BarredPot();
void Init_TestProb_Hydro_ParticleTest();
void Init_TestProb_Hydro_CDM_LSS();
void Init_TestProb_Hydro_Zeldovich();
void Init_TestProb_Hydro_EnergyPowerSpectrum();

void Init_TestProb_ELBDM_ExtPot();
void Init_TestProb_ELBDM_JeansInstabilityComoving();
void Init_TestProb_ELBDM_JeansInstabilityPhysical();
void Init_TestProb_ELBDM_Soliton();
void Init_TestProb_ELBDM_SelfSimilarHalo();
void Init_TestProb_ELBDM_VortexPairRotating();
void Init_TestProb_ELBDM_VortexPairLinear();
void Init_TestProb_ELBDM_IsolatedHalo();
void Init_TestProb_ELBDM_GaussianWavePacket();
void Init_TestProb_ELBDM_LSS();
void Init_TestProb_ELBDM_PlaneWave();
void Init_TestProb_ELBDM_Perturbation();
void Init_TestProb_ELBDM_HaloMerger();

void Init_TestProb_ELBDM_Halo_Stability_Test();
void Init_TestProb_ELBDM_Black_Hole_in_Halo();
void Init_TestProb_ELBDM_Black_Hole_in_Soliton();
void Init_TestProb_ELBDM_Correlation_Function();
void Init_TestProb_ELBDM_Soliton_Toy_Model();
void Init_TestProb_ELBDM_Halo_Stability_Test_No_Soliton();
void Init_TestProb_ELBDM_Halo_Stability_Test_Soliton_Substituted();



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize the target test problem
//
// Note        :  1. Use TESTPROB_ID to choose the target test problem
//                2. All test problem IDs are defined in "include/Typedef.h"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

// ******************************************************************
// add the new test problem IDs here
// ******************************************************************
   switch ( TESTPROB_ID )
   {
      case TESTPROB_NONE :                                                                                         break;

      case TESTPROB_HYDRO_BLAST_WAVE :                      Init_TestProb_Hydro_BlastWave();                       break;
      case TESTPROB_HYDRO_ACOUSTIC_WAVE :                   Init_TestProb_Hydro_AcousticWave();                    break;
      case TESTPROB_HYDRO_BONDI :                           Init_TestProb_Hydro_Bondi();                           break;
      case TESTPROB_HYDRO_CLUSTER_MERGER :                  Init_TestProb_Hydro_ClusterMerger();                   break;
      case TESTPROB_HYDRO_AGORA_ISOLATED_GALAXY :           Init_TestProb_Hydro_AGORA_IsolatedGalaxy();            break;
      case TESTPROB_HYDRO_CAUSTIC :                         Init_TestProb_Hydro_Caustic();                         break;
      case TESTPROB_HYDRO_SPHERICAL_COLLAPSE :              Init_TestProb_Hydro_SphericalCollapse();               break;
      case TESTPROB_HYDRO_KELVIN_HELMHOLTZ_INSTABILITY :    Init_TestProb_Hydro_KelvinHelmholtzInstability();      break;
      case TESTPROB_HYDRO_RIEMANN :                         Init_TestProb_Hydro_Riemann();                         break;
      case TESTPROB_HYDRO_JET :                             Init_TestProb_Hydro_Jet();                             break;
      case TESTPROB_HYDRO_PLUMMER :                         Init_TestProb_Hydro_Plummer();                         break;
      case TESTPROB_HYDRO_BARRED_POT :                      Init_TestProb_Hydro_BarredPot();                       break;
      case TESTPROB_HYDRO_GRAVITY :                         Init_TestProb_Hydro_Gravity();                         break;
      case TESTPROB_HYDRO_MHD_ABC :                         Init_TestProb_Hydro_MHD_ABC();                         break;
      case TESTPROB_HYDRO_MHD_ORSZAG_TANG_VORTEX :          Init_TestProb_Hydro_MHD_OrszagTangVortex();            break;
      case TESTPROB_HYDRO_MHD_LINEAR_WAVE :                 Init_TestProb_Hydro_MHD_LinearWave();                  break;
      case TESTPROB_HYDRO_JEANS_INSTABILITY :               Init_TestProb_Hydro_JeansInstability();                break;
      case TESTPROB_HYDRO_PARTICLE_EQUILIBRIUM_IC :         Init_TestProb_Hydro_ParEqmIC();                        break;
      case TESTPROB_HYDRO_PARTICLE_TEST :                   Init_TestProb_Hydro_ParticleTest();                    break;
      case TESTPROB_HYDRO_CDM_LSS :                         Init_TestProb_Hydro_CDM_LSS();                         break;
      case TESTPROB_HYDRO_ZELDOVICH :                       Init_TestProb_Hydro_Zeldovich();                       break;
      case TESTPROB_HYDRO_ENERGY_POWER_SPECTRUM :           Init_TestProb_Hydro_EnergyPowerSpectrum();             break;

      case TESTPROB_ELBDM_EXTPOT :                          Init_TestProb_ELBDM_ExtPot();                          break;
      case TESTPROB_ELBDM_JEANS_INSTABILITY_COMOVING :      Init_TestProb_ELBDM_JeansInstabilityComoving();        break;
//    case TESTPROB_ELBDM_JEANS_INSTABILITY_PHYSICAL :      Init_TestProb_ELBDM_JeansInstabilityPhysical();        break;
      case TESTPROB_ELBDM_SOLITON :                         Init_TestProb_ELBDM_Soliton();                         break;
      case TESTPROB_ELBDM_SELF_SIMILAR_HALO :               Init_TestProb_ELBDM_SelfSimilarHalo();                 break;
      case TESTPROB_ELBDM_VORTEX_PAIR_ROTATING :            Init_TestProb_ELBDM_VortexPairRotating();              break;
      case TESTPROB_ELBDM_VORTEX_PAIR_LINEAR :              Init_TestProb_ELBDM_VortexPairLinear();                break;
      case TESTPROB_ELBDM_ISOLATED_HALO :                   Init_TestProb_ELBDM_IsolatedHalo();                    break;
      case TESTPROB_ELBDM_GAUSSIAN_WAVE_PACKET :            Init_TestProb_ELBDM_GaussianWavePacket();              break;
      case TESTPROB_ELBDM_LSS :                             Init_TestProb_ELBDM_LSS();                             break;
      case TESTPROB_ELBDM_PLANE_WAVE :                      Init_TestProb_ELBDM_PlaneWave();                       break;
      case TESTPROB_ELBDM_PERTURBATION :                    Init_TestProb_ELBDM_Perturbation();                    break;
      case TESTPROB_ELBDM_HALO_MERGER :                     Init_TestProb_ELBDM_HaloMerger();                      break;
      case TESTPROB_ELBDM_HALO_STABILITY_TEST :             Init_TestProb_ELBDM_Halo_Stability_Test();             break;
      case TESTPROB_ELBDM_BLACK_HOLE_IN_HALO :              Init_TestProb_ELBDM_Black_Hole_in_Halo();              break;
      case TESTPROB_ELBDM_BLACK_HOLE_IN_SOLITON :           Init_TestProb_ELBDM_Black_Hole_in_Soliton();           break;
      case TESTPROB_ELBDM_CORRELATION_FUNCTION :            Init_TestProb_ELBDM_Correlation_Function();            break;
      case TESTPROB_ELBDM_SOLITON_TOY_MODEL :               Init_TestProb_ELBDM_Soliton_Toy_Model();               break;
      case TESTPROB_ELBDM_HALO_STABILITY_TEST_NO_SOLITON :  Init_TestProb_ELBDM_Halo_Stability_Test_No_Soliton();  break;
      case TESTPROB_ELBDM_HALO_STABILITY_TEST_SOLITON_SUBSTITUTED :  Init_TestProb_ELBDM_Halo_Stability_Test_Soliton_Substituted();  break;

      default: Aux_Error( ERROR_INFO, "unsupported TESTPROB_ID (%d) !!\n", TESTPROB_ID );
   } // switch( TESTPROB_ID )

} // FUNCTION : Init_TestProb

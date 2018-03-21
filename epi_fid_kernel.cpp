#include "Imaging3T/epi_base/epi_fid_zjw/epi_fid_kernel.h"

#include <cmath>

#include "Common/DHL/dhl.h"
#include "Common/MRGeometry/coordinate_transformation.h"
#include "Common/UProtocol/mrexam_uprotocol_interface.h"
#include "Common/UProtocol/mrexam_parameter_interface.h"

#include "SeqCtrl/seq_ctrl_define.h"
#include "SeqCtrl/SeqFW/SeqRTController/seq_exception.h"
#include "Common/MRGeometry/rotation_matrix.h"
#include "SeqCtrl/SeqFW/sequence_context.h"
#include "SeqCtrl/SeqFW/SeqSysAccessor/system_information_interface.h"
#include "SeqCtrl/SeqFW/system_info.h"


#include "inl/Common/seq_logger.h"
#include "inl/Common/SeqAuxil/grad_performance.h"
#include "inl/Common/SeqAuxil/slice_scan_ordering.h"
#include "Common/UProtocol/MRUProtAccessor/uprotocol_parameter_extraction.h"
#include "inl/Common/SeqModule/fill_time_stb.h"
#include "inl/Common/SeqModule/acquisition_index.h"
#include "inl/Common/SeqModule/fill_dhl.h"
#include "inl/Common/SeqModule/trigger_stb.h"

#include "Imaging3T/epi_base/Common/epi_excit_stb.h"
#include "Imaging3T/epi_base/Common/epi_prephase_stb.h"
#include "Imaging3T/epi_base/Common/epi_acqui_stb.h"
#include "Imaging3T/epi_base/Common/epi_spoiler_stb.h"
#include "Imaging3T/epi_base/Common/epi_kordering.h"
#include "ImagingSeq/epi_base/Common/epi_common_recon_define.h"

//disable harmless compiler warning 
#if defined(WIN64)||defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4482)  //nonstandard extension used: enum 
#endif

using namespace Umr;
NAME_SPACE_UMR_SEQUENCE_BEGIN

EPIFIDKernel::EPIFIDKernel()
: m_pKOrdering(NULL)
//, m_pKSpace(NULL)
//, m_pSliceArray(NULL)
, m_dB0_Tesla(0.0)
, m_dFov_Read_mm(0.0)
, m_dFov_Phase_mm(0.0)
, m_lStepsFromKspaceStartToCenterAtInterval(0)
, m_lStepsFromKspaceStartToCenterAtPPARef(0)
, m_dGpeBlipIntervalMoment_mTus_m(0.0)
, m_dGpeBlipAmplAtInterval_mT_m(0.0)
, m_dPreGpeBlipMaxMoment_mTus_m(0.0)
, m_dPreGpeBlipMaxAmpl_mT_m(0.0)
, m_dPreGpeBlipMomentRefDivInterval(0.001)
, m_dMaxwellDivSquareZ(0.0)
, m_lCurrEchoIndex(0)
, m_lPPASampleFactor(2)
, m_adRotMatrixLog2Phy(9,0.0)
, m_lNumberOfDummy(0)
, m_bIsTriggerActive(true)
, m_bIsNeedTriggerOut(true)
, m_bIsNeedDummy(true)
, m_bIsGpeActive(true)
, m_bIsPrepared(false)
, m_lReasonPrepareFailed(kSuccess)
, m_lStepsPPADummy(0)
, m_bIsCorrDeltaPhase(false)
, m_eScanStatus(kImageMode)
, m_eSceneMode(kDefaultMode)
, m_ePPAScanMode(kUnipolarMode)
, m_pEPIPNSCalc(NULL)
{

}

EPIFIDKernel::~EPIFIDKernel()
{
    if (NULL != m_pEPIPNSCalc)
    {
        delete m_pEPIPNSCalc;
        m_pEPIPNSCalc = NULL;
    }

}

bool EPIFIDKernel::Initialize()
{

    if (NULL == (m_spExcitSTB = shared_ptr<EPIExcitSTB>(new EPIExcitSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPI ExcitSTB";
        return false;
    }
    if (!m_spExcitSTB->Initialize())
    {
        SEQ_LOG_ERROR << "EPI ExcitSTB is failed to init";
        return false;
    }
    if (NULL == (m_spPrePhaSTB = shared_ptr<EPIPrePhaseSTB>(new EPIPrePhaseSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPI PrephaseSTB";
        return false;
    }
    if (!m_spPrePhaSTB->Initialize())
    {
        SEQ_LOG_ERROR << "EPI PrePhaSTB is failed to init";
        return false;
    }
    if (NULL == (m_spAcquiSTB = shared_ptr<EPIAcquiSTB>(new EPIAcquiSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPI AcquisitionSTB";
        return false;
    }
    if (!m_spAcquiSTB->Initialize())
    {
        SEQ_LOG_ERROR << "EPI AcquiSTB is failed to init";
        return false;
    }
    if (NULL == (m_spEPISpoilerSTB = shared_ptr<EPISpoilerSTB>(new EPISpoilerSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPISpoilerSTB";
        return false;
    }
    if (!m_spEPISpoilerSTB->Initialize())
    {
        SEQ_LOG_ERROR << "EPISpoilerSTB is failed to init";
        return false;
    }
    if (NULL == (m_spTriggerFillTimeSTB = shared_ptr<FillTimeSTB>(new FillTimeSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPI FillTimeSTBAhead";
        return false;
    }
    if (NULL == (m_spFillTimeSTB = shared_ptr<FillTimeSTB>(new FillTimeSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPI FillTimeSTB";
        return false;
    }
    if (NULL == (m_spExpandFillTimeSTB = shared_ptr<FillTimeSTB>(new FillTimeSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPI FillTimeSTB";
        return false;
    }
    if (NULL == (m_spEndFillTimeSTB = shared_ptr<FillTimeSTB>(new FillTimeSTB())))
    {
        SEQ_LOG_ERROR << "failed to create EPI FillTimeSTB";
        return false;
    }
    if (NULL == ( m_spTriggerStb = shared_ptr<TriggerStb>(new TriggerStb())))
    {
        SEQ_LOG_ERROR << "failed to create EPI TriggerStb";
        return false;
    }
    if (!m_spTriggerStb->Initialize())
    {
        SEQ_LOG_ERROR << "EPI TriggerStb is failed to init";
        return false;
    }
    m_pEPIPNSCalc = new EPIPNSCalc();
    if (!m_pEPIPNSCalc->Initial())
    {
        SEQ_LOG_ERROR << "EPIPNSCalc Initail() fail";
        return false;
    }

    return true;
}

bool EPIFIDKernel::PreSettingAndCreation(const IUprotocol& rUprot)
{
    //-----------------------------------------------------------------------
    int64_t eScanScene = 0;
    GET_U_PARA_OPTION(rUprot,"Seq.App.Scene",&eScanScene);
    m_eSceneMode = eSceneMode(eScanScene);
    if (kHeadMode == m_eSceneMode )
    {
        m_ePPAScanMode = kBipolarMode;
    }
    else
    {
        m_ePPAScanMode = kUnipolarMode;
    }
    //-------------------------------------------------------------------
    //int64_t lContext = fGetSeqContext(rUprot);

    if (Umr::ePPAMethod::kKSpaceBased == m_pKSpace->GetPPAMethod())
    {
        if (Umr::ePPARefenceMode::kSeparated != m_pKSpace->GetPPARefenceMode())
        {
            if (!Context4Prep_IsValidRange())
            {
                SEQ_LOG_ERROR << "EPI PPA should be ExternalRef Mode";
            }
            return false;
        }
        if (kUnipolarMode == m_ePPAScanMode )
        {
            m_lPPASampleFactor = 2;
        }
        else if (kBipolarMode == m_ePPAScanMode )
        {
            m_lPPASampleFactor = 1;
        }
        else
        {
            m_lPPASampleFactor = 2;
        }
    }
    double dFlipAngle = 0.0;
    GET_U_PARA_DOUBLE(rUprot,"Seq.Basic.FlipAngle",&dFlipAngle);
    //---------------------------------------
    if (!m_spExcitSTB->CreateRf( "waveform.excit_slrlin_90deg_4tbw_200"
                               , MrTime(3.5,_ms)
                               , dFlipAngle
                               , 0.0))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "Rf Failed to be Created";
        }
        return false;
    }

    m_spExcitSTB->SetGroPreRefActive(true);
    m_spExcitSTB->SetGsRephaseRampScale(2.0);
    m_spExcitSTB->SetGrPreRefRampScale(1.0);

    //---------------------------------------
    double dGradMinRiseTime(10.0);
    double dGradMaxAmpl_mT_m(20.0);
    SA_GetSysParam4Seq(kdGrad_MinRiseTime_us_mT_m_Absolute, dGradMinRiseTime);
    SA_GetSysParam4Seq(kdGrad_MaxAmpl_mT_m_Absolute, dGradMaxAmpl_mT_m);

    //double dSystemMaxAmplScale =  1.0;
    //GET_U_PARA_DOUBLE(rUprot,"Seq.App.SystemMaxAmplScale",&dSystemMaxAmplScale);        //  prot default 1.15
    //dGradMaxAmpl_mT_m = dGradMaxAmpl_mT_m*dSystemMaxAmplScale;

    GradPerformance AcqSTBGradPerformance((dGradMaxAmpl_mT_m-1),dGradMinRiseTime);
    m_spAcquiSTB->SetGradPerformance(AcqSTBGradPerformance);

    double dGpeBlipRampScale =1.0;
    GET_U_PARA_DOUBLE(rUprot,"Seq.App.GpeBlipRampScale",&dGpeBlipRampScale);           //  prot  default  1.5
    m_spAcquiSTB->SetGpeBlipRampScale(dGpeBlipRampScale);
    double dGroRampScaleModified = 1.0;
    GET_U_PARA_DOUBLE(rUprot,"Seq.App.GroRampScaleModified",&dGroRampScaleModified);   // prot  default  1.3
    m_spAcquiSTB->SetGroRampScale(dGroRampScaleModified);

    GET_U_PARA_BOOL(rUprot,"Seq.App.GpeActive",&m_bIsGpeActive);
    m_spPrePhaSTB->SetPreGpeActive(m_bIsGpeActive);
    m_spAcquiSTB->SetGpeBlipActive(m_bIsGpeActive);
    int64_t iPartialEchoMode = 0 ;
    GET_U_PARA_OPTION(rUprot,"Seq.KSpace.PartialEcho",&iPartialEchoMode);
    double dPartialEchoFactor = 0.0;
    if (ePartialEchoMode::kPartialOff != ePartialEchoMode(iPartialEchoMode))
    {
        GET_U_PARA_DOUBLE(rUprot,"Seq.App.PartialEchoFactor",&dPartialEchoFactor);
    }
    m_spAcquiSTB->SetPartialEchoFactor(dPartialEchoFactor);
    //---------------------------------------
    double dOverSamplePE = 0.0;
    GET_U_PARA_DOUBLE(rUprot,"Seq.KSpace.OverSamplingPE",&dOverSamplePE);
    GET_U_PARA_DOUBLE(rUprot,"Seq.GLI.CommonPara.FOVro",&m_dFov_Read_mm);
    GET_U_PARA_DOUBLE(rUprot,"Seq.GLI.CommonPara.FOVpe",&m_dFov_Phase_mm);

    SA_GetSysParam4Seq(kdMagnet_NomFieldStrengthB0_T, m_dB0_Tesla);

    m_dFov_Phase_mm = m_dFov_Phase_mm * (1+dOverSamplePE/100.0);

    // always prepare for biggest GpeBlip in prepare section
    int64_t lStepsFromKspaceStartToCenter = 0;
    if (eKOrderingMode::kLinearDescending == m_pKSpace->GetKOrderingMode())
    {
        lStepsFromKspaceStartToCenter = m_pKOrdering->GetMaxLinesIndexPE()
                                      - m_pKOrdering->GetCenterLineIndexPE();
    } 
    else if (eKOrderingMode::kLinearAscending == m_pKSpace->GetKOrderingMode())
    {
        lStepsFromKspaceStartToCenter = m_pKOrdering->GetCenterLineIndexPE()
                                      - m_pKOrdering->GetMinLinesIndexPE();
    } 
    else
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "EPI do not support other  kordering Mode";
        }
        return false;
    }

    if ( Umr::ePPAMethod::kNone == m_pKSpace->GetPPAMethod()) 
    {
        // m_pKOrdering->GetPPAFactorPE() always .equ. 1;
        m_dGpeBlipIntervalMoment_mTus_m = 1.0e6 / (LAMO_CONSTANT_1H*m_dFov_Phase_mm);

        m_lStepsFromKspaceStartToCenterAtPPARef = 0;

        m_lStepsFromKspaceStartToCenterAtInterval = lStepsFromKspaceStartToCenter;
    } 
    else if (Umr::ePPAMethod::kKSpaceBased == m_pKSpace->GetPPAMethod())
    {
        m_dGpeBlipIntervalMoment_mTus_m = 1.0e6 * m_pKOrdering->GetPPAFactorPE()
                                        / (LAMO_CONSTANT_1H*m_dFov_Phase_mm);

        m_lStepsFromKspaceStartToCenterAtPPARef = m_pKOrdering->GetCenterLineIndexPE() // ref start index is always min
                                                - m_pKOrdering->GetPPARefLinesStartIndexPE();

        m_lStepsFromKspaceStartToCenterAtInterval = lStepsFromKspaceStartToCenter/ m_pKOrdering->GetPPAFactorPE();
    }
    else
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "do not have this type Acc Method";
        }
        return false;
    }


    m_dPreGpeBlipMomentRefDivInterval = static_cast<double>( m_lStepsFromKspaceStartToCenterAtPPARef )
                                      / static_cast<double>( m_lStepsFromKspaceStartToCenterAtInterval * m_pKOrdering->GetPPAFactorPE() );

    if ( m_lStepsFromKspaceStartToCenterAtInterval < 0
      || m_lStepsFromKspaceStartToCenterAtPPARef < 0)
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "Kordering calculation is not right";
        }
        return false;
    }

    ////////////////////////when only Bipolar Mode  PPA set Dummy /////////////////////////////////////
    if (Umr::ePPAMethod::kKSpaceBased == m_pKSpace->GetPPAMethod())
    {
        if (kBipolarMode == m_ePPAScanMode)
        {
            m_lStepsPPADummy = m_lStepsFromKspaceStartToCenterAtInterval - m_lStepsFromKspaceStartToCenterAtPPARef*m_lPPASampleFactor;
        }
        else
        {
            m_lStepsPPADummy = 0;
        }
    }
    else
    {
        m_lStepsPPADummy = 0;
    }
    if (0 > m_lStepsPPADummy)
    {
        m_lStepsPPADummy = 0;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////


    if (eKOrderingMode::kLinearDescending == m_pKSpace->GetKOrderingMode())
    {
        m_spAcquiSTB->SetGpeBlipMoment(-fabs(m_dGpeBlipIntervalMoment_mTus_m));
    } 
    else if (eKOrderingMode::kLinearAscending == m_pKSpace->GetKOrderingMode())
    {
        m_spAcquiSTB->SetGpeBlipMoment(fabs(m_dGpeBlipIntervalMoment_mTus_m));
    } 
    else
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "PPA_Acc method do not exist";
        }
        return false;
    }

    //m_spAcquiSTB->SetGroRampScale(1.0);

    //---------------------------------------
    m_spAcquiSTB->SetKSpacePtr(m_pKSpace);
    m_spAcquiSTB->SetKOrderingPtr(m_pKOrdering);
    //---------------------------------------
    m_spPrePhaSTB->SetPreGroRampScale(1.0);
    m_spPrePhaSTB->SetPreGpeRampScale(2.0);

    Slice rSlice = (*m_pSliceArray)[0];
    Umr::PatientPosition ePatientPos = HFS; 
    SA_GetSysParam4Seq(klPatientRegist_Position,(int64_t &)(ePatientPos));
    if (!CalcRotMatrixFromLogToPhy(rSlice,ePatientPos,&m_adRotMatrixLog2Phy))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "CalcRotMatrixFromLogToPhy is failed";
        }
        return false;
    }

    bool bIsRampAcquiActive = false;
    GET_U_PARA_BOOL(rUprot,"Seq.App.FreePara7",&bIsRampAcquiActive);
    m_spAcquiSTB->SetRampSamplingState(bIsRampAcquiActive);

    int64_t lROSamples = 0;
    int64_t lBandWidth = 0;

    Mrexam::IParameter* pPara = (rUprot)["Seq.Basic.BandWidth"]; 
    pPara->GetValueOfLong(&lBandWidth);

    pPara = (rUprot)["Seq.KSpace.MatrixRO"]; 
    pPara->GetValueOfLong(&lROSamples);

    m_spAcquiSTB->SetBandWidth(lBandWidth);
    m_spAcquiSTB->SetROSamples(lROSamples);
    m_spAcquiSTB->SetRotMatrixLog2Phy(m_adRotMatrixLog2Phy);
    return true;
}

bool EPIFIDKernel::RelationSetting()
{
    //---------------------------------------
    m_spExcitSTB->SetGroPreRefMoment(- fabs(m_spAcquiSTB->GetGroPreMoment()));
    m_spPrePhaSTB->SetPreGroMoment(- fabs(m_spAcquiSTB->GetGroMoment()/2));

    //---------------------------------------
    if ( m_dPreGpeBlipMomentRefDivInterval > 1.0 )
    {
        m_dPreGpeBlipMaxMoment_mTus_m = - static_cast<double>(m_lStepsFromKspaceStartToCenterAtPPARef)
                                      * m_spAcquiSTB->GetGpeBlipMoment() / m_pKOrdering->GetPPAFactorPE();
    } 
    else
    {
        m_dPreGpeBlipMaxMoment_mTus_m = - static_cast<double>(m_lStepsFromKspaceStartToCenterAtInterval)
                                      * m_spAcquiSTB->GetGpeBlipMoment();
    }
    double dGpePreMoment_PPAKernel = 0.0;
    dGpePreMoment_PPAKernel = - static_cast<double>(m_lStepsFromKspaceStartToCenterAtPPARef)
                              * m_spAcquiSTB->GetGpeBlipMoment() / m_pKOrdering->GetPPAFactorPE();
    m_spPrePhaSTB->SetPreGpeMoment(m_dPreGpeBlipMaxMoment_mTus_m);
    m_spExcitSTB->SetGpePreMoment(dGpePreMoment_PPAKernel);
    ///////////////////////////////////////////////////////////////////////////
    MrTime ExcitSTBPreGroFlattopTime;
    CalcFlattopTimeForPreGro( fabs(m_spAcquiSTB->GetGroPreMoment()),
                              fabs(m_spAcquiSTB->GetGroAmpl_mT_m()),
                              m_spAcquiSTB->GetGroRampUpTime(),
                              & ExcitSTBPreGroFlattopTime);
    m_spExcitSTB->SetPreGroRampUpTime(m_spAcquiSTB->GetGroRampUpTime());
    m_spExcitSTB->SetPreGroFlattopTime(ExcitSTBPreGroFlattopTime);

    MrTime PrePhaSTBPreGroFlattopTime;
    CalcFlattopTimeForPreGro( fabs(m_spAcquiSTB->GetGroMoment()),
                              fabs(m_spAcquiSTB->GetGroAmpl_mT_m()),
                              m_spAcquiSTB->GetGroRampUpTime(),
                              & PrePhaSTBPreGroFlattopTime);
    m_spPrePhaSTB->SetPreGroRampUpTime(m_spAcquiSTB->GetGroRampUpTime());
    m_spPrePhaSTB->SetPreGroFlattopTime(PrePhaSTBPreGroFlattopTime);

    double dMinRiseTimeExcit = 10.0;
    if (ExcitSTBPreGroFlattopTime > MrTime(0,_10us))
    {
        SA_GetSysParam4Seq(kdGrad_MinRiseTime_us_mT_m_Absolute, dMinRiseTimeExcit);
    }
    else
    {
        dMinRiseTimeExcit = m_spAcquiSTB->GetGroRampUpTime()(_us)/m_spAcquiSTB->GetGroAmpl_mT_m();
    }
    double dSystemMaxGrad_mT_m = 20;
    SA_GetSysParam4Seq(kdGrad_MaxAmpl_mT_m_Absolute, dSystemMaxGrad_mT_m);
    double dMaxPreGroAmpl_mT_m = (abs(m_spAcquiSTB->GetGroAmpl_mT_m())+0.5) < (dSystemMaxGrad_mT_m-1) ? (abs(m_spAcquiSTB->GetGroAmpl_mT_m())+0.5):(dSystemMaxGrad_mT_m-1) ;
    GradPerformance ExcitPreGroPerformance(dMaxPreGroAmpl_mT_m,dMinRiseTimeExcit);

    double dMinRiseTimePrePhase = 10.0;
    if (PrePhaSTBPreGroFlattopTime > MrTime(0,_10us))
    {
        SA_GetSysParam4Seq(kdGrad_MinRiseTime_us_mT_m_Absolute, dMinRiseTimePrePhase);
    }
    else
    {
        dMinRiseTimePrePhase = m_spAcquiSTB->GetGroRampUpTime()(_us)/m_spAcquiSTB->GetGroAmpl_mT_m();
    }
    GradPerformance PrePhasePreGroPerformance(dMaxPreGroAmpl_mT_m,dMinRiseTimePrePhase);

    double dGradNomRiseTime(10.0);
    double dGradNomAmpl_mT_m(15.0);
    SA_GetSysParam4Seq(kdGrad_NomRiseTime_us_mT_m, dGradNomRiseTime);
    //SA_GetSysParam4Seq(kdGrad_NomAmpl_mT_m, dGradNomAmpl_mT_m);
    GradPerformance NomGradPerformance(dGradNomAmpl_mT_m,dGradNomRiseTime);

    vector<GradPerformance> vec_ExcitGradPerformance,vec_PrePhaseGradPerformance;
    vec_ExcitGradPerformance.resize(3);
    vec_ExcitGradPerformance[0] = ExcitPreGroPerformance;      //RO
    vec_ExcitGradPerformance[1] = NomGradPerformance;         //PE
    vec_ExcitGradPerformance[2] = NomGradPerformance;         //SS
    vec_PrePhaseGradPerformance.resize(2);
    vec_PrePhaseGradPerformance[0] = PrePhasePreGroPerformance;   //RO
    vec_PrePhaseGradPerformance[1] = NomGradPerformance;      //PE

    m_spExcitSTB->SetGradPerformance(vec_ExcitGradPerformance);
    m_spPrePhaSTB->SetGradPerformance(vec_PrePhaseGradPerformance);
    ///////////////////////////////////////////////////////////////////////////////////
    return true;
}

bool EPIFIDKernel::PostSetting(const IUprotocol& rUprot)
{
    //int64_t lContext = fGetSeqContext(rUprot);

    //---------------------------------------
    m_spAcquiSTB->SetGroNcoAdditionalPhase(m_spExcitSTB->GetRfPhase_degree());
    m_spAcquiSTB->SetGroNcoFreqShift(0);

    //---------------------------------------
    m_dGpeBlipAmplAtInterval_mT_m = m_spAcquiSTB->GetGpeBlipAmpl_mT_m();
    m_dPreGpeBlipMaxAmpl_mT_m = m_spPrePhaSTB->GetPreGpeAmpl_mT_m();

    //---------------------------------------
    if (!Context4Prep_IsValidRange())
    {
        if (!CalcMaxwellCoefficient(m_adRotMatrixLog2Phy,&m_dMaxwellDivSquareZ))
        {
            return false;
        }
    }

    if (m_bIsNeedDummy)
    {
        GET_U_PARA_LONG(rUprot,"Seq.App.NumberOfDummyScan",&m_lNumberOfDummy);
    }
    else
    {
        m_lNumberOfDummy = 0;
    }

    if (m_bIsTriggerActive)
    {
        m_spTriggerStb->SetTriggerMode(kSyncTriggerExport);
        m_spTriggerFillTimeSTB->SetFillTime(m_spTriggerStb->GetDuration());
    }

    //TE
    int64_t lTE_us = 0;
    GET_U_PARA_LONG(rUprot,"Seq.Basic.TE",&lTE_us);
    m_TE = MrTime(static_cast<double>(lTE_us),_us);

    m_FillTimeAtIntervalLineAcqui = m_TE
                                  - m_spExcitSTB->GetTEAfterRf()
                                  - m_spAcquiSTB->GetDuration()
                                    * 3
                                  - m_spPrePhaSTB->GetDuration()
                                  - m_spAcquiSTB->GetDuration()
                                    * (static_cast<double>( m_lStepsFromKspaceStartToCenterAtInterval) + 0.5);
    m_FillTimeAtPPARefLineAcqui = m_FillTimeAtIntervalLineAcqui
                                + m_spAcquiSTB->GetDuration() * (static_cast<double>
                                ( m_lStepsFromKspaceStartToCenterAtInterval
                                - m_lStepsFromKspaceStartToCenterAtPPARef ));

    if ( m_FillTimeAtIntervalLineAcqui < MrTime(0.0)
      || m_FillTimeAtPPARefLineAcqui < MrTime(0.0))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "filltime is negtive ";
        }
        return false;
    }

    m_FillTimeAtIntervalLineAcqui.RoundUp(_10us);
    m_FillTimeAtPPARefLineAcqui.RoundUp(_10us);

    //m_TimeForExpand = m_spAcquiSTB->GetDuration() 
    //                  * static_cast<double>(( ( m_pKOrdering->GetScanLinesPE() - m_lStepsFromKspaceStartToCenterAtInterval )
    //                  - ( m_pKOrdering->GetPPARefLinesPE() - m_lStepsFromKspaceStartToCenterAtPPARef )));

    //filltime_stb->prepare() always true
    m_spFillTimeSTB->SetFillTime( m_FillTimeAtIntervalLineAcqui );

    //Duration
    MrTime ADC_RFDeadTime = SA_GetDeadTime(kDT_ADC_RF);
    m_spEndFillTimeSTB->SetFillTime(ADC_RFDeadTime);
    m_Duration = m_spExcitSTB->GetDuration()
                 + m_spAcquiSTB->GetDuration()*3
                 + m_spFillTimeSTB->GetDuration()
                 + m_spPrePhaSTB->GetDuration()
                 + m_spAcquiSTB->GetDuration()*static_cast<double>(m_pKOrdering->GetScanLinesPE())
                 + m_spEPISpoilerSTB->GetDuration()
                 + m_spEndFillTimeSTB->GetDuration();

    MrTime PPAKernelDuration;
    PPAKernelDuration = m_spExcitSTB->GetDuration()
                      + m_spAcquiSTB->GetDuration()*3
                      + m_spPrePhaSTB->GetDuration()
                      + m_spAcquiSTB->GetDuration()*static_cast<double>(m_pKOrdering->GetPPARefLinesPE()*m_lPPASampleFactor+m_lStepsPPADummy)
                      + m_spEPISpoilerSTB->GetDuration();
                          
    m_TimeForExpand = m_Duration - m_spEndFillTimeSTB->GetDuration() - PPAKernelDuration;

    m_spExpandFillTimeSTB->SetFillTime( MrTime(static_cast<double>(abs(m_TimeForExpand(_ns))),_ns) );
    if (m_TimeForExpand.GetValue(_ns) < 0)
    {
        m_Duration = m_Duration - m_TimeForExpand;
    } 
    if(m_bIsTriggerActive)
    {
        m_Duration += m_spTriggerStb->GetDuration();
    }
    return true;
}

bool EPIFIDKernel::SetParaToUprotocol(IUprotocol& rUprot)
{
    bool bIRIP_OffCenterRO = false;
    GET_U_PARA_BOOL(rUprot,"Seq.App.IRIP_OffCenterRO",&bIRIP_OffCenterRO);
    SET_U_PARA_BOOL(rUprot,"IRIP.FromSeq.App.Use_IRIP_OffCenterRO",bIRIP_OffCenterRO);
    double dFov_Read_mm = 0.0;
    double dOffcenter_RO_Relative = 0.0;
    GET_U_PARA_DOUBLE(rUprot,"Seq.GLI.CommonPara.FOVro",&dFov_Read_mm);
    SliceArray SliceArrayForOffcenter = *m_pSliceArray;
    double dSliceOffcenterRO = SliceArrayForOffcenter[0].GetSliceOffCenterRO_mm();
    dOffcenter_RO_Relative = static_cast<double>( (dSliceOffcenterRO)/dFov_Read_mm/2 );
    SET_U_PARA_DOUBLE(rUprot,"IRIP.FromSeq.App.OffCenterRO_Relative",dOffcenter_RO_Relative);
    //APP
    SET_U_PARA_DOUBLE(rUprot,"ScanInfo.EchoSpacing",m_spAcquiSTB->GetEchoSpaceTime()(_ms,kDouble));

    //Maxwell
    SET_U_PARA_DOUBLE(rUprot,"IRIP.FromSeq.App.MaxwellDivSquareZ",m_dMaxwellDivSquareZ);

    //PPA_Average////////////////////////////////////////////////////////////////////////////////////
    int64_t lPPASeparatedAverage = 1;
    if (kBipolarMode==m_ePPAScanMode)
    {
        lPPASeparatedAverage = 2;
    }
    else
    {
        lPPASeparatedAverage = 1;
    }
    SET_U_PARA_LONG(rUprot,"Seq.PPA.SeparatedAverage",lPPASeparatedAverage);
    //////////////////////////////////////////////////////////////////////////////////////////////////

    //IRIP
    m_spAcquiSTB->SetRampSampInfo(rUprot);
    SET_U_PARA_LONG(rUprot,"IRIP.FromSeq.KSpace.FTLengthPE",m_pKOrdering->GetFTLengthPE());
    SET_U_PARA_LONG(rUprot,"IRIP.FromSeq.KSpace.PartialPELines",m_pKOrdering->GetPartialLinesPE());
    if (eKOrderingMode::kLinearAscending == m_pKSpace->GetKOrderingMode())
    {
        SET_U_PARA_BOOL(rUprot,"IRIP.FromSeq.KSpace.PartialPEReverse",false);
    } 
    else if (eKOrderingMode::kLinearDescending == m_pKSpace->GetKOrderingMode())
    {
        SET_U_PARA_BOOL(rUprot,"IRIP.FromSeq.KSpace.PartialPEReverse",true);
    }
    else
    {
        SEQ_LOG_ERROR << "kordering mode is not support";
        return false;
    }
    string sSeqName;
    sSeqName = GetDefSeqDisplayName(rUprot,false);
    SET_U_PARA_STRING(rUprot,"ScanInfo.SeqDisplayName",sSeqName);
    if (Umr::ePPAMethod::kNone == m_pKSpace->GetPPAMethod())
    {
        SET_U_PARA_LONG(rUprot,"IRIP.FromSeq.PPA.PPAFactorPE", 1);
    } 
    else
    {
        SET_U_PARA_LONG(rUprot,"IRIP.FromSeq.PPA.PPAFactorPE",static_cast<int64_t>(m_pKOrdering->GetPPAFactorPE()));
    }
    SET_U_PARA_LONG(rUprot,"IRIP.FromSeq.PPA.FirstRefLineIdxPE",m_pKOrdering->GetPPARefLinesStartIndexPE());
    SET_U_PARA_LONG(rUprot,"IRIP.FromSeq.PPA.RefLineLengthPE",m_pKOrdering->GetPPARefLinesPE());

    //SET_U_PARA_STRING(rUprot,"IRIP.PipeLineConfig.Procedure","UmrIripEPIProcedure");

    return true;
}

bool EPIFIDKernel::Prepare(const IUprotocol& rUprot)
{
    //int64_t lContext = fGetSeqContext(rUprot);
    m_bIsPrepared = false;
    if (!PreSettingAndCreation(rUprot))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "failed to set the relationship between STBs at begin";
        }
        return false;
    }

    if (!m_spAcquiSTB->Prepare(rUprot))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "AcquiSTB prepare failed";
        }
        m_lReasonPrepareFailed = kAcqSTB;
        return false;
    }

    if (!RelationSetting())
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "failed to set the relationship between STBs at end";
        }
        return false;
    }

    if (!m_spExcitSTB->Prepare(rUprot))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "ExcitSTB prepare failed";
        }
        m_lReasonPrepareFailed = kExcitSTB;
        return false;
    }

    if (!m_spPrePhaSTB->Prepare(rUprot))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "PrePhaSTB prepare failed";
        }
        return false;
    }
    if (!m_spEPISpoilerSTB->Prepare(rUprot))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "EPISpoilerSTB prepare failed";
        }
        return false;
    }
    if (!PostSetting(rUprot))
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "prepare not finished!";
        }
        m_lReasonPrepareFailed = kTE;
        return false;
    }

    if (!PrepareKernelGradInfo4RMS())
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "PrepareKernelGradInfo4RMS not finished!";
        }
        return false;
    }

    m_bIsPrepared=true;
    m_lReasonPrepareFailed = kSuccess;
    return true;
}

bool EPIFIDKernel::BeginCheck()
{
    return true;
}

bool EPIFIDKernel::EndCheck()
{
    return true;
}

void EPIFIDKernel::SetDHL( const IUprotocol& /*rUprot*/
                         , const Slice& rSlice
                         , const uint16_t& lCurrentSlice 
                         , const AcqIndex& rAcqIndex)
{
    DHL* pDHL = m_spAcquiSTB->GetADCPtr()->GetDHL();

    pDHL->RemoveCtrlFlag(DHL_PHASE_CORRECTION);
    pDHL->RemoveCtrlFlag(DHL_PPA_REFLINE);
    pDHL->RemoveCtrlFlag(DHL_READOUT_REVERSION);
    pDHL->RemoveCtrlFlag(DHL_FIRST_LINE_SLICE);
    pDHL->RemoveCtrlFlag(DHL_LAST_LINE_SLICE);
    pDHL->RemoveCtrlFlag(DHL_EPI_RAMPSAMPLING);
    pDHL->RemoveCtrlFlag(DHL_EPI_PHASE_CORRECTION_PEOFF);

    FillDhl4LciAndSlice(rAcqIndex,rSlice,pDHL);
    pDHL->SetLciSlice(lCurrentSlice);

    if (kScanStatus_Imaging == rAcqIndex.GetScanStatus())
    {
        if (m_lNumberOfDummy <= rAcqIndex.GetRepeatIndex())
        {
            pDHL->SetLciRepeat(rAcqIndex.GetRepeatIndex()-static_cast<uint16_t>(m_lNumberOfDummy));
        }
    }

    pDHL->SetPECenter(static_cast<uint16_t>(m_pKOrdering->GetCenterLineIndexPE()));
    pDHL->SetROCenter(static_cast<uint16_t>(m_spAcquiSTB->GetSamplesAfterRegrid()-m_pKSpace->GetMatrixRO()/2));
    if (m_spAcquiSTB->GetActualAdcSamples() != m_spAcquiSTB->GetSamplesAfterRegrid())
    {
        pDHL->AddCtrlFlag(DHL_EPI_RAMPSAMPLING);
    }
    // PPA dDummy///////////////////////////////////////////////////////////////
    int64_t m_lStepsPPADummy_temp = 0;
    if (m_eScanStatus == kPPAMode)
    {
        m_lStepsPPADummy_temp = m_lStepsPPADummy;
    }
    else
    {
        m_lStepsPPADummy_temp = 0;
    }
    //////////////////////////////////////////////////////////////////////////////
    //// CtrlFlag
    //if (KScanStatus_Prep == rAcqIndex.GetScanStatus())
    //{
    //    pDHL->AddCtrlFlag(DHL_PPA_REFLINE); // PPA Ref
    //} 
    if (m_spAcquiSTB->GetGroAmpl_mT_m() < 0.0)
    {
        pDHL->AddCtrlFlag(DHL_READOUT_REVERSION); // RO Reverse
    }

    // DataInfo
    if (m_eScanStatus == kPPAMode)
    {
        pDHL->AddCtrlFlag(DHL_PPA_REFLINE); // PPA Ref
        if (eKOrderingMode::kLinearAscending == m_pKSpace->GetKOrderingMode())
        {
            pDHL->SetLciPELine( static_cast<uint16_t>(
                                m_pKOrdering->GetPPARefLinesStartIndexPE()
                             + (m_lCurrEchoIndex-m_lStepsPPADummy)/m_lPPASampleFactor));
        } 
        else if (eKOrderingMode::kLinearDescending == m_pKSpace->GetKOrderingMode())
        {
            pDHL->SetLciPELine( static_cast<uint16_t>(
                                m_pKOrdering->GetPPARefLinesStartIndexPE()
                              + m_pKOrdering->GetPPARefLinesPE() - 1
                              - (m_lCurrEchoIndex-m_lStepsPPADummy)/m_lPPASampleFactor));
        }
    } 
    else if (m_eScanStatus == kImageMode)
    {
        pDHL->SetLciPELine(static_cast<uint16_t>(m_pKOrdering->GetLinesIndexPE(rAcqIndex.GetShotIndex(),m_lCurrEchoIndex)));
    }
    else if (m_eScanStatus == kCorrectMode)
    {
        pDHL->AddCtrlFlag(DHL_EPI_PHASE_CORRECTION_PEOFF); 
        pDHL->SetLciPELine(static_cast<uint16_t>(m_pKOrdering->GetLinesIndexPE(rAcqIndex.GetShotIndex(),m_lCurrEchoIndex)));
    }

    pDHL->SetLciSegment(static_cast<uint16_t>(m_lCurrEchoIndex-m_lStepsPPADummy_temp));

    pDHL->SetLciSPELine(0);

    if (m_eScanStatus == kCorrectMode || m_eScanStatus == kImageMode)
    {
        if(0 == rAcqIndex.GetShotIndex() && 0 == m_lCurrEchoIndex)
        {
            pDHL->AddCtrlFlag(DHL_FIRST_LINE_SLICE);
        }
        if(rAcqIndex.GetShotIndex() == m_pKOrdering->GetShots() - 1 
            && m_lCurrEchoIndex == m_pKOrdering->GetScanLinesPE() - 1)
        {
            pDHL->AddCtrlFlag(DHL_LAST_LINE_SLICE);
        }
    }
}

void EPIFIDKernel::SetPhaseCorrDHL( const IUprotocol& /*rUprot*/
                                  , const Slice& rSlice
                                  , const uint16_t& lCurrentSlice
                                  , const AcqIndex& rAcqIndex
                                  , const int64_t& lCurrEchoIndex)
{
    DHL* pDHL = m_spAcquiSTB->GetADCPtr()->GetDHL();

    pDHL->RemoveCtrlFlag(DHL_PHASE_CORRECTION);
    pDHL->RemoveCtrlFlag(DHL_PPA_REFLINE);
    pDHL->RemoveCtrlFlag(DHL_READOUT_REVERSION);
    pDHL->RemoveCtrlFlag(DHL_FIRST_LINE_SLICE);
    pDHL->RemoveCtrlFlag(DHL_LAST_LINE_SLICE);
	pDHL->RemoveCtrlFlag(DHL_EPI_RAMPSAMPLING);
    pDHL->RemoveCtrlFlag(DHL_EPI_PHASE_CORRECTION_PEOFF);

    FillDhl4LciAndSlice(rAcqIndex,rSlice,pDHL);
    pDHL->SetLciSlice(lCurrentSlice);

    if (kScanStatus_Imaging == rAcqIndex.GetScanStatus())
    {
        if (m_lNumberOfDummy <= rAcqIndex.GetRepeatIndex())
        {
            pDHL->SetLciRepeat(rAcqIndex.GetRepeatIndex()-static_cast<uint16_t>(m_lNumberOfDummy));
        }
    }

    pDHL->SetPECenter(static_cast<uint16_t>(m_pKOrdering->GetCenterLineIndexPE()));
    pDHL->SetROCenter(static_cast<uint16_t>(m_spAcquiSTB->GetSamplesAfterRegrid()-m_pKSpace->GetMatrixRO()/2));
    if (m_spAcquiSTB->GetActualAdcSamples() != m_spAcquiSTB->GetSamplesAfterRegrid())
    {
        pDHL->AddCtrlFlag(DHL_EPI_RAMPSAMPLING);
    }

    // CtrlFlag
    pDHL->AddCtrlFlag(DHL_PHASE_CORRECTION); // Phase Correct
    if (m_eScanStatus == kPPAMode)
    {
        pDHL->AddCtrlFlag(DHL_PPA_REFLINE); // PPA Ref
    } 
    if (m_eScanStatus == kCorrectMode)
    {
        pDHL->AddCtrlFlag(DHL_EPI_PHASE_CORRECTION_PEOFF); 
    }
    if (m_spAcquiSTB->GetGroAmpl_mT_m() < 0.0)
    {
        pDHL->AddCtrlFlag(DHL_READOUT_REVERSION); // RO Reverse
    }
    // DataInfo
    pDHL->SetLciPELine(static_cast<uint16_t>(m_pKOrdering->GetCenterLineIndexPE()));
    pDHL->SetLciSPELine(0);
    pDHL->SetLciSegment(static_cast<uint16_t>(lCurrEchoIndex)); // 3 ref scans
}

bool EPIFIDKernel::Execute( const IUprotocol& rUprot
                          , const AcqIndex& rAcqIndex)
{
    if (Context4Exec_IsGradCheck())
    {
        return Check(rUprot);
    }

    m_eScanStatus = GetScanStatus(rAcqIndex);
    m_eSampleMode = GetSampleMode(rAcqIndex);
    //map slice index to scan order
    int64_t lCurrentSlice = 0;
    lCurrentSlice = GetCurrentSlice(rUprot,m_pKSpace->GetMultiSliceMode(),m_pSliceArray->size(),rAcqIndex.GetSliceIndex());
    const Slice& rCurrentSlice = (*m_pSliceArray)[lCurrentSlice];

    if ( m_eScanStatus == kPPAMode)
    {
        m_spAcquiSTB->SetGpeBlipAmpl_mT_m(m_dGpeBlipAmplAtInterval_mT_m / static_cast<double>(m_pKOrdering->GetPPAFactorPE()));
        if ( m_dPreGpeBlipMomentRefDivInterval > 1.0 )
        {
            m_spPrePhaSTB->SetPreGpeAmpl_mT_m(m_dPreGpeBlipMaxAmpl_mT_m);
        } 
        else
        {
            m_spPrePhaSTB->SetPreGpeAmpl_mT_m(m_dPreGpeBlipMaxAmpl_mT_m * m_dPreGpeBlipMomentRefDivInterval);
        }

        m_spFillTimeSTB->SetFillTime(m_FillTimeAtPPARefLineAcqui);
        m_spExcitSTB->SetGpePreActive(false);
        if (m_bIsTriggerActive)
        {
            m_spTriggerFillTimeSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
        }
    }
    else 
    {
        m_spAcquiSTB->SetGpeBlipAmpl_mT_m(m_dGpeBlipAmplAtInterval_mT_m);
        if ( m_dPreGpeBlipMomentRefDivInterval > 1.0 )
        {
            m_spPrePhaSTB->SetPreGpeAmpl_mT_m(m_dPreGpeBlipMaxAmpl_mT_m / m_dPreGpeBlipMomentRefDivInterval);
        } 
        else
        {
            m_spPrePhaSTB->SetPreGpeAmpl_mT_m(m_dPreGpeBlipMaxAmpl_mT_m);
        }
        m_spFillTimeSTB->SetFillTime(m_FillTimeAtIntervalLineAcqui);

        m_spExcitSTB->SetGpePreActive(false);
        if (m_bIsTriggerActive)
        {
            if((m_lNumberOfDummy==rAcqIndex.GetRepeatIndex())&&(true == m_bIsNeedTriggerOut))
            {
                m_spTriggerStb->Execute(rUprot,rCurrentSlice,rAcqIndex);
                m_bIsNeedTriggerOut = false;
            }
            else
            {
                m_spTriggerFillTimeSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
            }
        }
    }
    //else
    //{
    //    if (!Context4Prep_IsValidRange())
    //    {
    //        SEQ_LOG_ERROR << "the scan type is not defined";
    //    }
    //    return false;
    //}

    if (eKOrderingMode::kLinearAscending == m_pKSpace->GetKOrderingMode())
    {
        m_spPrePhaSTB->SetPreGpeAmpl_mT_m( - fabs(m_spPrePhaSTB->GetPreGpeAmpl_mT_m()) );
    } 
    else if (eKOrderingMode::kLinearDescending == m_pKSpace->GetKOrderingMode())
    {
        m_spPrePhaSTB->SetPreGpeAmpl_mT_m( fabs(m_spPrePhaSTB->GetPreGpeAmpl_mT_m()) );
    }
    else
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "now only support ascend or descend kordering";
        }
        return false;
    }

    if (!IsPrepared())
    {
        SEQ_LOG_ERROR << "prepare failed";
        return false;
    }

    if (!BeginCheck())
    {
        SEQ_LOG_ERROR << "some variable reset failed!";
        return false;
    }

    //---------------------------------------
    //TODO (xianwang.jiang@united-imaging.com) refocus and excit_rf nco_freq

    //m_spExcitSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);

    m_spAcquiSTB->SetCompensationGroActive(false);    //GradCompensation is active when truely execute
    m_spAcquiSTB->SetCompensationGpeActive(false);
    m_spAcquiSTB->SetGpeBlipActive(false);

        //initial: negative amplitude

    if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
    {
        m_spAcquiSTB->ReverseGroAmpl();
    }
    if ((kScanStatus_Imaging == rAcqIndex.GetScanStatus()))
    {
        if(m_lNumberOfDummy>rAcqIndex.GetRepeatIndex()) 
        {
            m_spAcquiSTB->SetAdcActive(false);
        }
        else 
        {
            m_spAcquiSTB->SetAdcActive(true); 
        }
    }
    else if((KScanStatus_Prep == rAcqIndex.GetScanStatus()))
    {
        m_spAcquiSTB->SetAdcActive(true); 
    }

    int64_t lRefEchoNum = 0;
    while (lRefEchoNum)
    {
        m_spAcquiSTB->ReverseGroAmpl();
        m_spAcquiSTB->SetPeLineIndex(0);
        if (Context4Exec_IsNormal())
        {
            SetPhaseCorrDHL(rUprot,rCurrentSlice,static_cast<uint16_t>(lCurrentSlice),rAcqIndex,3-lRefEchoNum);
        }
        if (1 == lRefEchoNum)
        {
            m_spAcquiSTB->SetCompensationGpeActive(false);             // the end STB is not need to compensation
        }
        m_spAcquiSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
        --lRefEchoNum;
    }
    if (m_eScanStatus == kImageMode || m_eScanStatus == kCorrectMode)  //
    {
        m_spFillTimeSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
    } 

    if (m_eScanStatus == kCorrectMode)
    {
        m_spPrePhaSTB->SetPreGpeActive(false);
    }
    else
    {
        m_spPrePhaSTB->SetPreGpeActive(m_bIsGpeActive);
    }
    ///----------------------default----------------------------
    if (0 < m_spPrePhaSTB->GetPreGroAmpl_mT_m())
    {
        m_spPrePhaSTB->ReversePreGroAmpl();
    }
    //--------------------------- PPA ---------------------------------------------
    if (kBipolarMode == m_ePPAScanMode)
    {
        if ( kPPA_EvenMode == m_eSampleMode )
        {
            if (0 == m_lStepsPPADummy%2)
            {
                if (0 > m_spPrePhaSTB->GetPreGroAmpl_mT_m())
                {
                    m_spPrePhaSTB->ReversePreGroAmpl();
                }
            }
            else
            {
                if (0 < m_spPrePhaSTB->GetPreGroAmpl_mT_m())
                {
                    m_spPrePhaSTB->ReversePreGroAmpl();
                }
            }
        }
        else if(kPPA_OddMode == m_eSampleMode)
        {
            if (0 == m_lStepsPPADummy%2)
            {
                if (0 < m_spPrePhaSTB->GetPreGroAmpl_mT_m())
                {
                    m_spPrePhaSTB->ReversePreGroAmpl();
                }
            }
            else
            {
                if (0 > m_spPrePhaSTB->GetPreGroAmpl_mT_m())
                {
                    m_spPrePhaSTB->ReversePreGroAmpl();
                }
            }
        }
    }
    //--------------------------------------------------------------------
    //m_spPrePhaSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);

    int64_t lTotalEchoNum = 0;
	int64_t lTotalPeNum = 128;
    if (Umr::ePPAMethod::kNone == m_pKSpace->GetPPAMethod())
    {
        lTotalEchoNum = m_pKOrdering->GetScanLinesPE();
    } 
    else if (Umr::ePPAMethod::kKSpaceBased == m_pKSpace->GetPPAMethod())
    {
        if (m_eScanStatus == kPPAMode)
        {
            lTotalEchoNum = m_pKOrdering->GetPPARefLinesPE()*m_lPPASampleFactor+ m_lStepsPPADummy;   //m_lStepsPPADummy is not zeros when only biploar mode ;
        } 
        else if (m_eScanStatus ==kImageMode || m_eScanStatus == kCorrectMode)
        {
            lTotalEchoNum = m_pKOrdering->GetScanLinesPE();
        } 
        else
        {
            if (!Context4Prep_IsValidRange())
            {
                SEQ_LOG_ERROR << "no such Scan Status";
            }
            return false;
        }
    } 
    else
    {
        if (!Context4Prep_IsValidRange())
        {
            SEQ_LOG_ERROR << "acc method is not support";
        }
        return false;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
    {
        m_spAcquiSTB->ReverseGroAmpl();
    }

    if (kBipolarMode == m_ePPAScanMode)
    {
        if (kPPA_EvenMode == m_eSampleMode )
        {
            if (0 == m_lStepsPPADummy%2)
            {
                if (m_spAcquiSTB->GetGroAmpl_mT_m() < 0)
                {
                    m_spAcquiSTB->ReverseGroAmpl();
                }
            } 
            else
            {
                if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
                {
                    m_spAcquiSTB->ReverseGroAmpl();
                }
            }
        } 
        else if (kPPA_OddMode == m_eSampleMode)
        {
            if (0 == m_lStepsPPADummy%2)
            {
                if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
                {
                    m_spAcquiSTB->ReverseGroAmpl();
                }
            } 
            else
            {
                if (m_spAcquiSTB->GetGroAmpl_mT_m() < 0)
                {
                    m_spAcquiSTB->ReverseGroAmpl();
                }
            }
        }
    }
    //initial: negative amplitude      // scan b=0 && b average index is even  initial: positive amplitude
    if (kHeadMode == m_eSceneMode)
    {
        if (kB0_EvenMode == m_eSampleMode )
        {
            if (m_spAcquiSTB->GetGroAmpl_mT_m() < 0)
            {
                m_spAcquiSTB->ReverseGroAmpl();
            }
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //m_spAcquiSTB->SetGpeBlipActive(m_bIsGpeActive);
	for(m_lCurrPeIndex =0; m_lCurrPeIndex<lTotalPeNum; ++m_lCurrPeIndex)
	{
		m_spExcitSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
		//m_spFillTimeSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
		//m_spPrePhaSTB->SetPreGpeMoment(m_dPreGpeBlipMaxMoment_mTus_m*m_lCurrPeIndex/lTotalPeNum);
		//m_spPrePhaSTB->SetGradPerformance(vec_PrePhaseGradPerformance);
		m_spPrePhaSTB->SetPreGpeAmpl_mT_m(m_dPreGpeBlipMaxAmpl_mT_m*2*(m_lCurrPeIndex-lTotalPeNum/2)/lTotalPeNum);
		m_spPrePhaSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
		for(m_lCurrEchoIndex =0; m_lCurrEchoIndex<lTotalEchoNum; ++m_lCurrEchoIndex)
		{
			m_spAcquiSTB->ReverseGroAmpl();
			if (m_eScanStatus == kPPAMode)
			{
				if (eKOrderingMode::kLinearAscending == m_pKSpace->GetKOrderingMode())
				{
					m_spAcquiSTB->SetPeLineIndex( m_pKOrdering->GetPPARefLinesStartIndexPE()
											   + (m_lCurrEchoIndex-m_lStepsPPADummy)/m_lPPASampleFactor
												- m_pKOrdering->GetCenterLineIndexPE());
				} 
				else if (eKOrderingMode::kLinearDescending == m_pKSpace->GetKOrderingMode())
				{
					m_spAcquiSTB->SetPeLineIndex( m_pKOrdering->GetPPARefLinesStartIndexPE()
												+ m_pKOrdering->GetPPARefLinesPE() - 1
												- (m_lCurrEchoIndex-m_lStepsPPADummy)/m_lPPASampleFactor
												- m_pKOrdering->GetCenterLineIndexPE());
				} 
				if (kUnipolarMode == m_ePPAScanMode)
				{
					m_spAcquiSTB->SetCompensationGpeActive(false);         // PPA only use the Grad positive lines , compensation is not need
					if (0 < m_spAcquiSTB->GetGroAmpl_mT_m())
					{
						m_spAcquiSTB->SetAdcActive(true);
						m_spAcquiSTB->SetGpeBlipActive(false);
					}
					else
					{
						m_spAcquiSTB->SetAdcActive(false);
						m_spAcquiSTB->SetGpeBlipActive(m_bIsGpeActive);
					}
				} 
				else if(kBipolarMode == m_ePPAScanMode)      //m_lStepsPPADummy is not zeros when only biploar mode 
				{
					if (m_lCurrEchoIndex < m_lStepsPPADummy)
					{
						m_spAcquiSTB->SetAdcActive(false);
						m_spAcquiSTB->SetGpeBlipActive(false);
					}
					else
					{
						m_spAcquiSTB->SetCompensationGpeActive(true);          // control by the config file
						m_spAcquiSTB->SetGpeBlipActive(m_bIsGpeActive);
						m_spAcquiSTB->SetAdcActive(true);
					}
				}
			} 
			else if (m_eScanStatus == kImageMode)
			{
				m_spAcquiSTB->SetPeLineIndex( m_pKOrdering->GetLinesIndexPE(rAcqIndex.GetShotIndex(),m_lCurrEchoIndex)
											- m_pKOrdering->GetCenterLineIndexPE());
				if(m_lNumberOfDummy>rAcqIndex.GetRepeatIndex()) 
				{
					m_spAcquiSTB->SetAdcActive(false);
				}
				else 
				{
					m_spAcquiSTB->SetAdcActive(true); 
				}
				m_spAcquiSTB->SetGpeBlipActive(m_bIsGpeActive);
				m_spAcquiSTB->SetCompensationGpeActive(true);         // control by the config file
			}
			else if (m_eScanStatus == kCorrectMode)
			{
				m_spAcquiSTB->SetPeLineIndex(0);
				m_spAcquiSTB->SetGpeBlipActive(false);
				m_spAcquiSTB->SetCompensationGpeActive(true);         // control by the config file
				m_spAcquiSTB->SetAdcActive(true);
			}

			if (lTotalEchoNum - 1 == m_lCurrEchoIndex)
			{
				m_spAcquiSTB->SetGpeBlipActive(false);
				m_spAcquiSTB->SetCompensationGpeActive(false);        // the end STB is not need to compensation
			}
			if (Context4Exec_IsNormal())
			{
				SetDHL(rUprot,rCurrentSlice,static_cast<uint16_t>(lCurrentSlice),rAcqIndex);
			}
			m_spAcquiSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
		}
    if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
    {
        if (m_spEPISpoilerSTB->GetGroAmpl_mT_m() < 0)
        {
            m_spEPISpoilerSTB->ReverseSpoilerGro();
        }
    }
    else
    {
        if (m_spEPISpoilerSTB->GetGroAmpl_mT_m() > 0)
        {
            m_spEPISpoilerSTB->ReverseSpoilerGro();
        }
    }
    m_spEPISpoilerSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
    if (m_TimeForExpand.GetValue(_ns) < 0)
    {
        if (m_eScanStatus == kImageMode || m_eScanStatus == kCorrectMode)
        {
            m_spExpandFillTimeSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
        }
    } 
    else
    {
        if (m_eScanStatus == kPPAMode)
        {
            m_spExpandFillTimeSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
        }
    }
	}
    m_spEndFillTimeSTB->Execute(rUprot,rCurrentSlice,rAcqIndex);
    if (!EndCheck())
    {
        SEQ_LOG_ERROR << "some variable reset failed!";
        return false;
    }

    return true;

}
//////////////////////////////////////////////////////////////////////////
///defined by itself
//////////////////////////////////////////////////////////////////////////
void EPIFIDKernel::SetKOrderingPtr(EPIKOrdering* pKOrdering)
{
    m_pKOrdering = pKOrdering;
}
void EPIFIDKernel::SetNumberOfDummy(int64_t _var)
{
    m_lNumberOfDummy=_var;
}
void EPIFIDKernel::SetTriggerActive(bool _var)
{
    m_bIsTriggerActive=_var;
}
void EPIFIDKernel::SetNeedDummy(bool _var)
{
    m_bIsNeedDummy=_var;
}
//---------------------------------------
bool EPIFIDKernel::Check(const IUprotocol& rUprot)
{
    AcqIndex myAcqIndex;
    m_spAcquiSTB->SetGpeBlipAmpl_mT_m(m_dGpeBlipAmplAtInterval_mT_m);
    m_spPrePhaSTB->SetPreGpeAmpl_mT_m(m_dPreGpeBlipMaxAmpl_mT_m);
    m_spFillTimeSTB->SetFillTime(m_FillTimeAtIntervalLineAcqui);

    //---------------------------------------
    //TODO (xianwang.jiang@united-imaging.com) refocus and excit_rf nco_freq

    const Slice& rCurrentSlice = (*m_pSliceArray)[0];

    m_spExcitSTB->Execute(rUprot,rCurrentSlice,myAcqIndex);
    m_spAcquiSTB->SetCompensationGroActive(false);     //GradCompensation is not active when gradient check 
    m_spAcquiSTB->SetCompensationGpeActive(false);
    //initial: negative amplitude
    if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
    {
        m_spAcquiSTB->ReverseGroAmpl();
    }
    int64_t lRefEchoNum = 3;
    m_spAcquiSTB->SetGpeBlipActive(false);
    while (lRefEchoNum)
    {
        m_spAcquiSTB->ReverseGroAmpl();
        m_spAcquiSTB->Execute(rUprot,rCurrentSlice,myAcqIndex);
        --lRefEchoNum;
    }

    m_spFillTimeSTB->Execute(rUprot,rCurrentSlice,myAcqIndex);

    m_spPrePhaSTB->Execute(rUprot,rCurrentSlice,myAcqIndex);

    //initial: negative amplitude
    if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
    {
        m_spAcquiSTB->ReverseGroAmpl();
    }

    int64_t lTotalEchoNum = 0;
    lTotalEchoNum = m_pKOrdering->GetScanLinesPE();

    m_spAcquiSTB->SetGpeBlipActive(true);
    for(m_lCurrEchoIndex =0; m_lCurrEchoIndex<lTotalEchoNum; ++m_lCurrEchoIndex)
    {
        m_spAcquiSTB->SetPeLineIndex( m_pKOrdering->GetLinesIndexPE(0,m_lCurrEchoIndex)
                    - m_pKOrdering->GetCenterLineIndexPE());

        m_spAcquiSTB->ReverseGroAmpl();
        if (lTotalEchoNum - 1 == m_lCurrEchoIndex)
        {
            m_spAcquiSTB->SetGpeBlipActive(false);
        }
        m_spAcquiSTB->Execute(rUprot,rCurrentSlice,myAcqIndex);
    }
    if (m_spAcquiSTB->GetGroAmpl_mT_m() > 0)
    {
        if (m_spEPISpoilerSTB->GetGroAmpl_mT_m() < 0)
        {
            m_spEPISpoilerSTB->ReverseSpoilerGro();
        }
    }
    else
    {
        if (m_spEPISpoilerSTB->GetGroAmpl_mT_m() > 0)
        {
            m_spEPISpoilerSTB->ReverseSpoilerGro();
        }
    }
    m_spEPISpoilerSTB->Execute(rUprot,rCurrentSlice,myAcqIndex);

    return true;
}
MrTime EPIFIDKernel::GetDuration()const
{
    return m_Duration;
}
int64_t EPIFIDKernel::GetNumberOfDummy()const
{
    return m_lNumberOfDummy;
}
bool EPIFIDKernel::IsPrepared()const
{
    return m_bIsPrepared;
}
//---------------------------------------
double EPIFIDKernel::GetEnergy_J()const
{
    return m_spExcitSTB->GetEnergy_J();
}
bool EPIFIDKernel::HasExcitation()const
{
    return true;
}
MrTime EPIFIDKernel::GetExcitationTime()const
{
    MrTime myTime;
    myTime = m_spExcitSTB->GetExcitationTime();
    myTime.RoundUp(_10us);
    return myTime;
}
bool EPIFIDKernel::HasEcho()const
{
    return true;
}
const vector<MrTime> EPIFIDKernel::GetEchoesTime()const
{
    vector<MrTime> paTE;
    paTE.push_back(m_TE);
    return paTE;
}
MrTime EPIFIDKernel::GetEchoSpacing()const
{
    return m_spAcquiSTB->GetEchoSpaceTime();
}

bool EPIFIDKernel::CalcRotMatrixFromLogToPhy( const Slice& rSlice, Umr::PatientPosition ePatientPosture,vector<double>* padRotMatLog2Phy )
{
    if (NULL == padRotMatLog2Phy)
    {
        SEQ_THROW("Null pointer In CalcRotMatrixFromLogToPhy");
    }
    vector <double> adOrientation;
    double dRotationAngle = 0.0;
    adOrientation = rSlice.GetOrientation();
    dRotationAngle =rSlice.GetRotDegree();
    CalcRotationMatrixFromLogicToPhysical( adOrientation, dRotationAngle,
                                           ePatientPosture, padRotMatLog2Phy);
    return true;
}

bool EPIFIDKernel::CalcMaxwellCoefficient( vector<double> adRotMatrixLog2Phy,double* pdMaxwellCoeff )
{
    if (NULL == pdMaxwellCoeff)
    {
        SEQ_THROW("Null pointer In CalcMaxwellCoefficient");
    }
    double dMaxwellValue;
    dMaxwellValue = LAMO_CONSTANT_1H * 1.0e6 / 2.0 / m_dB0_Tesla
        * pow(m_spAcquiSTB->GetGroAmpl_mT_m() * 1.0e-3,2)
        * m_dFov_Phase_mm * 1.0e-3
        * ( (m_spAcquiSTB->GetEchoSpaceTime())(_s,kDouble)
        - (m_spAcquiSTB->GetGroRampUpTime())(_s,kDouble) * 4.0/3.0);
    double dCoeffDirection;
    dCoeffDirection = 0.25*(pow(adRotMatrixLog2Phy[1]*adRotMatrixLog2Phy[6],2) + pow(adRotMatrixLog2Phy[4]*adRotMatrixLog2Phy[6],2))    // 0.25  = 1/4
        + (pow(adRotMatrixLog2Phy[0]*adRotMatrixLog2Phy[7],2) + pow(adRotMatrixLog2Phy[3]*adRotMatrixLog2Phy[7],2)
        - adRotMatrixLog2Phy[3]*adRotMatrixLog2Phy[4]*adRotMatrixLog2Phy[6]*adRotMatrixLog2Phy[7]
    - adRotMatrixLog2Phy[0]*adRotMatrixLog2Phy[1]*adRotMatrixLog2Phy[6]*adRotMatrixLog2Phy[7]);
    *pdMaxwellCoeff = dMaxwellValue * dCoeffDirection;
    return true;
}

double EPIFIDKernel::CalcRMS_SeqKernel()
{
    double dRMS_Kernel = 0.0;
    int64_t lTotalEchoNum = m_pKOrdering->GetPPARefLinesPE()*m_lPPASampleFactor + m_lStepsPPADummy;
    dRMS_Kernel = m_spExcitSTB->CalcRORMS_STB() + m_spPrePhaSTB->CalcRORMS_STB() 
                  + m_spAcquiSTB->CalcRORMS_STB()*(lTotalEchoNum + 3)
                  + m_spEPISpoilerSTB->CalcRORMS_STB();
    return dRMS_Kernel;
}

double EPIFIDKernel::CalcRMS_PPAKernel()
{
    double dRMS_Kernel = 0.0;
    int64_t lTotalEchoNum = m_pKOrdering->GetPPARefLinesPE()*m_lPPASampleFactor+ m_lStepsPPADummy;
    dRMS_Kernel = m_spExcitSTB->CalcRORMS_STB() +  m_spPrePhaSTB->CalcRORMS_STB() 
                  + m_spAcquiSTB->CalcRORMS_STB()*(lTotalEchoNum + 3)
                  + m_spEPISpoilerSTB->CalcRORMS_STB();
    return dRMS_Kernel;
}

EReason4PrepareFailed EPIFIDKernel::GetReasonKernelPrepareFailed()
{
    return m_lReasonPrepareFailed;
}

sGradArray EPIFIDKernel::GetKernelGradInfo4RMS() const
{
    return m_vecGradInfo;
}

void EPIFIDKernel::PrepareAllSTBGradInfo4RMS( vector<sGradArray>* pvec_GradInfo )
{
    vector<sGradArray> vec_GradInfo;
    vec_GradInfo.clear();
    // queue GradInfor for kernel
    int64_t lSTBStartTime_us = 0;
    m_spExcitSTB->SetSTBStartTime4RMS(lSTBStartTime_us);              // excitation is first
    vec_GradInfo.push_back(m_spExcitSTB->GetGradInfo4RMS());
    lSTBStartTime_us += m_spExcitSTB->GetDuration()(_us);
    m_spAcquiSTB->SetSTBStartTime4RMS(lSTBStartTime_us);              // three lines for phase corr
    vec_GradInfo.push_back(m_spAcquiSTB->GetGradInfo4RMS());
    lSTBStartTime_us += m_spAcquiSTB->GetDuration()(_us);
    m_spAcquiSTB->SetSTBStartTime4RMS(lSTBStartTime_us);
    vec_GradInfo.push_back(m_spAcquiSTB->GetGradInfo4RMS());
    lSTBStartTime_us += m_spAcquiSTB->GetDuration()(_us);
    m_spAcquiSTB->SetSTBStartTime4RMS(lSTBStartTime_us);
    vec_GradInfo.push_back(m_spAcquiSTB->GetGradInfo4RMS());
    lSTBStartTime_us += m_spAcquiSTB->GetDuration()(_us) + m_spFillTimeSTB->GetDuration()(_us);
    m_spPrePhaSTB->SetSTBStartTime4RMS(lSTBStartTime_us);             //   PrePhaseSTB is after FillTime
    vec_GradInfo.push_back(m_spPrePhaSTB->GetGradInfo4RMS());
    lSTBStartTime_us += m_spPrePhaSTB->GetDuration()(_us);            // AcqSTB is need more lines according to kordering
    int64_t lSegments = m_pKOrdering->GetSegments();
    for (int64_t ii = 0; ii < lSegments; ii++)
    {
        m_spAcquiSTB->SetSTBStartTime4RMS(lSTBStartTime_us);
        vec_GradInfo.push_back(m_spAcquiSTB->GetGradInfo4RMS());
        lSTBStartTime_us += m_spAcquiSTB->GetDuration()(_us);
    }
    m_spEPISpoilerSTB->SetSTBStartTime4RMS(lSTBStartTime_us);
    vec_GradInfo.push_back(m_spEPISpoilerSTB->GetGradInfo4RMS());
    (*pvec_GradInfo) = vec_GradInfo;
}

bool EPIFIDKernel::PrepareKernelGradInfo4RMS()
{
    vector<sGradArray> vec_GradInfo4RMS;
    vec_GradInfo4RMS.clear();
    PrepareAllSTBGradInfo4RMS(&vec_GradInfo4RMS);

    if (vec_GradInfo4RMS.empty())
    {
        return false;
    }
    m_vecGradInfo = vec_GradInfo4RMS[0];
    for (int64_t index = 1;index < vec_GradInfo4RMS.size();index++)
    {
        m_vecGradInfo += vec_GradInfo4RMS[index];
    }
    return true;
}
Umr_Sequence::eScanStatus EPIFIDKernel::GetScanStatus(const AcqIndex& rAcqIndex)
{
    if(KScanStatus_Prep == rAcqIndex.GetScanStatus() && 
        ((!m_bIsCorrDeltaPhase)||(m_bIsCorrDeltaPhase && 0 < rAcqIndex.GetUsertIndex(3))))
    {
        return kPPAMode;
    }
    else if (kScanStatus_Imaging == rAcqIndex.GetScanStatus())
    {
        return kImageMode;
    }
    else if(KScanStatus_Prep == rAcqIndex.GetScanStatus() && m_bIsCorrDeltaPhase && 0==rAcqIndex.GetUsertIndex(3))
    {
        return kCorrectMode;
    }
    else
    {
        return kErrorStatus; // 
    }
}
eSampleMode EPIFIDKernel::GetSampleMode(const AcqIndex& rAcqIndex)
{
    if (KScanStatus_Prep == rAcqIndex.GetScanStatus() && 
        ((!m_bIsCorrDeltaPhase) && 0 == rAcqIndex.GetUsertIndex(3))||(m_bIsCorrDeltaPhase && 1 == rAcqIndex.GetUsertIndex(3)))
    {
        return kPPA_OddMode;
    }
    else if (KScanStatus_Prep == rAcqIndex.GetScanStatus() && 
        ((!m_bIsCorrDeltaPhase) && 1 == rAcqIndex.GetUsertIndex(3))||(m_bIsCorrDeltaPhase && 2 == rAcqIndex.GetUsertIndex(3)))
    {
        return kPPA_EvenMode;
    }
    else
    {
        return kOtherSampleMode;
    }
}

void EPIFIDKernel::SetIsCorrDeltaPhase(bool bIsExtraCorr)
{
    m_bIsCorrDeltaPhase = bIsExtraCorr;
}

bool EPIFIDKernel::PNSCheck(EPNSExceedType& eResult)
{
    eResult = kNoExceed;

    double dPNSLimit = 1.25;
    double dGroAmp_log = m_spAcquiSTB->GetGroAmpl_mT_m();
    double dGpeAmp_log = m_spAcquiSTB->GetGpeBlipAmpl_mT_m();
    MrTime GroRT_log= m_spAcquiSTB->GetGroRampUpTime();
    MrTime GpeRT_log= m_spAcquiSTB->GetGpeRampUpTime();
    const Slice& CurrSlice = (*m_pSliceArray)[0];
    vector<double> adRotMatrix = CurrSlice.GetRotMatrix()->GetRotMatrix();

    double dPNS_Acqui = m_pEPIPNSCalc->CalcPNS_Acqui(dGroAmp_log,dGpeAmp_log,GroRT_log,GpeRT_log,adRotMatrix);

    if (dPNS_Acqui>dPNSLimit)
    {
        eResult = kAcquiExceed;
    }

    return true;
}


NAME_SPACE_END

//disable harmless compiler warning 
#if defined(WIN64)||defined(_MSC_VER)
#pragma warning(pop)
#endif

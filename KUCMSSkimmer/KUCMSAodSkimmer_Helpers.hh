//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "KUCMSAodSVSkimmer.hh"
#include "KUCMSHelperFunctions.hh"

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//	Helper functions that depended on varibles from GetEntry()
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getRhIdx( uInt rhDetID ){

  for( int idx = 0; idx < ECALRecHit_ID->size(); idx++ ){ if( rhDetID == (*ECALRecHit_ID)[idx] ) return idx; }
  //std::cout << " -- !! no rhDetID to (*rhID)[idx] match !! ---------------------- " << std::endl;
  return -1;

}//<<>>int KUCMSAodSkimmer::getRhIdx( int rhDetID )

uInt KUCMSAodSkimmer::getLeadRhID( std::vector<uInt> recHitIds ){

  uInt result;
  float enr(0.0);
  if( recHitIds.size() < 1 || ECALRecHit_ID->size() < 1 ){ return 0; }
  for( auto id : recHitIds ){
    auto rhenr = (*ECALRecHit_energy)[getRhIdx(id)];
    if( rhenr > enr ){ enr = rhenr; result = id; }
  }//<<>>for (const auto recHit : recHits )

  return result;

}//>>>>EcalRecHit KUCMSAodSkimmer::getLeadRh( rhGroup recHitsi

float KUCMSAodSkimmer::clstrR9( std::vector<uInt> recHitIds ){

  if( recHitIds.size() < 1 || ECALRecHit_ID->size() < 1 ){ return -2.0; }
  auto leadRhID = getLeadRhID( recHitIds );
  auto leadRhEn = (*ECALRecHit_energy)[getRhIdx(leadRhID)];
  float sumRhEn(0);
  if( recHitIds.size() < 1 ){ return -99.0; }
  for ( auto id : recHitIds ){ sumRhEn +=  (*ECALRecHit_energy)[getRhIdx(id)]; }
  return sumRhEn > 0 ? leadRhEn/sumRhEn  : - 1.0;

}//<<>>float KUCMSAodSkimmer::clstrR9( vector<uInt> recHitIds )

std::vector<float> KUCMSAodSkimmer::getRhGrpEnergies( std::vector<uInt> rechitids ){

  std::vector<float> result;
  if( rechitids.size() < 1 || ECALRecHit_ID->size() < 1 ){ result.push_back(-9.0); return result; }
  for ( auto id : rechitids ){ result.push_back((*ECALRecHit_energy)[getRhIdx(id)]); }
  return result;

};//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpEnergies( std::vector<uInt> rechitids )

std::vector<float> KUCMSAodSkimmer::getRhGrpTimes( std::vector<uInt> rechitids ){

  std::vector<float> result;
  if( rechitids.size() < 1 || ECALRecHit_ID->size() < 1 ){ result.push_back(-99.0); return result; }
  for ( auto id : rechitids ){
    auto rhtime = (*ECALRecHit_time)[getRhIdx(id)];
    auto rh0tof = (*ECALRecHit_0TOF)[getRhIdx(id)];
    auto rhpvtof = (*ECALRecHit_pvTOF)[getRhIdx(id)];
    result.push_back(rhtime+rh0tof-rhpvtof); 
  }//<<>>for ( auto id : recHitIds )
  return result;

};//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpTimes( std::vector<uInt> rechitids )

std::vector<float> KUCMSAodSkimmer::getRhGrpEigenFromAngles( std::vector<uInt> rechitids ){

  //bool verbose = true;
  bool verbose = false;

  if( DEBUG || verbose ) std::cout << " Starting getRhGrpEigen_sph  " << std::endl;

  std::vector<float> emptyReturn(9,-9);
  std::vector<float> egwts;
  std::vector<float> rhetas, rhphis;
  std::vector<float> logwtvec, tresvec;
  auto nRecHits = rechitids.size();
  if( nRecHits < 5 || ECALRecHit_ID->size() < 1 ){ if( verbose ) std::cout << " ----  rechit collection has too few rechits" << std::endl; return emptyReturn; }
  float sumRhEn(0);
  for ( auto id : rechitids ){ sumRhEn +=  (*ECALRecHit_energy)[getRhIdx(id)]; }
  if( verbose ) std::cout << " --- EigenAngles sumRhEn : " << sumRhEn << std::endl;
  if( sumRhEn <= 0 ){ if( verbose ) std::cout << " ----  rechit collection has no energy" << std::endl; return emptyReturn; }
  if( DEBUG ) std::cout << "1a, ";
  std::cout << " getRhGrpEigenFromAngles : nRechits: " << nRecHits << " -------------------------- " << std::endl;
  for( uInt it(0); it < nRecHits; it++ ){

    const auto rhIDX = getRhIdx(rechitids[it]);
    //auto idinfo = DetIDMap[rechitids[it]];
    auto idinfo = timeCali->getDetIdInfo(rechitids[it]);
    auto isEB = idinfo.ecal == 0;//ECAL::EB;
    if( rhIDX == -1 ){ if( verbose ) std::cout << " ---- Bad idx !!!!! -- In getRhGrpEigen ---- " << std::endl; return emptyReturn; }
    if( not isEB ){ if( verbose ) std::cout << " ---- rechit group has EE members " << idinfo.ecal << std::endl; return emptyReturn; }

    std::cout << " rhIDX: " << rhIDX << std::endl;
    const auto rhEtaPos = idinfo.i2;//recHitPos.ieta();
    std::cout << " rhEtaPos: " << rhEtaPos << std::endl;
    rhetas.push_back((rhEtaPos>0)?rhEtaPos+84.5:rhEtaPos+85.5);
    const auto rhPhiPos = idinfo.i1;//recHitPos.iphi();
    std::cout << " rhPhiPos: " << rhPhiPos << std::endl;
    rhphis.push_back(rhPhiPos-0.5);
    auto rhenergy = (*ECALRecHit_energy)[rhIDX];
    auto logwt = std::max(0.0, 4.2 + log(rhenergy/sumRhEn));// cut at rh energy < 1.5% of cluster
    logwtvec.push_back(logwt);

  }//<<>>for( uInt it(0); it < rechits.size(); it++ )

  std::vector<float> detas, dphis, angles, redlogwtvec;
  auto meta = mean( rhetas, logwtvec );
  auto mphi = meanIPhi( rhphis, logwtvec );
  for( uInt it(0); it < rhetas.size(); it++ ){

    float deta = rhetas[it]-meta;
    float dphi = dIPhi( rhphis[it], mphi );
    if( dphi < 0.5 && deta < 0.5 ) continue;
    detas.push_back(deta);
    dphis.push_back(dphi);
    float angle = getAngle( deta, dphi );
    angles.push_back(angle);
    redlogwtvec.push_back(logwtvec[it]);

  }//<<>>for( uInt it(0); it < etas.size(); it++ )

  auto eigens = getRhGrpEigen( angles, redlogwtvec );
  if( verbose ) std::cout << " --- EigenAngles VAlue : " << eigens[2] << std::endl;

  auto phiCorrFactor = 0.8;
  std::cout << " mean eta: " << meta << " mean phi: " << mphi << std::endl; 
  auto sxx = var( detas, meta, redlogwtvec );
  auto syy = var( dphis, mphi, redlogwtvec, accum(redlogwtvec)/phiCorrFactor );
  auto sxy = cvar( detas, meta, dphis, mphi, redlogwtvec, accum(redlogwtvec)/std::sqrt(phiCorrFactor) );
  std::cout << " sxx: " << sxx << " syy: " << syy << " sxy: " << sxy << std::endl;
  auto smaj = (sxx + syy + std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
  auto smin = (sxx + syy - std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
  auto sang = std::atan((sxx-syy+std::sqrt(sq2(syy-sxx)+4.*sq2(sxy)))/(2.*sxy));

  //eigens[0] //0 geoeigan x vec
  //eigens[1] //1 geoeigan y vec
  //eigens[2] //2 geoeigan mag vec
  eigens.push_back(smaj);//3
  eigens.push_back(smin);//4
  eigens.push_back(sang);//5

  if( DEBUG ) std::cout << " Done" << std::endl;;
  return eigens;

}//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpEigen( std::vector<uInt> rechitids )

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////!!!!!!!!!!!!!!!!!11  need rhPosX, rhPosY, & rhPosZ to be saved for following functions   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

std::vector<float> KUCMSAodSkimmer::getLeadTofRhTime( std::vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ ){

  std::vector<float> result;
  if( ( recHitIds.size() < 1 ) || ( ECALRecHit_ID->size() < 1 ) ){ result.push_back(-99); return result; }
  auto lrhid = getLeadRhID(recHitIds);
  auto lrhidx = getRhIdx(lrhid);
  auto X = (*ECALRecHit_rhx)[lrhidx];
  auto Y = (*ECALRecHit_rhy)[lrhidx];
  auto Z = (*ECALRecHit_rhz)[lrhidx];
  const auto d_rh = hypo( X, Y, Z);
  const auto d_pv = hypo( X-vtxX, Y-vtxY, Z-vtxZ);
  const auto tof = (d_rh-d_pv)/SOL;
  for( int idx = 0; idx < ECALRecHit_time->size(); idx++ ){result.push_back((*ECALRecHit_time)[idx]-tof);}
  return result;

}//>>>>vector<float> KUCMSAodSkimmer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )


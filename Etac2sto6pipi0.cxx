#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "VertexFit/Helix.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "McTruth/McParticle.h"
#include "VertexFit/IVertexDbSvc.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TRandom.h" 
using CLHEP::Hep3Vector;
#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep2Vector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "ParticleID/ParticleID.h"
#include "VertexFit/KalmanKinematicFit.h"
#include <vector>
#include "Etac2STo6PiP0Alg/Comb.h"
#include "Etac2STo6PiP0Alg/Etac2STo6PiP0.h"
typedef std::vector<int> Vint;
typedef std::vector<Comb> VComb;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
const double mpi = 0.13497;
const double xmass[5] = 
{
	0.000511, 0.105658, 0.139570, 0.493677, 0.938272
}
;
bool debug =  false;
const double mka = xmass[3];
int Ncut00,Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6,Ncut7,Ncut8,Ncut9,Ncut10;
////////////////////////////////////////////////////////////////////
Etac2STo6PiP0::Etac2STo6PiP0(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) 
{  
	declareProperty("CmsEnergy", m_ecms = 3.686);
	declareProperty("Vr0cut", m_vr0cut=1.0);
	declareProperty("Vz0cut", m_vz0cut=10.0);
	declareProperty("EnergyThreshold", m_energyThreshold=0.02);
	declareProperty("EnergyThreshold_b", m_energyThreshold_b=0.025);
	declareProperty("EnergyThreshold_e", m_energyThreshold_e=0.05);
	declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
	declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
	declareProperty("GammaAngCut", m_gammaAngCut=20.0);
	declareProperty("idNo",m_idNo = 100443);
	//Declare the properties  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Etac2STo6PiP0::initialize()
{
	MsgStream log(msgSvc(), name());

	log << MSG::INFO << "in initialize()" << endmsg;

	StatusCode status;

	NTuplePtr nt_mct_angular(ntupleSvc(), "FILE1/mct_angular");
	if (nt_mct_angular) {m_tuple_angular = nt_mct_angular;}
	else
	{
		m_tuple_angular = ntupleSvc() -> book("FILE1/mct_angular", CLID_ColumnWiseTuple,"ks N-tuple");

		if(m_tuple_angular)
		{
			status = m_tuple_angular -> addItem("egam1",  m_egam1);
			status = m_tuple_angular -> addItem("egam2",  m_egam2);
			status = m_tuple_angular -> addItem("egam",  m_egam);


			status = m_tuple_angular -> addItem("epi0_psip",  m_epi0_psip);
			status = m_tuple_angular -> addItem("epi0_hc",  m_epi0_hc);
			status = m_tuple_angular -> addItem("ppi0_psi",  m_ppi0_psip);
			status = m_tuple_angular -> addItem("ppi0_hc",  m_ppi0_hc);

		}
		else
		{
			log << MSG::ERROR << "Cannot book N-tuple: " << long(m_tuple_angular) << endmsg;

			return StatusCode::FAILURE;
		}
	}

	NTuplePtr nt2(ntupleSvc(), "FILE1/chi2omg");
	if ( nt2 ) m_tuple = nt2;
	else 
	{
		m_tuple = ntupleSvc()->book ("FILE1/chi2omg", CLID_ColumnWiseTuple,"chi2omg");
		if ( m_tuple )    
		{
			status = m_tuple->addItem("cc",        m_cc, 0,500);
			status = m_tuple->addItem("vx0",       m_cc,  m_vx0);
			status = m_tuple->addItem("vy0",       m_cc,  m_vy0);
			status = m_tuple->addItem("vz0",       m_cc,  m_vz0);
			status = m_tuple->addItem("vr0",       m_cc,  m_vr0);
			status = m_tuple->addItem("rvz0",      m_cc,  m_rvz0);
			status = m_tuple->addItem("rvxy0",     m_cc,  m_rvxy0);
			status = m_tuple->addItem("costheta",  m_cc,  m_costheta);
			status = m_tuple->addItem("rvphi0",    m_cc,  m_rvphi0);

			status = m_tuple->addItem("nn",        m_nn, 0,200);
			status = m_tuple->addItem("time",      m_nn, m_time);
			status = m_tuple->addItem("costhetaN", m_nn, m_costhetaN);
			status = m_tuple->addItem("eraw",      m_nn, m_eraw);
			status = m_tuple->addItem("dang",      m_nn, m_dang);
			status = m_tuple->addItem("dang2",     m_dang2);
			status = m_tuple->addItem("dang_pi",     m_dang_pi);
			status = m_tuple->addItem("dang_proton",     m_dang_proton);

			status = m_tuple->addItem("ntotCh",  m_ntotCh );
			status = m_tuple->addItem("ntotNe",  m_ntotNe );

			status = m_tuple->addItem("run",  m_run );
			status = m_tuple->addItem("rec",  m_rec );
			status = m_tuple->addItem("indexmc",   m_idxmc, 0, 100);
			status = m_tuple->addItem("pdgid",     m_idxmc, m_pdgid);
			status = m_tuple->addItem("trkIdx",     m_idxmc, m_trkIdx);
			status = m_tuple->addItem("motheridx", m_idxmc, m_motheridx);

			status = m_tuple->addItem("nGam",      m_nGam);
			status = m_tuple->addItem("mpipirec",   m_mpipirec);
			status = m_tuple->addItem("mggrec",   m_mggrec);






			status = m_tuple->addItem("eGam",     m_eGam);
			status = m_tuple->addItem("eGam1",     m_eGam1);
			status = m_tuple->addItem("eGam2",     m_eGam2);
			status = m_tuple->addItem("eraw1",     m_eraw1);
			status = m_tuple->addItem("eraw2",     m_eraw2);
			status = m_tuple->addItem("eoverp1",     m_eoverp1);
			status = m_tuple->addItem("eoverp2",     m_eoverp2);
			status = m_tuple->addItem("ptot1",     m_ptot1);
			status = m_tuple->addItem("ptot2",     m_ptot2);
			status = m_tuple->addItem("depth1",     m_depth1);
			status = m_tuple->addItem("depth2",     m_depth2);


			status = m_tuple->addItem("npp",     m_npp);//add for PID
			status = m_tuple->addItem("npm",     m_npm);//add for PID
			status = m_tuple->addItem("npim",     m_npim);//add for PID
			status = m_tuple->addItem("npip",     m_npip);//add for PID


			/////for bkg
			status = m_tuple->addItem("metac",    m_metac);

			status = m_tuple->addItem("m3pi_a",       m_m3pi_a);
			status = m_tuple->addItem("m3pi_b",       m_m3pi_b);
			status = m_tuple->addItem("m3pi_c",       m_m3pi_c);
			status = m_tuple->addItem("m3pi_d",       m_m3pi_d);
			status = m_tuple->addItem("m3pi_e",       m_m3pi_e);
			status = m_tuple->addItem("m3pi_f",       m_m3pi_f);
			status = m_tuple->addItem("m3pi_g",       m_m3pi_g);
			status = m_tuple->addItem("m3pi_h",       m_m3pi_h);
			status = m_tuple->addItem("m3pi_i",       m_m3pi_i);

			status = m_tuple->addItem("P3pi11",      4,       P_P3pi11);
			status = m_tuple->addItem("P3pi22",      4,       P_P3pi22);
			status = m_tuple->addItem("P3pi33",      4,       P_P3pi33);
			status = m_tuple->addItem("P3pi44",      4,       P_P3pi44);
			status = m_tuple->addItem("P3pi55",      4,       P_P3pi55);
			status = m_tuple->addItem("P3pi66",      4,       P_P3pi66);
			status = m_tuple->addItem("P3pi77",      4,       P_P3pi77);
			status = m_tuple->addItem("P3pi88",      4,       P_P3pi88);
			status = m_tuple->addItem("P3pi99",      4,       P_P3pi99);

			status = m_tuple->addItem("mgam3pi_a",       m_mgam3pi_a);
			status = m_tuple->addItem("mgam3pi_b",       m_mgam3pi_b);
			status = m_tuple->addItem("mgam3pi_c",       m_mgam3pi_c);
			status = m_tuple->addItem("mgam3pi_d",       m_mgam3pi_d);
			status = m_tuple->addItem("mgam3pi_e",       m_mgam3pi_e);
			status = m_tuple->addItem("mgam3pi_f",       m_mgam3pi_f);
			status = m_tuple->addItem("mgam3pi_g",       m_mgam3pi_g);
			status = m_tuple->addItem("mgam3pi_h",       m_mgam3pi_h);
			status = m_tuple->addItem("mgam3pi_i",       m_mgam3pi_i);

			status = m_tuple->addItem("chi2",      m_chi2);
			status = m_tuple->addItem("chi2_all",      m_chi2_all);

			status = m_tuple->addItem("P2pi11",      4,       P_P2pi11);
			status = m_tuple->addItem("P2pi22",      4,       P_P2pi22);
			status = m_tuple->addItem("P2pi33",      4,       P_P2pi33);
			status = m_tuple->addItem("P2pi44",      4,       P_P2pi44);
			status = m_tuple->addItem("P2pi55",      4,       P_P2pi55);
			status = m_tuple->addItem("P2pi66",      4,       P_P2pi66);
			status = m_tuple->addItem("P2pi77",      4,       P_P2pi77);
			status = m_tuple->addItem("P2pi88",      4,       P_P2pi88);
			status = m_tuple->addItem("P2pi99",      4,       P_P2pi99);
			status = m_tuple->addItem("m68g",    m_m68g);
			status = m_tuple->addItem("m78g",    m_m78g);
			status = m_tuple->addItem("m2pia",    m_m2pia);
			status = m_tuple->addItem("m2pib",    m_m2pib);
			status = m_tuple->addItem("m2pic",    m_m2pic);
			status = m_tuple->addItem("m2pid",    m_m2pid);
			status = m_tuple->addItem("m2pie",    m_m2pie);
			status = m_tuple->addItem("m2pif",    m_m2pif);
			status = m_tuple->addItem("m2pig",    m_m2pig);
			status = m_tuple->addItem("m2pih",    m_m2pih);
			status = m_tuple->addItem("m2pii",    m_m2pii);

			status = m_tuple->addItem("mgam2pia",    m_mgam2pia);
			status = m_tuple->addItem("mgam2pib",    m_mgam2pib);
			status = m_tuple->addItem("mgam2pic",    m_mgam2pic);
			status = m_tuple->addItem("mgam2pid",    m_mgam2pid);
			status = m_tuple->addItem("mgam2pie",    m_mgam2pie);
			status = m_tuple->addItem("mgam2pif",    m_mgam2pif);
			status = m_tuple->addItem("mgam2pig",    m_mgam2pig);
			status = m_tuple->addItem("mgam2pih",    m_mgam2pih);
			status = m_tuple->addItem("mgam2pii",    m_mgam2pii);


			status = m_tuple->addItem("m2pi_Ks_new",    m_m2pi_Ks_new);

			status = m_tuple->addItem("chi2_pi0",    m_chi2_pi0);
			status = m_tuple->addItem("chi2_pi02gam",    m_chi2_pi02gam);
			status = m_tuple->addItem("chi2_2pi0",    m_chi2_2pi0);
			status = m_tuple->addItem("chi22g",    m_chi22g);
			status = m_tuple->addItem("chi23g",    m_chi23g);
			status = m_tuple->addItem("chi24g",    m_chi24g);




			status = m_tuple->addItem("m2pirec_bkg",    m_m2pirec_bkg);
			status = m_tuple->addItem("m3gam_bkg",    m_m3gam_bkg);
			status = m_tuple->addItem("mpi0rec_bkg",    m_mpi0rec_bkg);


			status = m_tuple->addItem("P_pip1",      4,       m_P_pip1);
			status = m_tuple->addItem("P_pim1",      4,       m_P_pim1);
			status = m_tuple->addItem("P_pip2",      4,       m_P_pip2);
			status = m_tuple->addItem("P_pim2",      4,       m_P_pim2);
			status = m_tuple->addItem("P_pip3",      4,       m_P_pip3);
			status = m_tuple->addItem("P_pim3",      4,       m_P_pim3);

			status = m_tuple->addItem("P_gam6",      4,       m_P_gam6);
			status = m_tuple->addItem("P_gam7",      4,       m_P_gam7);
			status = m_tuple->addItem("P_pi0",      4,       m_P_pi0);
			status = m_tuple->addItem("P_gam8",      4,       m_P_gam8);




                        /////for 3C Method -0713
			status = m_tuple->addItem("metac_3C",    m_metac_3C);
			status = m_tuple->addItem("chi2_3C",      m_chi2_3C);



		}
		else   
		{ 
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple) <<
				endmsg;
			return StatusCode::FAILURE;
		}
	}
	//
	//--------end of book--------
	//

	//--------------------------------------------------------------
	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Etac2STo6PiP0::execute() 
{
	StatusCode sc = StatusCode::SUCCESS;

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	SmartDataPtr<EvtRecTrackCol> evtRecTrackCol( eventSvc(), "/Event/EvtRec/EvtRecTrackCol" );
	if ( !evtRecTrackCol  ) {
		log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
		return StatusCode::FAILURE;

	}

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	if (!eventHeader)
	{
		log << MSG::FATAL << "Could not find Event Header" << endreq;
		return StatusCode::FAILURE;
	}
	if(debug)cout<<__LINE__<<endl;
	m_run = eventHeader->runNumber();
	m_rec = eventHeader->eventNumber();
	int runNo = m_run;
	int event = m_rec;
	if(m_rec%10000==0)cout<<"Run   "<<m_run<<"     Event   "<<m_rec<<endl; 
	//MC information
	if(debug)cout<<__LINE__<<endl;
	if (eventHeader->runNumber()<0)
	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),
				"/Event/MC/McParticleCol");
		int m_numParticle = 0;
		if (!mcParticleCol)
		{
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}
		else
		{
			bool jpsiDecay = false;
			bool psipDecay = false;
			int rootIndex = -1;
			bool strange =false;
			Event::McParticleCol::iterator iter_mc;
			for(iter_mc=mcParticleCol->begin(); iter_mc != mcParticleCol->end();
					iter_mc++)
			{
				if((*iter_mc)->primaryParticle()&&(*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()==11) { strange =true;}
				if ((*iter_mc)->primaryParticle()) continue;
				if (!(*iter_mc)->decayFromGenerator()) continue;
				if ((*iter_mc)->particleProperty()==m_idNo)
				{
					psipDecay = true;
					rootIndex = (*iter_mc)->trackIndex();
				}
				if (!psipDecay) continue;

				int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
				if(strange && ((*iter_mc)->mother()).particleProperty()!=m_idNo) mcidx--;//added
				int pdgid = (*iter_mc)->particleProperty();
				m_pdgid[m_numParticle] = pdgid;
				m_trkIdx[m_numParticle]=rootIndex;
				m_motheridx[m_numParticle] = mcidx;
				m_numParticle += 1;


				///@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				double mcpx = (*iter_mc)->initialFourMomentum().px(); //最初px
				double mcpy = (*iter_mc)->initialFourMomentum().py(); //最初py
				double mcpz = (*iter_mc)->initialFourMomentum().pz(); //最初pz
				double mcen = (*iter_mc)->initialFourMomentum().e();  //最初能量
				HepLorentzVector p4mc;
				p4mc.setPx(mcpx);
				p4mc.setPy(mcpy);
				p4mc.setPz(mcpz);
				p4mc.setE(mcen);
				double mmm = p4mc.mag();
				double ppp = p4mc.rho();


				int mctrue_id = (*iter_mc)->particleProperty();
				int mctrue_mothPID = ((*iter_mc)->mother()).particleProperty();





				///@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@










			}
		}
		m_idxmc = m_numParticle;
		if(debug)cout<<__LINE__<<endl;

		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		Vp4 mctrue_pi0_psip,mctrue_pi0_hc;
		Vint index_pi0_psip,index_pi0_hc;
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			int mctrue_id = (*iter_mc)->particleProperty();
			int mctrue_mothPID = ((*iter_mc)->mother()).particleProperty();

			if((mctrue_id == 111)&&(mctrue_mothPID==100443)){
				mctrue_pi0_psip.push_back((*iter_mc)->initialFourMomentum());
				index_pi0_psip.push_back(((*iter_mc)->mother()).trackIndex());
			}

			if((mctrue_id == 111)&&(mctrue_mothPID==10443)){
				mctrue_pi0_hc.push_back((*iter_mc)->initialFourMomentum());
				index_pi0_hc.push_back(((*iter_mc)->mother()).trackIndex());
			}
		}
		if(mctrue_pi0_psip.size()==1&&mctrue_pi0_hc.size()==1)
		{
			HepLorentzVector ecms_mother(3.686*0.011,0.,0.,3.686);
			HepLorentzVector p4pi0_psip = mctrue_pi0_psip[0];
			HepLorentzVector p4pi0_hc = mctrue_pi0_hc[0];

			m_epi0_psip = p4pi0_psip.e();
			m_epi0_hc = p4pi0_hc.e();
			m_ppi0_psip =sqrt( p4pi0_psip.px()*p4pi0_psip.px() + p4pi0_psip.py()*p4pi0_psip.py() + p4pi0_psip.pz()*p4pi0_psip.pz()    );
			m_ppi0_hc = sqrt( p4pi0_hc.px()*p4pi0_hc.px() + p4pi0_hc.py()*p4pi0_hc.py() + p4pi0_hc.pz()*p4pi0_hc.pz()     ); 

			//0713			m_tuple_angular->write();
		}
	}

	if(debug)cout<<__LINE__<<endl;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),EventModel::EvtRec::EvtRecTrackCol);
	//	if(evtRecEvent->totalCharged()>10)return StatusCode::SUCCESS;
	//	if(evtRecEvent->totalNeutral()>20)return StatusCode::SUCCESS;
	m_ntotCh = evtRecEvent->totalCharged();
	m_ntotNe = evtRecEvent->totalNeutral();
	//HepLorentzVector psip(0.011*m_ecms,0.0,0.0,m_ecms);
	HepLorentzVector psip(0.011*m_ecms,0.0,0.0,m_ecms);

	Vint iGood,iChrgp,iChrgn;;
	iGood.clear();
	iChrgp.clear();
	iChrgn.clear();
	int nCharge = 0;
	Hep3Vector xorigin(0,0,0);
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	bool status = vtxsvc->isVertexValid();
	if(status)
	{
		double* dbv = vtxsvc->PrimaryVertex();     //vertex[0] = vx; vertex[1]= vy; vertex[2] = vz;
		double* dbvsigma = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}
	int cc=0;
	for(int i = 0;i < evtRecEvent->totalCharged(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *mdcTrk =(*itTrk)->mdcTrack();
		if (!(*itTrk)->isMdcKalTrackValid()) continue;//MdcKalTrk
		double x0 = mdcTrk->x();
		double y0 = mdcTrk->y();
		double z0 = mdcTrk->z();
		double r0 = mdcTrk->r();
		double phi0 = mdcTrk->helix(1);
		double ptrk = mdcTrk->p();
		double chi2 = mdcTrk->chi2();
		double xv=xorigin.x();
		double yv=xorigin.y();
		double zv=xorigin.z();
		double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
		m_vx0[cc] = x0;
		m_vy0[cc] = y0;
		m_vz0[cc] = z0;
		m_vr0[cc] = Rxy;
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=vecipa[1];       
		double Vctheta =cos(mdcTrk->theta());
		m_costheta[cc]=Vctheta;
		m_rvxy0[cc]=Rvxy0;
		m_rvz0[cc]=Rvz0;
		m_rvphi0[cc]=Rvphi0;
		cc++;
		if(fabs(Rvz0) >=m_vz0cut ) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;
		if(fabs(Vctheta)>0.93)continue;
		iGood.push_back(i);
		nCharge += mdcTrk->charge();		
		if(mdcTrk->charge()>0)
			iChrgp.push_back(i);
		else
			iChrgn.push_back(i);
	}
	m_cc=cc;
	Ncut0++; 
	if(debug)cout<<__LINE__<<endl;
	int ntot=iGood.size();
	if((ntot != 6 )||(nCharge != 0)) return sc;
	Ncut00++;
	int nChp=iChrgp.size();
	int nChn=iChrgn.size();
	if((nChp != 3)||(nChn != 3)) return sc;
	Ncut1++;

	//------------------------------------------------------------------------------------------
	//emc
	double eraw[2],eoverp[2],ptrk[2],px_[2],py_[2],pz_[2],ptot[2];
	int goodshower[2];
	for(int i=0;i<2;i++){ goodshower[i]=0;eraw[i]=eoverp[i]=ptrk[i]=0.; px_[i]=py_[i]=pz_[i]=ptot[i]=0.; }
	for(int i_emc = 0; i_emc < ntot; i_emc++)
	{
		EvtRecTrackIterator  itTrk = evtRecTrackCol->begin() + iGood[i_emc];
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		ptrk[i_emc] = mdcTrk->p();
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		eraw[i_emc] = emcTrk->energy();
		eoverp[i_emc] = eraw[i_emc]/ptrk[i_emc];

		px_[i_emc] = mdcTrk->px();
		py_[i_emc] = mdcTrk->py();
		pz_[i_emc] = mdcTrk->pz();
		ptot[i_emc] = sqrt(px_[i_emc]*px_[i_emc] + py_[i_emc]*py_[i_emc]+ pz_[i_emc]*pz_[i_emc]);

		goodshower[i_emc]++;
	}

	//muc
	double depth[2],hits[2],layers[2],brlastlay[2],eclastlay[2];

	for(int i=0;i<ntot;i++){
		layers[i] = -2;
		EvtRecTrackIterator  itTrk = evtRecTrackCol->begin() + iGood[i];
		if((*itTrk)->isMucTrackValid()){
			RecMucTrack *evt_mu = (*itTrk)->mucTrack();
			depth[i] = (evt_mu)->depth();
			hits[i] = (evt_mu)->numHits();
			layers[i] = (evt_mu)->numLayers();
			brlastlay[i]= (evt_mu)->brLastLayer();
			eclastlay[i] = (evt_mu)->ecLastLayer();
		}
	}

	m_eraw1 = eraw[0];
	m_eraw2 = eraw[1];

	m_eoverp1 = eoverp[0];
	m_eoverp2 = eoverp[1];

	m_ptot1 = ptot[0];
	m_ptot2 = ptot[1];

	m_depth1 = depth[0];
	m_depth2 = depth[1];

	//------------------------------------------------------------------------------------------

	Vint iGam;
	iGam.clear();
	int getTime=0;
	int nn = 0;
	double dang = 999.;
	for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

		//double dang = 200.;
		for(int j = 0; j < evtRecEvent->totalCharged(); j++)
		{  
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();
			double angd = extpos.angle(emcpos);
			if(angd < dang)dang = angd;
		}

		double eraw = emcTrk->energy();
		int time= emcTrk->time();
		double theta = emcTrk->theta();
		double costheta=cos(theta);
		dang = dang * 180 / (CLHEP::pi);
		m_costhetaN[nn]=costheta;
		m_dang[nn] = dang;
		m_eraw[nn]= eraw;
		m_time[nn]=time;
		nn++;
		if(fabs(costheta)<0.8)
		{
			//m_energyThreshold = 0.025;
			m_energyThreshold = m_energyThreshold_b;
		}

		else if(fabs(costheta)< 0.92&&fabs(costheta)>0.86)
		{
			//m_energyThreshold = 0.05;
			m_energyThreshold = m_energyThreshold_e;
		}
		else
		{
			continue;
		}

		if(time>14||time<0)continue;
		if(eraw < m_energyThreshold) continue;
		iGam.push_back(i);
	}
	m_nn=nn;
	m_dang2 = dang;
	int nGam=iGam.size();
	//	if((nGam<4||(nGam>12))) return sc;
	if((nGam<3)||(nGam>20)) return sc;
	//if(nGam<3) return sc;
	m_nGam = nGam;
	Ncut2++;


	if(debug)cout<<__LINE__<<endl;
	Vp4 pGam;
	pGam.clear();
	for (int i = 0; i < nGam; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + iGam[i];
		RecEmcShower* emcTrk = (*itTrk) -> emcShower();

		double eraw = emcTrk -> energy();
		double phi  = emcTrk -> phi();
		double the  = emcTrk -> theta();

		HepLorentzVector ptrk;
		ptrk.setPx(eraw * sin(the) * cos(phi));
		ptrk.setPy(eraw * sin(the) * sin(phi));
		ptrk.setPz(eraw * cos(the));
		ptrk.setE(eraw);
		//ptrk = ptrk.boost(-boov);
		pGam.push_back(ptrk);
	}
	if(debug)cout<<__LINE__<<endl;

	Vint ikap,ipip;
	ikap.clear();
	ipip.clear();
	Vp4 pkap,ppip;
	pkap.clear();
	ppip.clear();

	if(debug)cout<<__LINE__<<endl;
	//-----------------------Test fot PID ------------------------//


	double mpip = xmass[2],mpim = xmass[2],mpm = xmass[4],mpp = xmass[4];//pi+ pi- p- p质量
	//  Vp4 ppip,ppm,ppim,ppp;//动量
	Vp4 ppm,ppim,ppp;//动量
	//  Vint ipip,ipm,ipim,ipp;//容器 存储事例信息？
	Vint ipm,ipim,ipp;//容器 存储事例信息？
	ppip.clear();
	ppm.clear();
	ppim.clear();
	ppp.clear();
	ipip.clear();
	ipm.clear();
	ipim.clear();
	ipp.clear();
	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nChp; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgp[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()); // use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon() |pid->onlyProton());    // seperater Pion/Kaon
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		double m_prob_pi = pid->probPion();
		double m_prob_K = pid->probKaon();
		double m_prob_p = pid->probProton();
		//if((pid->probPion()<pid->probKaon())||(pid->probPion()<pid->probProton())) continue;
		if(!((m_prob_pi>m_prob_p)&&(m_prob_pi>m_prob_K)&&(m_prob_pi>0.001 )  )) continue;

		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su 
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);//PID can set to e //这里是pi+介子 pion
		ipip.push_back(iChrgp[i]);
		HepLorentzVector ptrkp;
		ptrkp.setPx(mdcKalTrk->px());
		ptrkp.setPy(mdcKalTrk->py());
		ptrkp.setPz(mdcKalTrk->pz());
		double p3 = ptrkp.mag();
		ptrkp.setE(sqrt(p3*p3+mpip*mpip));
		ppip.push_back(ptrkp);
	}
	for(int i = 0; i < nChn; i++)//chrn指  带负电径迹
	{


		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgn[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()); // use PID sub-system
		pid->identify(pid->onlyPion() | pid->onlyKaon() |pid->onlyProton());    // seperater Pion/Kaon
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		double m_prob_p = pid->probProton();
		double m_prob_pi = pid->probPion();
		double m_prob_K = pid->probKaon();
		if(!((m_prob_pi>m_prob_p)&&(m_prob_pi>m_prob_K)&&(m_prob_pi>0.001 )  )) continue;
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack su
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);//PID can set to e//pion指 pi介子
		ipim.push_back(iChrgn[i]);
		HepLorentzVector ptrkm;
		ptrkm.setPx(mdcKalTrk->px());
		ptrkm.setPy(mdcKalTrk->py());
		ptrkm.setPz(mdcKalTrk->pz());
		double p3 = ptrkm.mag();
		ptrkm.setE(sqrt(p3*p3+mpim*mpim));
		ppim.push_back(ptrkm);
	}
	int npim = ipim.size();
	int npip = ipip.size();

	m_npip=npip;//pi+ 数
	m_npim=npim;//pi- 数
	//     cout<<"m_npp "<<npp  <<endl;
	//   cout<<"m_npm "<<npm  <<endl;
	//   cout<<"m_npip "<<npip  <<endl;
	//   cout<<"m_npim "<<npim  <<endl;

	// if(npim<2||npp<1||npip<1) return sc;
	if(!(npim==3&& npip ==3)) return sc;//add  0721 
	Ncut3++;
	if(debug)cout<<__LINE__<<endl;

	//-----------------PID END--------------------------------------------------------------------------------//    
	//pi0 list

	KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
	Vp4 ppi0fit_list;
	Vdouble chi_pi0fit;
	Vdouble m_pi0fit;

	ppi0fit_list.clear();
	chi_pi0fit.clear();
	m_pi0fit.clear();

	Vint i1_pi0,i2_pi0;
	i1_pi0.clear();
	i2_pi0.clear();

	int igoodpi0=0;
	long pi0ntmother[100];
	for(int i=0;i<100;i++)
	{
		pi0ntmother[i]=0;

	}

	for(int i=0;i<nGam-1;i++)
	{
		for(int j=i+1;j<nGam;j++)
		{
			double m_pi0=(pGam[i]+pGam[j]).m();



			if(m_pi0<0.08||m_pi0>0.2) continue;
			HepLorentzVector ppi0fit2;
			RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
			RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();

			kmfit->init();
			kmfit->setBeamPosition(xorigin);//!!!!!!!!!//------------------------------------------------------
			kmfit->AddTrack(0, 0.0, g1Trk);
			kmfit->AddTrack(1, 0.0, g2Trk);
			kmfit->AddResonance(0, 0.1349766, 0, 1);
			kmfit->setChisqCut(9999,0.05);//-------------------------------------------------------------------
			kmfit->setIterNumber(3);//-------------------------------------------------------------------------
			bool   oksq    = kmfit->Fit();
			double pi0chi2 = kmfit->chisq();
			if ( pi0chi2 > 20 ) continue;
			ppi0fit2= kmfit->pfit(0)+kmfit->pfit(1);
			ppi0fit_list.push_back(ppi0fit2);
			chi_pi0fit.push_back(pi0chi2);
			i1_pi0.push_back(iGam[i]);
			i2_pi0.push_back(iGam[j]);

			m_pi0fit.push_back(m_pi0);

			igoodpi0++;
		}
	}

	int pi0_num=ppi0fit_list.size();
	if(debug)cout<<__LINE__<<endl;
	if( pi0_num<1 )return sc;



	if(debug)cout<<__LINE__<<endl;
	Ncut4++; 

	/////////////////////-------------ETA end----------------/////////
	//                Test vertex                        /////////////

	////////// add 0721//////////
	double chisum=999.;
	WTrackParameter wxip,wxim,wlamb1,wlamb2,wlambda1,wlambda,wtmp1,wtmp2;
	Vp4 hc_list;
	hc_list.clear();
	double chi = 9999.,dltmin=999.;
	int index_pim01 = -1, index_pim02 = -1, index_pp = -1;


	HepPoint3D vx(0.0, 0.0, 0.0);//初始化顶点
	HepSymMatrix Evx(3, 0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx*bx;
	Evx[1][1] = by*by;
	Evx[2][2] = bz*bz;

	VertexParameter vxpar;
	vxpar.setVx(vx);
	vxpar.setEvx(Evx);

	VertexFit* vtxfit = VertexFit::instance();
	double deltmin=999.;

	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);

	RecMdcKalTrack *pip1Trk = (*(evtRecTrkCol->begin()+ipip[0]))->mdcKalTrack();
	WTrackParameter wvpip1Trk = WTrackParameter(xmass[2], pip1Trk->getZHelix(), pip1Trk->getZError());

	RecMdcKalTrack *pim1Trk = (*(evtRecTrkCol->begin()+ipim[0]))->mdcKalTrack();
	WTrackParameter wvpim1Trk = WTrackParameter(xmass[2], pim1Trk->getZHelix(), pim1Trk->getZError());

	RecMdcKalTrack *pip2Trk = (*(evtRecTrkCol->begin()+ipip[1]))->mdcKalTrack();
	WTrackParameter wvpip2Trk = WTrackParameter(xmass[2], pip2Trk->getZHelix(), pip2Trk->getZError());

	RecMdcKalTrack *pim2Trk = (*(evtRecTrkCol->begin()+ipim[1]))->mdcKalTrack();
	WTrackParameter wvpim2Trk = WTrackParameter(xmass[2], pim2Trk->getZHelix(), pim2Trk->getZError());

	RecMdcKalTrack *pip3Trk = (*(evtRecTrkCol->begin()+ipip[2]))->mdcKalTrack();
	WTrackParameter wvpip3Trk = WTrackParameter(xmass[2], pip3Trk->getZHelix(), pip3Trk->getZError());

	RecMdcKalTrack *pim3Trk = (*(evtRecTrkCol->begin()+ipim[2]))->mdcKalTrack();
	WTrackParameter wvpim3Trk = WTrackParameter(xmass[2], pim3Trk->getZHelix(), pim3Trk->getZError());


	vtxfit->init();
	vtxfit->AddTrack(0, wvpip1Trk);
	vtxfit->AddTrack(1, wvpim1Trk);
	vtxfit->AddTrack(2, wvpip2Trk);
	vtxfit->AddTrack(3, wvpim2Trk);
	vtxfit->AddTrack(4, wvpip3Trk);
	vtxfit->AddTrack(5, wvpim3Trk);
	vtxfit->AddVertex(0, vxpar,0,1,2,3,4,5);
	if (!vtxfit->Fit(0)) return sc;
	vtxfit->Swim(0);

	Ncut5++; 
	WTrackParameter wpip1 = vtxfit->wtrk(0);
	WTrackParameter wpim1 = vtxfit->wtrk(1);

	WTrackParameter wpip2 = vtxfit->wtrk(2);
	WTrackParameter wpim2 = vtxfit->wtrk(3);

	WTrackParameter wpip3 = vtxfit->wtrk(4);
	WTrackParameter wpim3 = vtxfit->wtrk(5);

	HepPoint3D xorigin1 = vtxfit->vx(0);//add
	HepSymMatrix xem1 = vtxfit->Evx(0);//add

	if(debug)cout<<__LINE__<<endl;
	//Apply Kinematic fit  for Signal ------------
	double chisq_t = 999999.; int ig1 = -1; int ig2 = -1; int ig3 = -1; int ig4 = -1; double chisq; double pi0chi2;
	HepLorentzVector ppi01fit, petafit; int index_pi01 = -1, index_gam = -1;
	for(int i=0;i<pi0_num;i++)
	{
		for(int j=0 ;j< nGam;j++)
		{
			if(i1_pi0[i]==iGam[j]) continue;
			if(i2_pi0[i]==iGam[j]) continue;


			RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +i1_pi0[i])) -> emcShower();
			RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +i2_pi0[i])) -> emcShower();
			RecEmcShower *g3Trk = (*(evtRecTrkCol -> begin() +iGam[j])) -> emcShower();

			kmfit -> init();
			kmfit -> AddTrack(0, wpip1);
			kmfit -> AddTrack(1, wpim1);
			kmfit -> AddTrack(2, wpip2);
			kmfit -> AddTrack(3, wpim2);
			kmfit -> AddTrack(4, wpip3);
			kmfit -> AddTrack(5, wpim3);
			kmfit -> AddTrack(6, 0.0, g1Trk);
			kmfit -> AddTrack(7, 0.0, g2Trk);
			kmfit -> AddTrack(8, 0.0, g3Trk);  //
			kmfit -> AddResonance(0, 0.1349766, 6, 7);
			kmfit -> AddFourMomentum(1, psip);

			bool oksq = kmfit -> Fit();
			if(oksq){
				double chi2=kmfit->chisq();
				double chi2_all=chi2+chi_pi0fit[i];
				if (chi2_all < chisq_t)
				{
					chisq_t = chi2_all;
					chisq = chi2;

					index_pi01 = i;
					index_gam = j;
				}
			}
		}
	}
	if(chisq>200)return sc;
	if(index_pi01==-1||index_gam==-1 )return sc;
	if(debug)cout<<__LINE__<<endl;
	Ncut6++;
	RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +i1_pi0[index_pi01])) -> emcShower();
	RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +i2_pi0[index_pi01])) -> emcShower();
	RecEmcShower *g3Trk = (*(evtRecTrkCol -> begin() +iGam[index_gam])) -> emcShower();
	if(debug)cout<<__LINE__<<endl;
	kmfit -> init();
	kmfit -> AddTrack(0, wpip1);
	kmfit -> AddTrack(1, wpim1);
	kmfit -> AddTrack(2, wpip2);
	kmfit -> AddTrack(3, wpim2);
	kmfit -> AddTrack(4, wpip3);
	kmfit -> AddTrack(5, wpim3);
	kmfit -> AddTrack(6, 0.0, g1Trk);
	kmfit -> AddTrack(7, 0.0, g2Trk);
	kmfit -> AddTrack(8, 0.0, g3Trk);  //
	kmfit -> AddResonance(0, 0.1349766, 6, 7);
	kmfit -> AddFourMomentum(1, psip);
	if(!kmfit->Fit()) return sc;
	Ncut7++;

	if(debug)cout<<__LINE__<<endl;

	m_chi2 = chisq;
	m_chi2_all = chisq_t;

	//////=======**************=======
	////////////////   赋予四动量
	//////=======**************=======

	HepLorentzVector pip1 = kmfit->pfit(0);
	HepLorentzVector pim1 = kmfit->pfit(1);

	HepLorentzVector pip2 = kmfit->pfit(2);
	HepLorentzVector pim2 = kmfit->pfit(3);

	HepLorentzVector pip3 = kmfit->pfit(4);
	HepLorentzVector pim3 = kmfit->pfit(5);

	HepLorentzVector gam6 = kmfit->pfit(6);
	HepLorentzVector gam7 = kmfit->pfit(7);
	HepLorentzVector Ppi0 = kmfit->pfit(6)+kmfit->pfit(7);

	HepLorentzVector gam8 = kmfit->pfit(8);

	m_m68g=(gam6+gam8).m();
	m_m78g=(gam7+gam8).m();

	HepLorentzVector P2pi1 = ( pip1 + pim1   )  ;
	HepLorentzVector P2pi2 = ( pip2 + pim2   )  ;
	HepLorentzVector P2pi3 = ( pip3 + pim3   )  ;

	HepLorentzVector P2pi4 = ( pip1 + pim2   )  ;
	HepLorentzVector P2pi5 = ( pip2 + pim1   )  ;
	HepLorentzVector P2pi6 = ( pip2 + pim3   )  ;

	HepLorentzVector P2pi7 = ( pip3 + pim1   )  ;
	HepLorentzVector P2pi8 = ( pip1 + pim3   )  ;
	HepLorentzVector P2pi9 = ( pip3 + pim2   )  ;

	P_P2pi11[0]=P2pi1.px();P_P2pi11[1]=P2pi1.py(); P_P2pi11[2]=P2pi1.pz(); P_P2pi11[3]=P2pi1.e();
	P_P2pi22[0]=P2pi2.px();P_P2pi22[1]=P2pi2.py(); P_P2pi22[2]=P2pi2.pz(); P_P2pi22[3]=P2pi2.e();
	P_P2pi33[0]=P2pi3.px();P_P2pi33[1]=P2pi3.py(); P_P2pi33[2]=P2pi3.pz(); P_P2pi33[3]=P2pi3.e();

	P_P2pi44[0]=P2pi4.px();P_P2pi44[1]=P2pi4.py(); P_P2pi44[2]=P2pi4.pz(); P_P2pi44[3]=P2pi4.e();
	P_P2pi55[0]=P2pi5.px();P_P2pi55[1]=P2pi5.py(); P_P2pi55[2]=P2pi5.pz(); P_P2pi55[3]=P2pi5.e();
	P_P2pi66[0]=P2pi6.px();P_P2pi66[1]=P2pi6.py(); P_P2pi66[2]=P2pi6.pz(); P_P2pi66[3]=P2pi6.e();

	P_P2pi77[0]=P2pi7.px();P_P2pi77[1]=P2pi7.py(); P_P2pi77[2]=P2pi7.pz(); P_P2pi77[3]=P2pi7.e();
	P_P2pi88[0]=P2pi8.px();P_P2pi88[1]=P2pi8.py(); P_P2pi88[2]=P2pi8.pz(); P_P2pi88[3]=P2pi8.e();
	P_P2pi99[0]=P2pi9.px();P_P2pi99[1]=P2pi9.py(); P_P2pi99[2]=P2pi9.pz(); P_P2pi99[3]=P2pi9.e();

	m_m2pia=P2pi1.m();
	m_m2pib=P2pi2.m();
	m_m2pic=P2pi3.m();
	m_m2pid=P2pi4.m();
	m_m2pie=P2pi5.m();
	m_m2pif=P2pi6.m();
	m_m2pig=P2pi7.m();
	m_m2pih=P2pi8.m();
	m_m2pii=P2pi9.m();

	m_mgam2pia=(P2pi1+gam8).m();
	m_mgam2pib=(P2pi2+gam8).m();
	m_mgam2pic=(P2pi3+gam8).m();
	m_mgam2pid=(P2pi4+gam8).m();
	m_mgam2pie=(P2pi5+gam8).m();
	m_mgam2pif=(P2pi6+gam8).m();
	m_mgam2pig=(P2pi7+gam8).m();
	m_mgam2pih=(P2pi8+gam8).m();
	m_mgam2pii=(P2pi9+gam8).m();
	//////=======**************=======
	////////////////   赋予四动量END
	//////=======**************=======
	//////////////2021  0928


	double tmKs1 = fabs( P2pi1.m() - 0.497611 );
	double tmKs2 = fabs( P2pi2.m() - 0.497611 );
	double tmKs3 = fabs( P2pi3.m() - 0.497611 );
	double tmKs4 = fabs( P2pi4.m() - 0.497611 );
	double tmKs5 = fabs( P2pi5.m() - 0.497611 );
	double tmKs6 = fabs( P2pi6.m() - 0.497611 );
	double tmKs7 = fabs( P2pi7.m() - 0.497611 );
	double tmKs8 = fabs( P2pi8.m() - 0.497611 );
	double tmKs9 = fabs( P2pi9.m() - 0.497611 );
	if( tmKs1 < tmKs2 && tmKs1 <tmKs3 && tmKs1 < tmKs4 &&tmKs1 < tmKs5 &&tmKs1 < tmKs6 &&tmKs1 < tmKs7 &&tmKs1 < tmKs8 &&tmKs1 < tmKs9 )
	{
		m_m2pi_Ks_new = P2pi1.m();
	}
	else if( tmKs2 < tmKs1 && tmKs2 <tmKs3 && tmKs2 < tmKs4 &&tmKs2 < tmKs5 &&tmKs2 < tmKs6 &&tmKs2 < tmKs7 &&tmKs2 < tmKs8 &&tmKs2 < tmKs9)
	{
		m_m2pi_Ks_new =P2pi2.m();
	}
	else if( tmKs3 < tmKs1 && tmKs3 <tmKs2 && tmKs3 < tmKs4 &&tmKs3 < tmKs5 &&tmKs3 < tmKs6 &&tmKs3 < tmKs7 &&tmKs3 < tmKs8 &&tmKs3 < tmKs9)
	{
		m_m2pi_Ks_new =P2pi3.m();
	}
	else if( tmKs4 < tmKs1 && tmKs4 <tmKs2 && tmKs4 < tmKs3 &&tmKs4 < tmKs5 &&tmKs4 < tmKs6 &&tmKs4 < tmKs7 &&tmKs4 < tmKs8 &&tmKs4 < tmKs9)
	{
		m_m2pi_Ks_new =P2pi4.m();
	}
	else if( tmKs5 < tmKs1 && tmKs5 <tmKs2 && tmKs5 < tmKs3 &&tmKs5 < tmKs4 &&tmKs5 < tmKs6 &&tmKs5 < tmKs7 &&tmKs5 < tmKs8 &&tmKs5 < tmKs9)
	{
		m_m2pi_Ks_new =P2pi5.m();
	}
	else if( tmKs6 < tmKs1 && tmKs6 <tmKs2 && tmKs6 < tmKs3 &&tmKs6 < tmKs4 &&tmKs6 < tmKs5 &&tmKs6 < tmKs7 &&tmKs6 < tmKs8 &&tmKs6 < tmKs9)
	{
		m_m2pi_Ks_new =P2pi6.m();
	}
	else if( tmKs7 < tmKs1 && tmKs7 <tmKs2 && tmKs7 < tmKs3 &&tmKs7 < tmKs4 &&tmKs7 < tmKs5 &&tmKs7 < tmKs6 &&tmKs7 < tmKs8 &&tmKs7 < tmKs9)
	{
		m_m2pi_Ks_new =P2pi7.m();
	}
	else if( tmKs8 < tmKs1 && tmKs8 <tmKs2 && tmKs8 < tmKs3 &&tmKs8 < tmKs4 &&tmKs8 < tmKs5 &&tmKs8 < tmKs6 &&tmKs8 < tmKs7 &&tmKs8 < tmKs9)
	{
		m_m2pi_Ks_new =P2pi8.m();
	}
	else if( tmKs9 < tmKs1 && tmKs9 <tmKs2 && tmKs9 < tmKs3 &&tmKs9 < tmKs4 &&tmKs9 < tmKs5 &&tmKs9 < tmKs6 &&tmKs9 < tmKs7 &&tmKs9 < tmKs8)
	{
		m_m2pi_Ks_new = P2pi9.m();
	}


	if(debug)cout<<__LINE__<<endl;
	if(debug)cout<<__LINE__<<endl;

	m_metac =( kmfit->pfit(0)+kmfit->pfit(1)+ kmfit->pfit(2)+kmfit->pfit(3)+kmfit->pfit(4)+kmfit->pfit(5)+ kmfit->pfit(6)+kmfit->pfit(7)).m() ;

	if(debug) cout<<__LINE__<<endl;






	HepLorentzVector P3pi1 = ( P2pi1  + Ppi0   )  ;
	HepLorentzVector P3pi2 = ( P2pi2  + Ppi0   )  ;
	HepLorentzVector P3pi3 = ( P2pi3  + Ppi0   )  ;
	HepLorentzVector P3pi4 = ( P2pi4  + Ppi0   )  ;
	HepLorentzVector P3pi5 = ( P2pi5  + Ppi0   )  ;
	HepLorentzVector P3pi6 = ( P2pi6  + Ppi0   )  ;
	HepLorentzVector P3pi7 = ( P2pi7  + Ppi0   )  ;
	HepLorentzVector P3pi8 = ( P2pi8  + Ppi0   )  ;
	HepLorentzVector P3pi9 = ( P2pi9  + Ppi0   )  ;
	m_m3pi_a=P3pi1.m();
	m_m3pi_b=P3pi2.m();
	m_m3pi_c=P3pi3.m();
	m_m3pi_d=P3pi4.m();
	m_m3pi_e=P3pi5.m();
	m_m3pi_f=P3pi6.m();
	m_m3pi_g=P3pi7.m();
	m_m3pi_h=P3pi8.m();
	m_m3pi_i=P3pi9.m();


	P_P3pi11[0]=P3pi1.px();P_P3pi11[1]=P3pi1.py(); P_P3pi11[2]=P3pi1.pz();P_P3pi11[3]=P3pi1.e();
	P_P3pi22[0]=P3pi2.px();P_P3pi22[1]=P3pi2.py(); P_P3pi22[2]=P3pi2.pz();P_P3pi22[3]=P3pi2.e();
	P_P3pi33[0]=P3pi3.px();P_P3pi33[1]=P3pi3.py(); P_P3pi33[2]=P3pi3.pz();P_P3pi33[3]=P3pi3.e();

	P_P3pi44[0]=P3pi4.px();P_P3pi44[1]=P3pi4.py(); P_P3pi44[2]=P3pi4.pz();P_P3pi44[3]=P3pi4.e();
	P_P3pi55[0]=P3pi5.px();P_P3pi55[1]=P3pi5.py(); P_P3pi55[2]=P3pi5.pz();P_P3pi55[3]=P3pi5.e();
	P_P3pi66[0]=P3pi6.px();P_P3pi66[1]=P3pi6.py(); P_P3pi66[2]=P3pi6.pz();P_P3pi66[3]=P3pi6.e();

	P_P3pi77[0]=P3pi7.px();P_P3pi77[1]=P3pi7.py(); P_P3pi77[2]=P3pi7.pz();P_P3pi77[3]=P3pi7.e();
	P_P3pi88[0]=P3pi8.px();P_P3pi88[1]=P3pi8.py(); P_P3pi88[2]=P3pi8.pz();P_P3pi88[3]=P3pi8.e();
	P_P3pi99[0]=P3pi9.px();P_P3pi99[1]=P3pi9.py(); P_P3pi99[2]=P3pi9.pz();P_P3pi99[3]=P3pi9.e();


	m_mgam3pi_a=(P3pi1+gam8 ).m();
	m_mgam3pi_b=(P3pi2+gam8 ).m();
	m_mgam3pi_c=(P3pi3+gam8 ).m();
	m_mgam3pi_d=(P3pi4+gam8 ).m();
	m_mgam3pi_e=(P3pi5+gam8 ).m();
	m_mgam3pi_f=(P3pi6+gam8 ).m();
	m_mgam3pi_g=(P3pi7+gam8 ).m();
	m_mgam3pi_h=(P3pi8+gam8 ).m();
	m_mgam3pi_i=(P3pi9+gam8 ).m();




	if(debug) cout<<__LINE__<<endl;

	if(debug) cout<<__LINE__<<endl;

	//  *** @@@@@@@@@  Selection for 2pirec === jpsi @@@@@@@@@  **** //
	//  *** @@@@@@@@@  挑选: pipi反冲(jpsi) @@@@@@@@@  **** //
	HepLorentzVector pipirec1 = (psip- pip1 -pim1   )  ;
	HepLorentzVector pipirec2 = (psip- pip2 -pim2   )  ;
	HepLorentzVector pipirec3 = (psip- pip3 -pim3   )  ;


	HepLorentzVector pipirec4 = (psip- pip1 -pim2   )  ;
	HepLorentzVector pipirec5 = (psip- pip2 -pim1   )  ;
	HepLorentzVector pipirec6 = (psip- pip2 -pim3   )  ;

	HepLorentzVector pipirec7 = (psip- pip3 -pim1   )  ;
	HepLorentzVector pipirec8 = (psip- pip1 -pim3   )  ;
	HepLorentzVector pipirec9 = (psip- pip3 -pim2   )  ;


	double tem1 = fabs(pipirec1.m() - 3.0969);
	double tem2 = fabs(pipirec2.m() - 3.0969);
	double tem3 = fabs(pipirec3.m() - 3.0969);
	double tem4 = fabs(pipirec4.m() - 3.0969);
	double tem5 = fabs(pipirec5.m() - 3.0969);
	double tem6 = fabs(pipirec6.m() - 3.0969);
	double tem7 = fabs(pipirec7.m() - 3.0969);
	double tem8 = fabs(pipirec8.m() - 3.0969);
	double tem9 = fabs(pipirec9.m() - 3.0969);
	if(tem1 < tem2 && tem1 <tem3 && tem1 < tem4 &&tem1 < tem5 &&tem1 < tem6 &&tem1 < tem7 &&tem1 < tem8 &&tem1 < tem9 )
	{
		m_m2pirec_bkg = pipirec1.m();
	}
	else if(tem2 < tem1 && tem2 <tem3 && tem2 < tem4 &&tem2 < tem5 &&tem2 < tem6 &&tem2 < tem7 &&tem2 < tem8 &&tem2 < tem9)
	{
		m_m2pirec_bkg = pipirec2.m();
	}
	else if(tem3 < tem1 && tem3 <tem2 && tem3 < tem4 &&tem3 < tem5 &&tem3 < tem6 &&tem3 < tem7 &&tem3 < tem8 &&tem3 < tem9)
	{
		m_m2pirec_bkg = pipirec3.m();
	}
	else if(tem4 < tem1 && tem4 <tem2 && tem4 < tem3 &&tem4 < tem5 &&tem4 < tem6 &&tem4 < tem7 &&tem4 < tem8 &&tem4 < tem9)
	{
		m_m2pirec_bkg = pipirec4.m();
	}
	else if(tem5 < tem1 && tem5 <tem2 && tem5 < tem3 &&tem5 < tem4 &&tem5 < tem6 &&tem5 < tem7 &&tem5 < tem8 &&tem5 < tem9)
	{
		m_m2pirec_bkg = pipirec5.m();
	}
	else if(tem6 < tem1 && tem6 <tem2 && tem6 < tem3 &&tem6 < tem4 &&tem6 < tem5 &&tem6 < tem7 &&tem6 < tem8 &&tem6 < tem9)
	{
		m_m2pirec_bkg = pipirec6.m();
	}
	else if(tem7 < tem1 && tem7 <tem2 && tem7 < tem3 &&tem7 < tem4 &&tem7 < tem5 &&tem7 < tem6 &&tem7 < tem8 &&tem7 < tem9)
	{
		m_m2pirec_bkg = pipirec7.m();
	}
	else if(tem8 < tem1 && tem8 <tem2 && tem8 < tem3 &&tem8 < tem4 &&tem8 < tem5 &&tem8 < tem6 &&tem8 < tem7 &&tem8 < tem9)
	{
		m_m2pirec_bkg = pipirec8.m();
	}
	else if(tem9 < tem1 && tem9 <tem2 && tem9 < tem3 &&tem9 < tem4 &&tem9 < tem5 &&tem9 < tem6 &&tem9 < tem7 &&tem9 < tem8)
	{
		m_m2pirec_bkg = pipirec9.m();
	}




	m_mpi0rec_bkg =(psip - kmfit ->pfit(6) - kmfit->pfit(7)  ).m();
	m_m3gam_bkg =( kmfit ->pfit(6) +  kmfit->pfit(7) + kmfit ->pfit(8) ).m();


	m_P_pip1[0]=pip1.px();
	m_P_pip1[1]=pip1.py();
	m_P_pip1[2]=pip1.pz();
	m_P_pip1[3]=pip1.e();

	m_P_pim1[0]=pim1.px();
	m_P_pim1[1]=pim1.py();
	m_P_pim1[2]=pim1.pz();
	m_P_pim1[3]=pim1.e();

	m_P_pip2[0]=pip2.px();
	m_P_pip2[1]=pip2.py();
	m_P_pip2[2]=pip2.pz();
	m_P_pip2[3]=pip2.e();

	m_P_pim2[0]=pim2.px();
	m_P_pim2[1]=pim2.py();
	m_P_pim2[2]=pim2.pz();
	m_P_pim2[3]=pim2.e();

	m_P_pip3[0]=pip3.px();
	m_P_pip3[1]=pip3.py();
	m_P_pip3[2]=pip3.pz();
	m_P_pip3[3]=pip3.e();

	m_P_pim3[0]=pim3.px();
	m_P_pim3[1]=pim3.py();
	m_P_pim3[2]=pim3.pz();
	m_P_pim3[3]=pim3.e();


	m_P_gam6[0]=gam6.px();
	m_P_gam6[1]=gam6.py();
	m_P_gam6[2]=gam6.pz();
	m_P_gam6[3]=gam6.e();

	m_P_gam7[0]=gam7.px();
	m_P_gam7[1]=gam7.py();
	m_P_gam7[2]=gam7.pz();
	m_P_gam7[3]=gam7.e();

	m_P_gam8[0]=gam8.px();
	m_P_gam8[1]=gam8.py();
	m_P_gam8[2]=gam8.pz();
	m_P_gam8[3]=gam8.e();

	m_P_pi0[0]=Ppi0.px();
	m_P_pi0[1]=Ppi0.py();
	m_P_pi0[2]=Ppi0.pz();
	m_P_pi0[3]=Ppi0.e();




	////////////////////////////////////////////////////////////////////////////////////////
	//**
	//*******************************************************************************************
	//*********************   3-C Method  *******************************
	//*******************************************************************************************
	//**0713
	KalmanKinematicFit *kmfit3C = KalmanKinematicFit::instance() ;

	RecEmcShower *Gam1Trk = (*(evtRecTrkCol -> begin() +i1_pi0[index_pi01])) -> emcShower();
	RecEmcShower *Gam2Trk = (*(evtRecTrkCol -> begin() +i2_pi0[index_pi01])) -> emcShower();
	RecEmcShower *Gam3Trk = (*(evtRecTrkCol -> begin() +iGam[index_gam])) -> emcShower();
	if(debug)cout<<__LINE__<<endl;
	kmfit3C -> init();
	kmfit3C -> AddTrack(0, wpip1);
	kmfit3C -> AddTrack(1, wpim1);
	kmfit3C -> AddTrack(2, wpip2);
	kmfit3C -> AddTrack(3, wpim2);
	kmfit3C -> AddTrack(4, wpip3);
	kmfit3C -> AddTrack(5, wpim3);
	kmfit3C -> AddTrack(6, 0.0, Gam1Trk);
	kmfit3C -> AddTrack(7, 0.0, Gam2Trk);
	kmfit3C -> AddResonance(0, 0.1349766, 6, 7);
        kmfit3C -> AddMissTrack(8, 0.0, Gam3Trk);  //
	kmfit3C -> AddFourMomentum(1, psip);
	bool oksq3C = kmfit3C->Fit();
	if (oksq3C && kmfit3C->chisq()>0) {


		if(debug)cout<<__LINE__<<endl;

		m_chi2_3C = kmfit3C->chisq();


		//////=======**************=======
		////////////////   赋予四动量
		//////=======**************=======

		m_metac_3C =( kmfit3C->pfit(0)+kmfit3C->pfit(1)+ kmfit3C->pfit(2)+kmfit3C->pfit(3)+kmfit3C->pfit(4)+kmfit3C->pfit(5)+ kmfit3C->pfit(6)+kmfit3C->pfit(7)).m() ;

	}

	//*******************************************************************************************
	//********************* END   3-C Method  *******************************
	//*******************************************************************************************
	//**
	////////////////////////////////////////////////////////////////////////////////////////










	/////=========**********  for possible BKG. ****===========/ ////
	double chisq_pi0=9999;
	for(int i=0;i<pi0_num;i++)
	{

		RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +i1_pi0[i])) -> emcShower();
		RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +i2_pi0[i])) -> emcShower();

		kmfit -> init();
		kmfit -> AddTrack(0, wpip1);
		kmfit -> AddTrack(1, wpim1);
		kmfit -> AddTrack(2, wpip2);
		kmfit -> AddTrack(3, wpim2);
		kmfit -> AddTrack(4, wpip3);
		kmfit -> AddTrack(5, wpim3);
		kmfit -> AddTrack(6, 0.0, g1Trk);
		kmfit -> AddTrack(7, 0.0, g2Trk);
		kmfit -> AddResonance(0, 0.1349766, 6, 7);
		kmfit -> AddFourMomentum(1, psip);
		bool oksq = kmfit->Fit();
		if(debug) cout<<__LINE__<<endl;

		if( oksq){
			double chi2 = kmfit->chisq();
			if(chi2 < chisq_pi0 ){

				chisq_pi0 = chi2;
			}
		}
		if(debug) cout<<__LINE__<<endl;
	}
	m_chi2_pi0 =chisq_pi0;
	if(debug) cout<<__LINE__<<endl;

	double chisq_pi02gam=9999;
	for(int i=0;i<pi0_num;i++)
	{
		for(int j=0 ;j< nGam-1;j++)
		{
			for(int k=j+1 ;k< nGam;k++)
			{
				if(i1_pi0[i]==iGam[j]) continue;
				if(i2_pi0[i]==iGam[j]) continue;

				if(i1_pi0[i]==iGam[k]) continue;
				if(i2_pi0[i]==iGam[k]) continue;

				RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +i1_pi0[i])) -> emcShower();
				RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +i2_pi0[i])) -> emcShower();
				RecEmcShower *g3Trk = (*(evtRecTrkCol -> begin() +iGam[j])) -> emcShower();
				RecEmcShower *g4Trk = (*(evtRecTrkCol -> begin() +iGam[k])) -> emcShower();

				kmfit -> init();
				kmfit -> AddTrack(0, wpip1);
				kmfit -> AddTrack(1, wpim1);
				kmfit -> AddTrack(2, wpip2);
				kmfit -> AddTrack(3, wpim2);
				kmfit -> AddTrack(4, wpip3);
				kmfit -> AddTrack(5, wpim3);
				kmfit -> AddTrack(6, 0.0, g1Trk);
				kmfit -> AddTrack(7, 0.0, g2Trk);
				kmfit -> AddTrack(8, 0.0, g3Trk);  
				kmfit -> AddTrack(9, 0.0, g4Trk);  
				kmfit -> AddResonance(0, 0.1349766, 6, 7);
				kmfit -> AddFourMomentum(1, psip);
				bool oksq = kmfit->Fit();
				if(debug) cout<<__LINE__<<endl;

				if( oksq){
					double chi2 = kmfit->chisq();
					if(chi2 < chisq_pi02gam ){

						chisq_pi02gam = chi2;
					}
				}
				if(debug) cout<<__LINE__<<endl;
			}
		}
	}
	m_chi2_pi02gam =chisq_pi02gam;
	if(debug) cout<<__LINE__<<endl;

	double chisq_2pi0=9999;
	for(int i=0;i< pi0_num-1;i++)
	{
		for(int j=i+1 ;j< pi0_num;j++)
		{

			if(i1_pi0[i]==i1_pi0[j]) continue;
			if(i2_pi0[i]==i1_pi0[j]) continue;

			if(i1_pi0[i]==i2_pi0[j]) continue;
			if(i2_pi0[i]==i2_pi0[j]) continue;

			RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +i1_pi0[i])) -> emcShower();
			RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +i2_pi0[i])) -> emcShower();
			RecEmcShower *g3Trk = (*(evtRecTrkCol -> begin() +i1_pi0[j])) -> emcShower();
			RecEmcShower *g4Trk = (*(evtRecTrkCol -> begin() +i2_pi0[j])) -> emcShower();

			kmfit->init();
			kmfit -> AddTrack(0, wpip1);
			kmfit -> AddTrack(1, wpim1);
			kmfit -> AddTrack(2, wpip2);
			kmfit -> AddTrack(3, wpim2);
			kmfit -> AddTrack(4, wpip3);
			kmfit -> AddTrack(5, wpim3);
			kmfit -> AddTrack(6, 0.0, g1Trk);
			kmfit -> AddTrack(7, 0.0, g2Trk);
			kmfit -> AddTrack(8, 0.0, g3Trk);  
			kmfit -> AddTrack(9, 0.0, g4Trk);  
			kmfit -> AddResonance(0, 0.135, 6, 7);
			kmfit -> AddResonance(1, 0.135, 8, 9);
			kmfit -> AddFourMomentum(2, psip);
			bool oksq = kmfit->Fit();

			if(debug) cout<<__LINE__<<endl;

			if( oksq){
				double chi2 = kmfit->chisq();
				if(chi2 < chisq_2pi0 ){

					chisq_2pi0 = chi2;
				}
			}
			if(debug) cout<<__LINE__<<endl;
		}
	}
	m_chi2_2pi0 =chisq_2pi0;
	if(debug) cout<<__LINE__<<endl;





	if(debug)cout<<__LINE__<<endl;

	double chisq_t2g=9999;
	double chisq_t3g=9999;
	double chisq_t4g=9999;

	for(int k=0;k<nGam-1;k++)
	{
		for(int m=k+1;m<nGam;m++)
		{

			RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +iGam[k])) -> emcShower();
			RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +iGam[m])) -> emcShower();
			kmfit -> init();
			kmfit -> AddTrack(0, wpip1);
			kmfit -> AddTrack(1, wpim1);
			kmfit -> AddTrack(2, wpip2);
			kmfit -> AddTrack(3, wpim2);
			kmfit -> AddTrack(4, wpip3);
			kmfit -> AddTrack(5, wpim3);
			kmfit -> AddTrack(6, 0.0, g1Trk);
			kmfit -> AddTrack(7, 0.0, g2Trk);
			kmfit -> AddFourMomentum(0, psip);

			bool oksq = kmfit -> Fit();
			if(oksq){
				double chi2=kmfit->chisq();
				if (chi2 < chisq_t2g)
				{
					chisq_t2g = chi2;
				}
			}
		}
	}
	m_chi22g=chisq_t2g;
	if(debug)cout<<__LINE__<<endl;



	for(int j=0;j<nGam-2;j++)
	{
		for(int k=j+1;k<nGam-1;k++)
		{
			for(int m=k+1;m<nGam;m++)
			{

				RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +iGam[j])) -> emcShower();
				RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +iGam[k])) -> emcShower();
				RecEmcShower *g3Trk = (*(evtRecTrkCol -> begin() +iGam[m])) -> emcShower();
				kmfit -> init();
				kmfit -> AddTrack(0, wpip1);
				kmfit -> AddTrack(1, wpim1);
				kmfit -> AddTrack(2, wpip2);
				kmfit -> AddTrack(3, wpim2);
				kmfit -> AddTrack(4, wpip3);
				kmfit -> AddTrack(5, wpim3);
				kmfit -> AddTrack(6, 0.0, g1Trk);
				kmfit -> AddTrack(7, 0.0, g2Trk);
				kmfit -> AddTrack(8, 0.0, g3Trk);  //
				kmfit -> AddFourMomentum(0, psip);

				bool oksq = kmfit -> Fit();
				if(oksq){
					double chi2=kmfit->chisq();
					if (chi2 < chisq_t3g)
					{
						chisq_t3g = chi2;
					}
				}
			}
		}
	}
	m_chi23g=chisq_t3g;
	if(debug)cout<<__LINE__<<endl;

	for(int i=0;i<nGam-3;i++)
	{
		for(int j=i+1;j<nGam-2;j++)
		{
			for(int k=j+1;k<nGam-1;k++)
			{
				for(int m=k+1;m<nGam;m++)
				{
					RecEmcShower *g1Trk = (*(evtRecTrkCol -> begin() +iGam[i])) -> emcShower();
					RecEmcShower *g2Trk = (*(evtRecTrkCol -> begin() +iGam[j])) -> emcShower();
					RecEmcShower *g3Trk = (*(evtRecTrkCol -> begin() +iGam[k])) -> emcShower();
					RecEmcShower *g4Trk = (*(evtRecTrkCol -> begin() +iGam[m])) -> emcShower();
					kmfit -> init();
					kmfit -> AddTrack(0, wpip1);
					kmfit -> AddTrack(1, wpim1);
					kmfit -> AddTrack(2, wpip2);
					kmfit -> AddTrack(3, wpim2);
					kmfit -> AddTrack(4, wpip3);
					kmfit -> AddTrack(5, wpim3);
					kmfit -> AddTrack(6, 0.0, g1Trk);
					kmfit -> AddTrack(7, 0.0, g2Trk);
					kmfit -> AddTrack(8, 0.0, g3Trk);  //
					kmfit -> AddTrack(9, 0.0, g4Trk);  //
					kmfit -> AddFourMomentum(0, psip);

					bool oksq = kmfit -> Fit();
					if(oksq){
						double chi2 = kmfit->chisq();
						if (chi2 < chisq_t4g)
						{
							chisq_t4g = chi2;
						}
					}
				}
			}
		}
	}
	m_chi24g=chisq_t4g;
	if(debug)cout<<__LINE__<<endl;

	if(debug)cout<<__LINE__<<endl;

	int Gam[3];
	Gam[0] = i1_pi0[index_pi01];
	Gam[1] = i2_pi0[index_pi01];
	Gam[2] = iGam[index_gam];

	int good_pi[6] = {ipip[0],ipip[1],ipip[2],ipim[0],ipim[1],ipim[2]};
	//	int good_proton[2] = {ipp[0],ipm[0]};

	double dang_pi_tmp=9999;
	double dang_proton_tmp=9999;
	for(int i=0; i<3; i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + Gam[i];
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

		for(int j=0; j<6; j++)
		{
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + good_pi[j];
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();
			double angd = extpos.angle(emcpos);
			if(angd < dang_pi_tmp)dang_pi_tmp = angd;
		}
	}
	m_dang_pi = dang_pi_tmp * 180 / (CLHEP::pi);

	m_tuple->write();

	return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Etac2STo6PiP0::finalize() 
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	cout<<"total number:         "<<Ncut0<<endl;
	cout<<"ntot==4, nCharge==0:  "<<Ncut00<<endl;
	cout<<"nChp = 2, nChn = 2:   "<<Ncut1<<endl;
	cout<<"nGam>=4:              "<<Ncut2<<endl;
	cout<<"pass PID:             "<<Ncut3<<endl;
	cout<<"pi0 list:             "<<Ncut4<<endl;
	cout<<"VertexFit:            "<<Ncut5<<endl;
	cout<<"4C:                   "<<Ncut6<<endl;
	cout<<"6C:         "<<Ncut7<<endl;
	return StatusCode::SUCCESS;
}
// *****************************************************************
// ** A macro to create correlated Gaussian-distributed variables **
// *****************************************************************

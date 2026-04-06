#ifndef __BUILDMETATOWERS_H__
#define __BUILDMETATOWERS_H__

//////////////////////////////////////////////////////////////////////////
//									//
//			Build Meta Towers				//
//									//
// 	Does what is says on the tin, take in calorimeters, applies a   //
// 	vertex based correction and combines the towers to make a	//
// 	single "shifted em+hadronic calorimeter" w/ HCAL bins		//
// 	Works with vandy skimmer data and Tower info			//
//									//
//	Authors: 	Ben Kimelman, Skaydi 				//
//	First commit: 	9 March 2026					//
//	This commit: 	10 March 2026					//
//	version: 	v1.0						//
//									//
//	Notes on this commit: First commit				//
//									//
//////////////////////////////////////////////////////////////////////////


//c++ classes

#include <math.h>
#include <vector>
#include <array>
#include <string>

//Vandy classes
#include <vandyclasses/Tower.h>
#include <vadyclasses/EventInfo.h>

//Calo classes 
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfo.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

//G4 objects
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>

//phool 
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/PHNodeOperation.h>

class BuildMetaTowers
{
	public:
		BuildMetaTowers(
				float radius = 1.245 /*average inner hcal radius*/, 
				std::string input="VandyClass" /*Vandy Class, TowerInfo or Array*/
			       ) 
		{
			//assume that if radius given is small it is given in m, else
			if(radius > 50 ) this->R = radius;
			if(radius < 50 ) this->R = radius * 100;
			this->T	= input;
		}
		~BuildMetaTowers(){};
		
		//struct for Array input 
		struct TowerArrayEntry
		{
			double Energy;
			double phi; //center of tower phi 
			double eta; //center of tower eta
		};
		enum CALO
		{
			META 	= 0,
			EMCAL 	= 1,
			IHCAL	= 2,
			OHCAL	= 3
		};
		//Common Methods
		void RunMetaTowerBuilder(float zVTX)
		{
			for( int i = 0; i < 24; i++ ) shiftedetaEdges.at(i)=calculateEtaShift(etaEdges.at(i), zVTX);
			for( int i = 0; i < 24; i++ )
			{
				double eta =10.; 
				if(i < (int) shiftedetaEdges.size() - 1 )
					eta = 0.5*(shiftedEtaEdges.at(i),  shiftedEtaEdges.at(i+1));
				else 
					eta = 0.5*(shiftedEtaEdges.at(i),  1.1);
				for( int j = 0; j < 64; j++)
				{
					double pt =100.; 
					TowerArrayEntry* tower = 
						new TowerArrayEntry* {0., shiftedEtaEdges.at(i), phiEdges.at(i)};
					if(i < (int) phiEdges.size() - 1 ) 
						phi = 0.5*(phiEdges.at(i),  phiEdges.at(i+1));
					else 	
						phi = 0.5*(phiEdges.at(i),  phiEdges.at(0));
					if(phi > M_PI) 	phi += -2*M_PI;
					if(phi < -M_PI) phi += 	2*M_PI;
					tower->eta = eta;
					tower->phi = phi;
					int index  = CalculateIndex(tower);
					MetaTower->at(index)=tower;
				}

			}
			for(int i=0; i<(int)EMReTowers->size(); i++) AddMetaTower(EMReTowers->at(i));
			for(int i=0; i<(int)IHCaTowers->size(); i++) AddMetaTower(IHCaTowers->at(i));
			for(int i=0; i<(int)OHCaTowers->size(); i++) AddMetaTower(OHCaTowers->at(i));
			return;

		}
		//Fun4All load in 
		void GetFun4AllTowers( 
				TowerInfoContainerv2 CAL, 
				RawTowerGeomContainer_Cylinderv1 geom, 
				std::array<TowerArrayEntry*, 1536>* Towers
				double zVTX,
				CALO calo
				)
		{
			for(int t = 0; t < (int)CAL.size(); t++)
			{
				auto key = CAL.encode_key(n) ;
				auto f4tower 	= CAL.get_tower_at_channel(t);
				int phibin 	= CAL.getTowerPhiBin(key);
				int etabin 	= CAL.getTowerEtabin(key);
				double phi	= geom.get_phicenter(phibin);
				double eta	= geom.get_etaecenter(etabin);
				double E	= f4tower->get_energy();
				eta 		= CalculateEtaShift(eta, zVTX, calo);
				TowerArrayEntry* tower =
					new TowerArrayEntry {E, phi, eta};
				int index 	= CalculateIndex(tower, true);
				Towers->at(index) = tower;
			}
			return;
		}	
		void GetEMCALTowers(
				TowerInfoContainerv2 EMCAL, 
				RawTowerGeomContainer_Cylinderv1 EMCAL_geom, 
				bool isRetower=false,
				double zVTX = 0.
				) //if using Fun4AllTower objects
		{
			if(!isRetower){
				std::cout<<"Please use the retowered EMCAL" <<std::endl;
				return;
			}
			CALO cal = CALO::EMCAL;
			GetFun4AllTowers(EMCAL, EMCAL_geom, EMReTowers, zVTX, cal);
			return;
		}
		void GetHCALTowers(
				TowerInfoContainerv2 HCAL, 
				RawTowerGeomContainer_Cylinderv1 HCAL_geom, 
				bool outer=false,
				double zVTX = 0.
				)
		{
			std::array<TowerEntryArray*, 1536>* HCaTowers = NULL;
			CALO cal = CALO::OHCAL;
			if(outer){
				HCaTowers 	= OHCaTowers;
				cal 		= CALO::OHCAL;
			}
			else{
				HCaTowers 	= IHCaTowers; 	
				cal 		= CALO::IHCAL;
			}
			GetFun4AllTowers(HCAL, HCAL_geom, HCaTowers, zVTX, cal);
			return;
		}

			//if using Fun4AllTower objects
		//Vandy Class load in 
		void GetEMCALTowers(
				std::vector<Tower*> EMCAL,
				bool isRetower=false,
				double zVTX = 0.
				)
		{
			if(!isRetower){
				std::cout<<"Please use retowered"<<std::endl;
				return;
			}
			CALO calo	= CALO::EMCAL;
			GetVandyTowers(EMCAL, EMReTowers, zVTX, calo);
			return;
		}
		void GetHCALTowers(
				std::vector<Tower*> HCAL,
				bool outer=false, 
				double zVTX = 0.
				);
		{
			if(outer)
			{
				CALO calo	= CALO::OHCAL;
				GetVandyTowers(HCAL, OHCaTowers, zVTX, calo);
			}
			else if(!outer)
			{
				CALO calo	= CALO::IHCAL;
				GetVandyTowers(HCAL, IHCaTowers, zVTX, calo);
			}
			return;
		}
		void GetVandyTowers(
				std::vector<Tower*> CAL,
				std::array<TowerArrayEntry*, 1536>* output_array
				double zVTX,
				CALO calo 
				)
		{
			for(int i = 0; i<(int)CAL.size(); i++)
			{
				Tower* tw 	= CAL.at(i);
				double E	= tw->e();
				double px 	= tw->px();
				double py	= tw->py();
				double pz 	= tw->pz();
				double p	= std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
				double phi 	= std::atan2(px, py);
				double eta	= std::atanh(pz/p);
				eta 		= CalculateEtaShift(eta, zVTX, calo);
				TowerArrayEntry* tower = new TowerArrayEntry {E, phi, eta};
				int index	= CalculateIndex(tower, true);
				output_array->at(i) = tower;
			}
			return;
		}	

		//Just basic array 
		void GetEMCALTowers(
				std::vector<TowerArrayEntry*> EMCAL,
			       	bool isRetower = false,
				double zVTX = 0.
				)
		{
			CALO calo = CALO::EMCAL;
			for(int i = 0; i <(int)EMCAL->size(); i++) 
			{
				EMReTowers->at(i) 	= EMCAL->at(i);
				EMReTowers->at(i)->eta 	= calculateEtaShift(EMReTower->at(i)->eta, zVTX, calo);
			}
			return;
		}

		void GetHCALTowers(
				std::vector<TowerArrayEntry*> HCAL,
				bool outer=false,
				double zVTX = 0.
				)
		{
			CALO calo = CALO::OHCAL;
			std::array<TowerArrayEntry*>* HC = OHCaTowers;
			if(!outer){
			       	calo 	= CALO::IHCAL;
				HC	= IHCaTowers;
			for(int i = 0; i <(int)HCAL->size(); i++) 
			{
				HC->at(i)->eta 	= calculateEtaShift(HCAL->at(i)->eta, zVTX, calo);
			}
			return;
		}
		//just returns configs
		void setRadius	( float r = 1.245 ) { this->R = r ; } ;
	       	float getRadius	( ) { return this->R; };
		
		void 		setInput
				( std::string input = "VandyClass" ) { this->T = input; }; 
		std::string 	getInput() { return this->T; }; 
		
		//return the final object
		std::array<TowerArrayEntry*, 1536>* getMetaTowers() { return this->MetaTowers; };

	private:
		//Private variables
		float 		R { 1245 }; //IHCAL half radius in cm
		std::string 	T { "VandyClass" }; //What input Type to expect
		const float 	R_OHCAL_Outer {269};
		const float 	R_OHCAL_Inner {182};
		const float 	R_IHCAL_Outer {137};
		const float 	R_IHCAL_Inner {116};
		const float 	R_EMCAL_Outer {116};
		const float	R_EMCAL_Inner { 90};	
		const float	R_OHCAL_MID { 1/2.*(R_OHCAL_Outer + R_OHCAL_Inner) };	
		const float	R_IHCAL_MID { 1/2.*(R_IHCAL_Outer + R_IHCAL_Inner) };	
		const float	R_EMCAL_MID { 1/2.*(R_EMCAL_Outer + R_EMCAL_Inner) };	
		std::array <TowerArrayEntry*, 1536>* MetaTowers = new std::array <TowerArrayEntry*, 1536> {}; 
		std::array <TowerArrayEntry*, 1536>* EMReTowers = new std::array <TowerArrayEntry*, 1536> {}; 
		std::array <TowerArrayEntry*, 1536>* IHCaTowers = new std::array <TowerArrayEntry*, 1536> {}; 
		std::array <TowerArrayEntry*, 1536>* OHCaTowers = new std::array <TowerArrayEntry*, 1536> {};
		
		//HCAL eta-phi physical geometry
		const std::array <double, 24> etaEdges 	
			{ -1.1000000, -1.0083333, -0.91666667, -0.82500000, -0.73333333, -0.64166667, 
				-0.55000000, -0.45833333, -0.36666667, -0.27500000, -0.18333333, -0.091666667, 
				1.1102230e-16, 0.091666667, 0.18333333, 0.27500000, 0.36666667, 0.45833333, 
				0.55000000, 0.64166667, 0.73333333, 0.82500000, 0.91666667, 1.0083333, 1.1000000 };

		const std::array <double, 64> phiEdges
			{ -0.053619781, 0.044554989, 0.14272976, 0.24090453, 0.33907930, 0.43725407, 
				0.53542884, 0.63360361, 0.73177838, 0.82995315, 0.92812792, 1.0263027, 
				1.1244775, 1.2226522, 1.3208270, 1.4190018, 1.5171765, 1.6153513, 1.7135261, 
				1.8117009, 1.9098756, 2.0080504, 2.1062252, 2.2043999, 2.3025747, 2.4007495, 
				2.4989242, 2.5970990, 2.6952738, 2.7934486, 2.8916233, 2.9897981, 3.0879729, 
				3.1861476, 3.2843224, 3.3824972, 3.4806720, 3.5788467, 3.6770215, 3.7751963, 
				3.8733710, 3.9715458, 4.0697206, 4.1678953, 4.2660701, 4.3642449, 4.4624197, 
				4.5605944, 4.6587692, 4.7569440, 4.8551187, 4.9532935, 5.0514683, 5.1496431, 
				5.2478178, 5.3459926, 5.4441674, 5.5423421, 5.6405169, 5.7386917, 5.8368664, 
				5.9350412, 6.0332160, 6.1313908, 6.2295655 };	
		
		std::array<double ,24> shiftedetaEdges {};

		//Private Methods
		double calculateEtaShift(double eta, double zVtx)
		{
			double z 	= R*std::sinh(eta);
			double zshift	= z-zVtx;
			zshift		= zshift / R;
			double etashift	= std::asinh(zshift);
			return etashift;
		}
		double calculateEtaShift(double eta, double zVtx, CAlO calo)
		{
			double r = this->R;
			if( 		calo == CALO::EMCAL ) 	r = R_EMCAL_MID;
			else if( 	calo == CALO::IHCAL )	r = R_IHCAL_MID;
			else if(	calo == CALO::OHCAL ) 	r = R_OHCAL_MID;
			else 					r = this->R;
			double z 	= R*std::sinh(eta);
			double zshift	= z-zVtx;
			zshift		= zshift / R;
			double etashift	= std::asinh(zshift);
			return etashift;
		}

			
		int calculateIndex(TowerArrayEntry tower, bool useSlant) 
		{
			int etabin 	= 0;
			int phibin	= 0;
			etabin 		= LookupBin(tower.eta, useSlant, true);
			phibin		= LookupBin(tower.phi, useSlant, false);
			int index = etabin * 64 + phibin ;
			return index;
		}
		int LookupBin(double towervalue, bool useSlant, bool isEta)
		{
			std::vector<double> tempBottom;
			int binN=0;
			if(isEta && useSlant)
			{
				for(int i = 0; i<(int) shiftedEtaEdges.size(); i++){
					tempBottom.push_back(shiftedEtaEdges.at(i));
				}
			}
			else if (isEta)
			{
				for(int i = 0; i<(int) EtaEdges.size(); i++){
					tempBottom.push_back(EtaEdges.at(i));
				}
			}
			else 
			{
				for(int i = 0; i<(int) PhiEdges.size(); i++){
					tempBottom.push_back(PhiEdges.at(i));
				}
			}
			for(int i = 0; i<(int)tempBottom.size(); i++)
			{
				if(towervalue >= tempBottom.at(i)){
					binN = i;
					continue;
				}
				if(towervalue < tempBottom.at(i)){
					break;
				}
			}
			return binN;
		}
		void decodeIndex(int index, int* ebr, int* pbr) 
		{
			
			int etabin	= index / 64;
			int phibin	= index % 64;
			ebr = &etabin;
			pbr = &phibin;
			return;
		}
		
		void shiftBinMins(double zVtx)
		{
			for(int i= 0; i<(int) etaEdges.size(); i++)
			{
				shiftedetaEdges[i]=calculateEtaShift(etaEdges[i], zVtx);
			}
			return;
		}
		void shiftEMCAL(std::vector<TowerArrayEntry*> emcal, bool isRetower, float zVTX)
		{
			if(!isRetower || (int) emcal.size() > 1536)
			{
				std::vector<TowerArrayEntry*> rtEM = RetowerEMCAL(emcal);
				emcal.clear();
				emcal = rtEM;
			}
			for(int i = 0; i<(int) emcal.size(); i++)
			{
				TowerArrayEntry* nTower = shiftTower(emcal.at(i), CALO::EMC, zVTX);
				int N = calculateIndex( *nTower );
				EMReTowers->at(N) = nTower;
			}
			return; 
		}
	
		void shiftIHCAL(std::vector<TowerArrayEntry*> ihcal, float zVTX);
		{
			for(int i = 0; i<(int) ihcal.size(); i++)
			{
				TowerArrayEntry* nTower = shiftTower(ihcal.at(i), CALO::IHC, zVTX);
				int N = calculateIndex( *nTower );
				IHCaTowers->at(N) = nTower;
			}
			return; 
		}
		void shiftOHCAL(std::vector<TowerArrayEntry*> ohcal, float zVTX);
		{
			for(int i = 0; i<(int) ohcal.size(); i++)
			{
				TowerArrayEntry* nTower = shiftTower(ohcal.at(i), CALO::OHC, zVTX);
				int N = calculateIndex( *nTower );
				OHCaTowers->at(N) = nTower;
			}
			return; 
		}
		void addMetaTower(TowerArrayEntry tower)
		{
			double eta_val 	= tower.eta;
			double phi_val	= tower.phi;
			int etaBin 	= findEtaBin(eta_val);
			int phiBin	= findPhiBin(phi_val);
			int N 		= calculateIndex(tower);
			if (!MetaTowers->at(N)) MetaTowers->at(N) = 
				new TowerArrayEntry* { 0, getPhiCenter(phiBin), getEtaCenter(etaBin) };
			MetaTowers->at(N)->E += tower.E;
			return;
		}
		TowerArrayEntry* shiftTower(TowerArrayEntry* tower, CALO calo, float zVTX) 
		{
			
			tower->eta = CalculateEtaShift(tower->eta, zVTX, calo);
			return tower; 
		}

};
			
#endif

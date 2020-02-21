#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"


class FilterAK8Jet : public edm::stream::EDFilter<> {
	public:
		explicit FilterAK8Jet(const edm::ParameterSet&);


	private:
		virtual bool filter(edm::Event&, const edm::EventSetup&) override;

		edm::Handle<pat::JetCollection> jetsAK8;
		edm::EDGetTokenT<pat::JetCollection> jetsAK8Token;
		edm::Handle<pat::JetCollection> jetsAK4;
		edm::EDGetTokenT<pat::JetCollection> jetsAK4Token;
};


FilterAK8Jet::FilterAK8Jet(const edm::ParameterSet& iConfig):
	jetsAK8Token (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8")))
	, jetsAK4Token (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK4")))
{
}

bool FilterAK8Jet::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	using namespace edm;
	using namespace std;

	iEvent.getByToken(jetsAK8Token, jetsAK8);
        iEvent.getByToken(jetsAK4Token, jetsAK4);

	size_t j = 0;
	for(size_t i=0; i<jetsAK8->size(); i++) {
		if(jetsAK8->at(i).pt()>170 && abs(jetsAK8->at(i).eta())<2.4) j++;
	}

        size_t j1 = 0;
        size_t j2 = 0;
	double NHF, NEMF, NumConst, CHF, CHM, CEMF;
	bool tightJetID;
	for(size_t i = 0;i<jetsAK4->size();i++){
		if(jetsAK4->at(i).pt()>30 && abs(jetsAK4->at(i).eta())<2.4){
			NHF  = jetsAK4->at(i).neutralHadronEnergyFraction();
			NEMF = jetsAK4->at(i).neutralEmEnergyFraction();
			CHF  = jetsAK4->at(i).chargedHadronEnergyFraction();
			CEMF = jetsAK4->at(i).chargedEmEnergyFraction();
			NumConst = jetsAK4->at(i).chargedMultiplicity()+jetsAK4->at(i).neutralMultiplicity();
			CHM      = jetsAK4->at(i).chargedMultiplicity();
			tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(jetsAK4->at(i).eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jetsAK4->at(i).eta())>2.4) && abs(jetsAK4->at(i).eta())<=2.7;
			if(tightJetID && jetsAK4->at(i).pt()>40) j1++;
                        if(tightJetID && jetsAK4->at(i).pt()>30) j2++;
		}
	}


	if (j>1 || (j1>0 && j2>1)) {return true;} else return false;

}
//define this as a plug-in
DEFINE_FWK_MODULE(FilterAK8Jet);

#include "G4EmDNAPhysics_option2and4.hh"
#include "G4PhysicsConstructorFactory.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4EmParameters.hh"
#include "G4VEmModel.hh"
//end of warning
#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4PhysicsListHelper.hh"
//particles
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

//process and model
#include "G4DNAElectronSolvation.hh"
#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAELSEPAElasticModel.hh"
#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"

// Warning : the following is needed in order to use EM Physics builders
// e+
#include "G4Positron.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
// gamma
#include "G4Gamma.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"
// end of warning

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4DummyModel.hh"
#include "G4EmConfigurator.hh"
//need G4PhysicsConstructorFactory
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option2and4);  //let the physicsConstructor in factory

//Constructor
G4EmDNAPhysics_option2and4::G4EmDNAPhysics_option2and4(G4int ver) : G4VPhysicsConstructor("G4EmDNAPhysics_option2and4"), verbose(ver)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetFluo(true);  
    param->SetAuger(true);  
    param->SetAugerCascade(true);  
    param->SetDeexcitationIgnoreCut(true);
    param->ActivateDNA();
    param->SetVerbose(verbose);
    //enum sub-type of electro-magetic in G4BUilderType.hh([0,bUnknown],[1, bTransportation], [2, bElectromagnetic], [3, bEmExtra],[4, bDecay], [5, bHadronElastic], [6, bHadronInelastic], [7,bStopping], [8, bIons])
    SetPhysicsType(bElectromagnetic);

}

//Destructor
G4EmDNAPhysics_option2and4::~G4EmDNAPhysics_option2and4()
{

}

// Construct Particles 
void G4EmDNAPhysics_option2and4::ConstructParticle()
{
    //bosons
    G4Gamma::Gamma();

    // leptons
    G4Electron::Electron();
    G4Positron::Positron();
  
    // baryons
    G4Proton::Proton();

    G4GenericIon::GenericIonDefinition();

    G4DNAGenericIonsManager * genericIonsManager;
    genericIonsManager=G4DNAGenericIonsManager::Instance();
    genericIonsManager->GetIon("alpha++");
    genericIonsManager->GetIon("alpha+");
    genericIonsManager->GetIon("helium");
    genericIonsManager->GetIon("hydrogen");

}

void G4EmDNAPhysics_option2and4::ConstructProcess()
{
    if(verbose > 1) 
    {
        G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
    }
    /*********************************************************************
    *               Iterate over all particle species                    *
    *********************************************************************/
    G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();
    auto myParticleIterator = GetParticleIterator();  // Get the particle iterator
    myParticleIterator->reset();
    while((*myParticleIterator)())
    {
        G4ParticleDefinition *particle = myParticleIterator->value(); // Get the particle Definition
        G4String particleName = particle->GetParticleName();   // Get the particle Name;particle->GetName() is obsolete
        /**********************************************************************************
        *                                 e-                                              *
        *         Elastic  ionisation excitation vibratexcitation attachment solvation 
        ***********************************************************************************/
        if(particleName == "e-")
        {
        /******************************************************************************************************
        //                               Solvation
        //         subexcitation electrons (<10 eV) occur, reach to 25meV
        //  G4DNAModelSubType(0,1,2,3,4,5):Unknown, fRichie1994eSolvation,fTerrisol1990eSolvation,
        //   fMeesungnoen2002eSolvation, fKreipl2009eSolvation,fMeesungnoenSolid2002eSolvation 
        *********************************************************************************************************/
            G4DNAElectronSolvation *solvation = new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
            // use macro file control as: /process/dna/e-SolvationSubType Meesungnoen2002
            //candidate: Ritchie1994, Terrisol1990,Kreipl2009,Meesungnoen2002_amorphous,Meesungnoen2002
            auto therm = G4DNASolvationModelFactory::GetMacroDefinedModel();//G4DNAOneStepThermalizationModel* therm = new G4DNAOneStepThermalizationModel(); 
            therm->SetHighEnergyLimit(10.*eV);// the limit of Champion Elastic model:7.4 eV , ELSEPA model:10eV,
            solvation->SetEmModel(therm);
            ph->RegisterProcess(solvation, particle);
        
        /*********************************************************************************************************
        //                              Elastic scattering
        // at low energy is dominant
        // Model: Champion, ELSEPA, UeharaScreenedRutherford, CPA100
        // will affect the distribution  of particles
        **********************************************************************************************************/
            G4DNAElastic *theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
            theDNAElasticProcess->SetEmModel(new G4DummyModel(), 1);   // use dummy model,will be modified
            ph->RegisterProcess(theDNAElasticProcess, particle);
        //Have two methods to attach A model to B model: this, and option8

        // Excitation
            G4DNAExcitation* theDNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");
            theDNAExcitationProcess->SetEmModel(new G4DummyModel(), 1);
            ph->RegisterProcess(theDNAExcitationProcess, particle);

            // Ionisation
            G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");
            theDNAIonisationProcess->SetEmModel(new G4DNAEmfietzoglouIonisationModel());
            G4DNABornIonisationModel* higherEnergyIonisation = new G4DNABornIonisationModel(); //initiate Born Ionisation model
            higherEnergyIonisation->SetLowEnergyLimit(10.*keV); // set the min limit to 10 keV ( the maximum limit of Emfietzoglou)
            theDNAIonisationProcess->SetEmModel(higherEnergyIonisation);
            ph->RegisterProcess(theDNAIonisationProcess, particle);

            //use the default Born Vibexcitation and attachment model
            // Vibrational excitation
            ph->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), particle);
            
            // Attachment
            ph->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), particle); 
        }
        /**************************************************************************************
        *                                 Protons                                              *
        ****************************************************************************************/
        else if(particleName == "proton")
        {
            ph->RegisterProcess(new G4DNAElastic("proton_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("proton_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);
       /*
            //  In option2 
            G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("proton_G4DNAIonisation");

            G4VEmModel* mod1 = new G4DNARuddIonisationExtendedModel();
            mod1->SetLowEnergyLimit(0*eV);
            mod1->SetHighEnergyLimit(500*keV);

            G4DNABornIonisationModel* mod2 = new G4DNABornIonisationModel();
            mod2->SetLowEnergyLimit(500*keV);
            mod2->SetHighEnergyLimit(100*MeV);
            mod2->SelectFasterComputation(true);

            theDNAIonisationProcess->SetEmModel(mod1);
            theDNAIonisationProcess->SetEmModel(mod2);

            ph->RegisterProcess(theDNAIonisationProcess, particle);
	    
            ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);
       */
        }
        /****************************************************************************************
        *                                 Hydrogen                                              *
        ****************************************************************************************/
        else if (particleName == "hydrogen") 
        {
            ph->RegisterProcess(new G4DNAElastic("hydrogen_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("hydrogen_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), particle);
        }
        /****************************************************************************************
        *                                 alpha                                                 *
        ****************************************************************************************/
        else if(particleName == "alpha")
        {
            ph->RegisterProcess(new G4DNAElastic("alpha_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("alpha_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), particle);
        }
        /****************************************************************************************
        *                                 alpha+                                                *
        ****************************************************************************************/
        else if(particleName == "alpha+")
        {
            ph->RegisterProcess(new G4DNAElastic("alpha+_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("alpha+_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), particle);
        }
        /****************************************************************************************
        *                                 helium                                                *
        ****************************************************************************************/
        else if(particleName == "helium")
        {
            ph->RegisterProcess(new G4DNAElastic("helium_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("helium_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), particle);
        }
        /****************************************************************************************
        *                                 GenericIon by Z. Francis(HZE)                         *
        ****************************************************************************************/
        else if(particleName == "GenericIon")
        {
            ph->RegisterProcess(new G4DNAIonisation("GenericIon_G4DNAIonisation"), particle);
        }
        //==========================================================================================
        // Warning : the following particles and processes are needed by EM Physics builders. They
        // are taken from the default Livermore Physics list. These particles are currently not
        // handled by Geant4-DNA.
        //==========================================================================================

        /****************************************************************************************
        *                                Positrons(e+)                                          *
        ****************************************************************************************/
        else if(particleName == "e+")
        {
            // the same as G4EmStandardPhysics_option3
            G4eMultipleScattering* msc = new G4eMultipleScattering();
            msc->SetStepLimitType(fUseDistanceToBoundary);
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.2, 100*um);      

            ph->RegisterProcess(msc, particle);
            ph->RegisterProcess(eIoni, particle);
            ph->RegisterProcess(new G4eBremsstrahlung(), particle);
            ph->RegisterProcess(new G4eplusAnnihilation(), particle);
        }
        /****************************************************************************************
        *                               Photons(gamma)                                          *
        ****************************************************************************************/
        else if(particleName == "gamma")
        {
             // photoelectric effect - Livermore model only
            G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
            thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
            ph->RegisterProcess(thePhotoElectricEffect, particle);

            // Compton scattering - Livermore model only
            G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
            theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
            ph->RegisterProcess(theComptonScattering, particle);

            // gamma conversion - Livermore model below 80 GeV
            G4GammaConversion* theGammaConversion = new G4GammaConversion();
            theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
            ph->RegisterProcess(theGammaConversion, particle);

            // default Rayleigh scattering is Livermore
            G4RayleighScattering* theRayleigh = new G4RayleighScattering();
            ph->RegisterProcess(theRayleigh, particle);}
            // Warning : end of particles and processes are needed by EM Physics builders 
        }
        /****************************************************************************************
        *                           Optionally handle deexcitation process                      *
        ****************************************************************************************/
        G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
        G4LossTableManager::Instance()->SetAtomDeexcitation(de);
        /****************************************************************************************
        *                 First methods:Set different models in different energies              *
        ****************************************************************************************/
        // use the G4EmConfigurator EM models for regions
        G4EmConfigurator *em_config = G4LossTableManager::Instance()->EmConfigurator();
        G4VEmModel* mod;
        /***************************************************************************************
        *   Elastic scattering. Use Champion(> 10keV) and Uehara(<10keV)
        ****************************************************************************************/
        mod = new G4DNAChampionElasticModel();
        mod->SetActivationLowEnergyLimit(10.*keV);
        em_config->SetExtraEmModel(
                                    "e-",              // Particle name
                                    "e-_G4DNAElastic", // Process name
                                    mod,               // G4VEmModel pointer
                                    "World",          // Region name // in the whole region
                                    10.*keV,           // Lower energy limit
                                    1.*MeV             // Upper energy limit
                                  );

        // Use Uehara below 10 keV
        mod = new G4DNAUeharaScreenedRutherfordElasticModel();
        em_config->SetExtraEmModel(
                                    "e-",              // Particle name
                                    "e-_G4DNAElastic", // Process name
                                    mod,               // G4VEmModel pointer
                                    "World",          // Region name
                                    0.,                // Lower energy limit
                                    10.*keV            // Upper energy limit
                                  );                                 
        /***************************************************************************************
        *   Excitation and Ionisation. Use Born(> 10keV) and Emfietzoglou (<10keV)
        ****************************************************************************************/
        mod = new G4DNABornExcitationModel();
        mod->SetActivationLowEnergyLimit(10.*keV);
        em_config->SetExtraEmModel(
                                    "e-",
                                    "e-_G4DNAExcitation",
                                    mod,
                                    "World",
                                    10.*keV,
                                    1.*MeV
                                   );

        mod = new G4DNAEmfietzoglouExcitationModel();
        em_config->SetExtraEmModel(
                                    "e-",
                                    "e-_G4DNAExcitation",
                                    mod,
                                    "World",
                                     0.,
                                    10.*keV
                                  );
    //-----------------------------------------------------------------------------------
    	mod = new G4DNABornIonisationModel();
        mod->SetActivationLowEnergyLimit(10.*keV);
        em_config->SetExtraEmModel(
                                    "e-",
                                    "e-_G4DNAIonisation",
                                    mod,
                                    "World",
                                    10.*keV,
                                    1.*MeV
                                  );

    mod = new G4DNAEmfietzoglouIonisationModel();
    em_config->SetExtraEmModel(
                                "e-",
                                "e-_G4DNAIonisation",
                                mod,
                                "World",
                                0.,
                                10.*keV
                              );

}

#ifndef G4EMDNAPHYSICS_OPTION2AND4_HH
#define G4EMDNAPHYSICS_OPTION2AND4_HH

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4EmDNAPhysics_option2and4 : public G4VPhysicsConstructor
{
    public:
        G4EmDNAPhysics_option2and4(G4int ver = 0);
        virtual ~G4EmDNAPhysics_option2and4();
        virtual void ConstructParticle();
        virtual void ConstructProcess();
    private:
        G4int verbose;

};
#endif

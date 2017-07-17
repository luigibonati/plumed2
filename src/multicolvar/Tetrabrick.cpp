/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR TETRAHEDRAL ORDER PARAMETER 
/*
Work in progress...

*/
//+ENDPLUMEDOC


class Tetrabrick : public MultiColvarBase {
private:
//  double nl_cut;
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit Tetrabrick(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(Tetrabrick,"TETRABRICK")

void Tetrabrick::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); 
}

Tetrabrick::Tetrabrick(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     switchingFunction.set(sw,errors);
     if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else { 
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  coordination of central atom and those within %s\n",( switchingFunction.description() ).c_str() );
  // Set the link cell cutoff
  setLinkCellCutoff( switchingFunction.get_dmax() );
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();
}

double Tetrabrick::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   // --- Calculate the tetrahedral order parameter ---
   
   // Define output quantities
   double tetra=0;
   vector<Vector> deriv(getNumberOfAtoms());
   Tensor virial;
   // Define temp quantities
   double d2i, d2j;								//square distances
   double cos; 									//cosine term
   double sw_i, sw_j, df_i, df_j;				//switch functions and derivatives
   
   // Loop on nearest neighbors and load distances di and dj
   for(unsigned i=1;i<(myatoms.getNumberOfAtoms()-1);++i){
      Vector& di=myatoms.getPosition(i); 		//relative position of atom i (with respect to k)
      if ( (d2i=di[0]*di[0])<rcut2 && 
           (d2i+=di[1]*di[1])<rcut2 &&
           (d2i+=di[2]*di[2])<rcut2) {
         for(unsigned j=i+1;j<myatoms.getNumberOfAtoms();++j){
            Vector& dj=myatoms.getPosition(j); 	//relative position of atom j (with respect to k)  
            if ( (d2j=dj[0]*dj[0])<rcut2 && 
                 (d2j+=dj[1]*dj[1])<rcut2 &&
                 (d2j+=dj[2]*dj[2])<rcut2) {
			// compute cosine term	
			cos = (di[0]*dj[0]+di[1]*dj[1]+di[2]*dj[2])/(sqrt(d2i)*sqrt(d2j)) +1./3.;
			// compute switching functions
			sw_i = switchingFunction.calculateSqr( d2i, df_i );
			sw_j = switchingFunction.calculateSqr( d2j, df_j );

            tetra += cos*cos*sw_i*sw_j;
			}
         }
      }
   }
   
   tetra=1 - 3./8.*tetra; 
   cout << tetra << endl;
   return tetra;
}

}
}


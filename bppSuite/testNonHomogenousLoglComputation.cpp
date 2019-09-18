/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encoutaged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Phyl/Model/Codon/YNGP_M2.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <iostream>


using namespace bpp;
using namespace std;


void printModelParameters(TreeLikelihood* tl)
{
  ParameterList parameters = tl->getParameters();
  for (size_t i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  cout << "\n" << endl;
}


int main() 
{
    try
    {
        // process tree
        TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.01):0.02,C:0.03,D:0.04);");
        Tree* ttree = dynamic_cast<Tree*>(tree);

        // process sequence data
        map<string, string> alphabetParams;
        alphabetParams["alphabet"] = "Codon(letter=DNA)";
        alphabetParams["genetic_code"] = "Standard";
        const Alphabet* nucAlphabet = SequenceApplicationTools::getAlphabet(alphabetParams, "", false);
        const CodonAlphabet* alphabet = dynamic_cast<const CodonAlphabet*>(nucAlphabet);
        unique_ptr<GeneticCode> gCode;
        gCode.reset(SequenceApplicationTools::getGeneticCode(alphabet->getNucleicAlphabet(), "Standard"));
        VectorSiteContainer sites(alphabet);
        sites.addSequence(BasicSequence("A", "AAATGGCTGTGCACGTCT", alphabet));
        sites.addSequence(BasicSequence("B", "AACTGGATCTGCATGTCT", alphabet));
        sites.addSequence(BasicSequence("C", "ATCTGGACGTGCACGTGT", alphabet));
        sites.addSequence(BasicSequence("D", "CAACGGGAGTGCGCCTAT", alphabet));

        // set partition A and feed it to the RELAX model with k=1
        map<string,string> params;
        params["model1"] = "YNGP_M2(kappa=2.0,omega0=0.1,omega2=2.0,theta1=0.5,theta2=0.8,frequencies=F0)";
        params["model2"] = "YNGP_M2(kappa=YNGP_M2.kappa_1,omega0=0.01,omega2=4,theta1=YNGP_M2.theta1_1,theta2=YNGP_M2.theta2_1,frequencies=F0)";
        params["nonhomogeneous"]="general";
        params["nonhomogeneous.number_of_models"] = "2";
        params["nonhomogeneous.stationarity"] = "yes";
        params["site.number_of_paths"] = "2";                               // the 3rd path mapping omega3 in the branches under chatacter states 0 and 1 is imlies the the other two paths
        params["site.path1"] = "model1[YN98.omega_1]&model2[YN98.omega_1]"; // map omega1 in the branches under character state 0 (=model1) to omega1 in the branches under character state 1 (=model2) 
        params["site.path2"] = "model1[YN98.omega_2]&model2[YN98.omega_2]";
        params["model1.nodes_id"] = "0";
        params["model2.nodes_id"] = "1,2,3,4";
        ConstantRateDistribution* rdist = new ConstantRateDistribution();
        MixedSubstitutionModelSet* DoubleM2Model = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), dynamic_cast<const SiteContainer*>(&sites), params));
        RNonHomogeneousMixedTreeLikelihood* DoubleM2TreeLikelihood = new RNonHomogeneousMixedTreeLikelihood(*ttree, dynamic_cast<const SiteContainer&>(sites), DoubleM2Model, rdist, true, false);
        DoubleM2TreeLikelihood->initialize();
        double DoubleM2LogLikelihood = -1*DoubleM2TreeLikelihood->getValue(); 

        // make sure you receive the same likelihood as YNGP_M2 (simple site model)
        map<string,string> m2params;
        m2params["model"] = "YNGP_M2(kappa=2.0,omega0=0.1,omega2=2.0,theta1=0.5,theta2=0.8,frequencies=F0)";
        m2params["nonhomogeneous"] = "no";
        TransitionModel* M2Model = PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), dynamic_cast<const SiteContainer*>(&sites), m2params);
        RHomogeneousMixedTreeLikelihood* M2TreeLikelihood = new RHomogeneousMixedTreeLikelihood(*ttree, dynamic_cast<const SiteContainer&>(sites), M2Model, rdist, true, false);
        M2TreeLikelihood->initialize();
        double M2LogLikelihood = -1*M2TreeLikelihood->getValue();
        if (abs(DoubleM2LogLikelihood - M2LogLikelihood) > 0.0001) // 2.9.19 bug here - caused by different transition probabilities, but was also occuring in the former version - need to check that it happens in M3 and report to Laurent
        {
            cout << "Error! Non homogenous M2 with paths restrcitions yields different likelihood than homogenous M2 model" << endl;
			cout << "NonHomogenous M2 Log Likelihood: " << DoubleM2LogLikelihood << endl;
            vector<double> DoubleM2LoglBySite = DoubleM2TreeLikelihood->getLikelihoodForEachSite();
            double DoubleM2TransitionRate_1 = dynamic_cast<MixedSubstitutionModel*>(DoubleM2Model->getModel(0))->Qij(0,0);
            double DoubleM2TransitionRate_2 = dynamic_cast<MixedSubstitutionModel*>(DoubleM2Model->getModel(1))->Qij(0,0);
            double DoubleM2TransitionProb = DoubleM2TreeLikelihood->getTransitionProbabilitiesPerRateClass(0,0)[0][0][0];
            double DoubleM2Pij_1 = (DoubleM2Model->getModel(0))->getPij_t(0)(0,0);
            double DoubleM2Pij_2 = (DoubleM2Model->getModel(1))->getPij_t(0)(0,0);
			printModelParameters(DoubleM2TreeLikelihood);
			cout << "Honogenous M2 Log Likelihood: " << M2LogLikelihood << endl;
            vector<double> M2LoglBYSite = M2TreeLikelihood->getLikelihoodForEachSite();
            double M2TransitionRate = dynamic_cast<MixedSubstitutionModel*>(M2Model)->Qij(0,0);
            double M2TransitionProb = M2TreeLikelihood->getTransitionProbabilitiesPerRateClass(0,0)[0][0][0];
            double M2Pij = M2Model->getPij_t(0)(0,0);
			printModelParameters(M2TreeLikelihood);
            return 1;
        }
        
        // free resources
        delete rdist;
        delete M2Model;
        delete M2TreeLikelihood;
        delete DoubleM2Model;
        delete DoubleM2TreeLikelihood;
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}